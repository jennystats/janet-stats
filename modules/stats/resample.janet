# resample.janet - bootstrap (percentile + BCa) and permutation
# tests. R conventions decoded from the boot 1.3-32 source
# (boot.ci / bca.ci / perc.ci / norm.inter / usual.jack).
(import typed :prefix "")
(import ./dist :prefix "")
(import ./num)

(defn arrays-only
  "Composition guard: the statistic loops (user fns), RNG draws and
  order statistics need janet-side values - a packed handle unpacks
  first (actionable). Resolves :impl first."
  [opts & args]
  (num/check-impl-landed opts)
  (each x args
    (unless (indexed? x)
      (error "stats: composition path consumes janet arrays; unpack packed handles first")))
  nil)

(defn/typed mean-plain {:args [:any] :ret :number} [xs]
  "Plain mean (the 1.0-force keeps int data honest)."
  (var s 0)
  (each v xs (+= s v))
  # force double: exact / on int data yields janet rationals, which are
  # NOT `=`-comparable to doubles despite printing identically
  (/ s (* 1.0 (length xs))))

#==------------------------------------------------------------------------==#
# norm.inter - R's quantile convention for bootstrap replicates (NOT
# type-7): integer rk -> exact order statistic; else normal-score
# interpolation between adjacent order statistics.
#==------------------------------------------------------------------------==#

(defn/typed norm-inter {:args [:any :number] :ret :number} [t alpha]
  "R's norm.inter quantile convention for bootstrap replicates (NOT
  type-7): integer rank -> exact order statistic; else normal-score
  interpolation between adjacent order statistics."
  (def sorted (sorted t))
  (def r (length sorted))
  (when (< r 2) (error "stats: norm-inter needs >= 2 replicates"))
  (def rk (* (inc r) alpha))
  (def k (math/floor rk))
  (cond
    (= k 0) (get sorted 0)
    (= k r) (get sorted (dec r))
    (= k rk) (get sorted (dec k))
    (let [z1 (qnorm alpha 0 1 {})
          z2 (qnorm (/ k (inc r)) 0 1 {})
          z3 (qnorm (/ (inc k) (inc r)) 0 1 {})
          tk (get sorted (dec k))
          tk1 (get sorted k)]
      (+ tk (* (/ (- z1 z2) (- z3 z2)) (- tk1 tk))))))

#==------------------------------------------------------------------------==#
# bootstrap - R boot(stype="i") resampling; math/rng is deterministic
# per seed string, so (bootstrap ... :seed "s") replays identically.
#==------------------------------------------------------------------------==#

(defn bootstrap
  "Resample bootstrap (R boot stype='i'): stat is a two-argument
  statistic (fn [data indices]) applied to resampled index vectors.
  Returns {:data :stat :t0 :t :B} - t0 the statistic on all indices,
  t the B replicates. Opts: {:B 2000 :seed \"stats\" :impl}; the
  replicate set is deterministic per :seed. Throws on :B not a
  positive integer or n < 2. Composition (user-fn loop - no kernel
  path)."
  [data stat &opt opts]
  (default opts {})
  # the replicate loop is bound by the
  # USER statistic fn and the RNG - no kernel shape applies; the stat
  # fn itself may be any landed stats fn at its own resolution
  (arrays-only opts data)
  (def b (get opts :B 2000))
  (unless (and (int? b) (pos? b)) (error "stats: :B must be a positive integer"))
  (def n (length data))
  (when (< n 2) (error "stats: bootstrap needs >= 2 observations"))
  (def rng (math/rng (string (get opts :seed "stats"))))
  (def t0 (stat data (seq [i :range [0 n]] i)))
  (def t (array/new b))
  (for _ 0 b
    (def idx (array/new n))
    (for j 0 n (array/push idx (math/rng-int rng n)))
    (array/push t (stat data idx)))
  {:data data :stat stat :t0 t0 :t t :B b})

(defn jackknife
  "Delete-one jackknife: stat (fn [data indices]) evaluated on every
  leave-one-out subset. Returns {:theta-hat :values :L} - L the
  (n-1)(t-hat - t-minus-i) influence values (boot's usual.jack),
  which feed boot-ci's :bca (:L opts key). Same user-fn binding as
  bootstrap."
  [data stat &opt opts]
  (default opts {})
  # same user-fn binding as bootstrap
  (arrays-only opts data)
  (def n (length data))
  (when (< n 2) (error "stats: jackknife needs >= 2 observations"))
  (def full (seq [i :range [0 n]] i))
  (def t-hat (stat data full))
  (def values (array/new n))
  (def lvals (array/new n))
  (for i 0 n
    (def minus (filter |(not= $ i) full))
    (def tm (stat data minus))
    (put values i tm)
    (put lvals i (* (dec n) (- t-hat tm))))
  {:theta-hat t-hat :values values :L lvals})

(defn/typed fin-numbers {:args [:any] :ret :array} [t]
  "The finite, non-NaN replicates of t."
  (filter |(and (number? $) (= $ $) (< (math/abs $) math/inf)) t))

(defn boot-ci
  "Bootstrap confidence interval from a bootstrap record.
  Opts: {:conf 0.95 :type :percentile|:bca :L :impl}; :bca needs
  :L - the jackknife influence values (jackknife's :L). Returns a
  [lo hi] tuple: :percentile by R's norm.inter convention, :bca by
  the R bca.ci adjustment (acceleration from the L sums). Throws on
  fewer than 2 finite replicates or an unsupported :type."
  [boot-out &opt opts]
  (default opts {})
  # order statistics (norm-inter sorts) are composition; the bca
  # influence sums (sl2/sl3) kernel under :num (vdot family)
  (num/check-impl-landed opts)
  (def conf (get opts :conf 0.95))
  (def type-kw (get opts :type :percentile))
  (def t (fin-numbers (boot-out :t)))
  (def t0 (boot-out :t0))
  (def r (length t))
  (when (< r 2) (error "stats: boot-ci needs >= 2 finite replicates"))
  (def alpha [(/ (- 1 conf) 2) (/ (+ 1 conf) 2)])
  (case type-kw
    :percentile
    [(norm-inter t (alpha 0)) (norm-inter t (alpha 1))]
    :bca
    (do
      (def lvals (get opts :L))
      (unless lvals (error "stats: :bca needs :L (jackknife influence values)"))
      (var less 0)
      (each v t (when (< v t0) (++ less)))
      (def w (qnorm (/ less r) 0 1 {}))
      (def sl2
        (if (= :num (num/resolve-impl opts))
          (let [k (num/kernels)
                lv (num/stage lvals k)]
            ((k 'vdot) lv lv))
          (do
            (var s 0)
            (each v lvals (+= s (* v v)))
            s)))
      (def sl3
        (if (= :num (num/resolve-impl opts))
          (let [k (num/kernels)
                lv (num/stage lvals k)]
            ((k 'vdot) lv ((k 'vmap2) :mul lv lv)))
          (do
            (var s 0)
            (each v lvals (+= s (* v v v)))
            s)))
      (def a (/ sl3 (* 6 (math/pow sl2 1.5))))
      (def adj
        (map (fn [al]
               (def za (qnorm al 0 1 {}))
               (def wz (+ w za))
               (pnorm (+ w (/ wz (- 1 (* a wz)))) 0 1 {}))
             alpha))
      [(norm-inter t (adj 0)) (norm-inter t (adj 1))])
    (errorf "stats: boot-ci type %q not supported" type-kw)))

#==------------------------------------------------------------------------==#
# permutation test - two-sample difference of means. Exact full
# enumeration over all C(n1+n2, n1) splits (R combn arithmetic
# order); MC with the Davison-Hinkley add-one p = (#+1)/(B+1).
#==------------------------------------------------------------------------==#

(defn/typed comb-count {:args [:number :number] :ret :number} [n k]
  "C(n,k) via exact arithmetic, returned as a double."
  # C(n,k) = prod_{i=1..k} (n+1-i)/i via exact arithmetic, returned as a
  # double (ratios are not `=`-comparable to int/double; counts are
  # exact in doubles far beyond our enumeration limits)
  (var res 1)
  (for i 1 (inc k)
    (set res (* res (/ (- (inc n) i) i))))
  (math/round res))

(defn/typed next-comb {:args [:any :number] :ret :any} [comb n]
  "Next k-combination in lexicographic order; nil when exhausted.
  Internal to perm-test-two."
  (def k (length comb))
  (var i (dec k))
  (while (and (>= i 0) (= (get comb i) (- n (- k i))))
    (-- i))
  (if (< i 0)
    nil
    (let [out (array ;comb)]
      (put out i (inc (get comb i)))
      (for j (inc i) k (put out j (+ (get out (dec j)) 1)))
      out)))

(defn/typed perm-split-stat {:args [:any :any :number :number] :ret :number}
  [all comb n1 n2]
  "The split mean-difference statistic over one combination
  (R combn arithmetic order)."
  # sum over comb (ascending) + complement (ascending) - same arithmetic
  # order as the R combn oracle
  (var s1 0)
  (each i comb (+= s1 (get all i)))
  (var s2 0)
  (def k (length comb))
  (var ci 0)
  (for i 0 (+ n1 n2)
    (if (and (< ci k) (= (get comb ci) i))
      (++ ci)
      (+= s2 (get all i))))
  (- (/ s1 n1) (/ s2 n2)))

(defn perm-test-two
  "Two-sample permutation test on the difference of means. Exact
  full enumeration when C(n1+n2, n1) <= 20000 (or opts
  {:exact true}), else Monte Carlo. Returns {:statistic :p-value
  :exact :B :alternative}; the MC p uses the Davison-Hinkley
  add-one convention (cnt+1)/(B+1). Opts: {:exact nil :B 999 :seed
  \"stats\" :alternative :two-sided|:greater|:less :impl}. Throws on
  empty samples or :B not a positive integer."
  [a b &opt opts]
  (default opts {})
  # the exact branch enumerates tiny
  # per-combination sums (combinatorial - staging per comb would be
  # catastrophic); the MC branch's per-permutation masked sums would
  # need O(n) staging per permutation - same class as the pure loop
  # at stats-regime n (the kruskal scatter-add reasoning)
  (arrays-only opts a b)
  (def alt (get opts :alternative :two-sided))
  (def n1 (length a))
  (def n2 (length b))
  (when (or (< n1 1) (< n2 1)) (error "stats: perm-test-two needs non-empty samples"))
  (def ntot (+ n1 n2))
  (def all @[;a ;b])
  (def tobs (- (mean-plain a) (mean-plain b)))
  (def max-exact 20000)
  (def n-perm (comb-count ntot n1))
  (def exact? (if (nil? (get opts :exact)) (<= n-perm max-exact) (get opts :exact)))
  (if exact?
    (do
      (var cnt 0)
      (var comb (seq [i :range [0 n1]] i))
      (var done false)
      (while (not done)
        (def tp (perm-split-stat all comb n1 n2))
        (def ge
          (case alt
            :two-sided (>= (math/abs tp) (math/abs tobs))
            :greater (>= tp tobs)
            :less (<= tp tobs)))
        (when ge (++ cnt))
        (def nxt (next-comb comb ntot))
        (if nxt (set comb nxt) (set done true)))
      {:statistic tobs :p-value (/ cnt n-perm) :exact true :B n-perm
       :alternative alt})
    (do
      (def bb (get opts :B 999))
      (unless (and (int? bb) (pos? bb)) (error "stats: :B must be a positive integer"))
      (def rng (math/rng (string (get opts :seed "stats"))))
      (var cnt 0)
      (def pool @[;a ;b])
      (for _ 0 bb
        # reservoir-style subset draw: mark first n1 selected positions
        (def sel (array/new-filled ntot false))
        (var chosen 0)
        (while (< chosen n1)
          (def i (math/rng-int rng ntot))
          (unless (get sel i)
            (put sel i true)
            (++ chosen)))
        (var s1 0)
        (var s2 0)
        (for i 0 ntot
          (if (get sel i) (+= s1 (get pool i)) (+= s2 (get pool i))))
        (def tp (- (/ s1 n1) (/ s2 n2)))
        (def ge
          (case alt
            :two-sided (>= (math/abs tp) (math/abs tobs))
            :greater (>= tp tobs)
            :less (<= tp tobs)))
        (when ge (++ cnt)))
      # Davison-Hinkley add-one (documented OUR convention)
      {:statistic tobs :p-value (/ (inc cnt) (* 1.0 (inc bb))) :exact false :B bb
       :alternative alt})))
