# stats/tests.janet - classical tests. R conventions.
# Conformance: the recorded R 4.5.0 fixture CSVs under tests/fixtures/stats.
# t-tests delegate their moment sums to the descriptive kernels
# (mean/var/sd accept resident handles - t-test-one/two are
# RESIDENT under :num); chisq-gof's cell sums are kernel-backed
# (vmap2/vreduce, no validation pass - the NaN domain contract);
# wilcoxon-sr/mann-whitney put their STATISTIC sums through the
# kernels (indicator-masked rank sums - sorts and exact DPs stay
# janet, the documented composition); chisq-test/fisher/kruskal are
# composition (small tables, scatter-adds - recorded reasoning).
(import typed :prefix "")
(import ./utils :prefix "")
(import ./dist :prefix "")
(import ./descriptive :prefix "")
(import ./ranksum :prefix "")
(import ./num)

(defn arrays-only
  "Composition guard: resolves/validates the impl, then asserts
  janet-array inputs - sorts, table glue and exact DPs need
  janet-side values; a packed handle unpacks first (actionable)."
  [opts & args]
  (num/check-impl-landed opts)
  (each x args
    (unless (indexed? x)
      (error "stats: composition path consumes janet arrays; unpack packed handles first")))
  nil)

(defn/typed df-check {:args [:number] :ret nil} [df]
  "Assert df >= 1. Internal."
  (when (< df 1) (errorf "stats: df must be >= 1, got %q" df)))

#==------------------------------------------------------------------------==#
# t-tests + confidence intervals
#==------------------------------------------------------------------------==#

(defn t-stat-one
  "The one-sample t statistic (m = the precomputed mean). Internal."
  [xs mu m &opt opts]
  (default opts {})
  # one-sample t statistic; opts threads the impl to the variance sum
  (def n (num/count xs))
  (when (< n 2) (error "stats: t-test needs >= 2 observations"))
  (* (- m mu) (math/sqrt (/ n (variance xs m opts)))))

(defn t-test-one
  "One-sample t test, R t.test(x, mu=): slot 2 is the null mean
  (default 0; pass nil), opts table trailing. Returns {:estimate
  :statistic :parameter :p-value :conf-interval :stderr} - the
  R key set, two-sided p. Opts: {:conf-level 0.95 :impl}; xs may
  be a packed handle under :num (the resident t-test). Throws on
  a non-number mu or n < 2."
  [xs mu &opt opts]
  (default opts {})
  # the mu-slot trap (public - R's t.test second positional), guarded
  # actionably; the remedy: (default mu 0) AFTER the guard, the
  # wilcoxon-sr/mad pattern
  (when (and mu (not (number? mu)))
    (error "stats: t-test-one expected a number (the null mean) in slot 2 - pass nil for the default 0, opts table trailing"))
  (default mu 0)
  # R t.test(x, mu=): list estimate statistic p-value ci + df.
  # the moment sums ride the descriptive kernels (xs may be a
  # packed handle under :num - the resident t-test)
  (num/check-impl-landed opts xs)
  (def m (mean xs opts))
  (def n (num/count xs))
  (def df (- n 1))
  (df-check df)
  (def t (t-stat-one xs mu m opts))
  # two-sided p (the default; alternative in opts when needed - v1: two-sided)
  (def p (* 2 (- 1 (pt (math/abs t) df {}))))
  (def conf (get opts :conf-level 0.95))
  (def tcrit (qt (- 1 (* 0.5 (- 1 conf))) df {}))
  (def halfw (* tcrit (sd xs m opts) (math/sqrt (/ n))))
  {:estimate m :statistic t :parameter df :p-value p
   :conf-interval [(- m halfw) (+ m halfw)]
   :stderr (/ (sd xs m opts) (math/sqrt n))})

(defn t-test-two
  "Two-sample t test, R t.test: Welch by default, pooled with
  {:var-equal true}. Returns {:estimate-diff :statistic
  :parameter :p-value :conf-interval :mean-x :mean-y}, two-sided
  p. Opts: {:var-equal false :conf-level 0.95 :impl}; resident
  handles ok. Throws on df < 1."
  [xs ys &opt opts]
  (default opts {})
  # Welch (default, var.equal=false) or pooled (var.equal=true); R t.test.
  # moment sums via the descriptive kernels (resident handles ok)
  (num/check-impl-landed opts xs)
  (def nx (num/count xs))
  (def ny (num/count ys))
  (def mx (mean xs opts))
  (def my (mean ys opts))
  (def vx (variance xs mx opts))
  (def vy (variance ys my opts))
  (def var-equal (get opts :var-equal false))
  (def df (if var-equal
           (+ nx ny -2)
           (let [a (/ vx nx)
                 b (/ vy ny)
                 se2 (+ a b)]
             (/ (* se2 se2)
                (+ (/ (* a a) (- nx 1)) (/ (* b b) (- ny 1)))))))
  (df-check df)
  (def se (if var-equal
           (math/sqrt (* (/ (+ (* (- nx 1) vx) (* (- ny 1) vy)) df)
                         (+ (/ nx) (/ ny))))
           (math/sqrt (+ (/ vx nx) (/ vy ny)))))
  (def t (/ (- mx my) se))
  (def p (* 2 (- 1 (pt (math/abs t) df {}))))
  (def conf (get opts :conf-level 0.95))
  (def tcrit (qt (- 1 (* 0.5 (- 1 conf))) df {}))
  (def halfw (* tcrit se))
  {:estimate-diff (- mx my) :statistic t :parameter df :p-value p
   :conf-interval [(- (- mx my) halfw) (+ (- mx my) halfw)]
   :mean-x mx :mean-y my})

(defn t-test-paired
  "Paired t test, R t.test(x, y, paired=TRUE): the one-sample t on
  the differences. Under :num the diffs are the vmap2 :sub of the
  staged inputs and flow to t-test-one as a handle - never
  unpacked. Opts pass through to t-test-one."
  [xs ys &opt opts]
  (default opts {})
  # R t.test(x, y, paired=TRUE): one-sample t on the differences.
  # resident: under :num the diffs are the vmap2 :sub of the
  # staged inputs and flow to t-test-one as a handle - never unpacked
  (if (and (= :num (num/resolve-impl opts))
           (not (or (number? xs) (number? ys))))
    (let [k (num/kernels)]
      (t-test-one ((k 'vmap2) :sub (num/stage xs k) (num/stage ys k)) 0 opts))
    (do
      (num/check-impl-landed opts xs)
      (def ds (seq [i :range [0 (length xs)]] (- (xs i) (ys i))))
      (t-test-one ds 0 opts))))

#==------------------------------------------------------------------------==#
# chi-square: GOF + contingency
#==------------------------------------------------------------------------==#

(defn chisq-gof
  "Chi-squared goodness of fit. expected = :uniform (calculated)
  or a vector of expected counts. Returns {:statistic :parameter
  :p-value}. Kernel-backed cell sums (d^2/e via vmap2+vreduce,
  resident handles chain; vector+num does NO validation pass -
  e <= 0 gives inf/NaN, the NaN domain contract; the pure path
  keeps its actionable throws). Opts: {:impl}; throws on n < 2
  or a length mismatch."
  [observed expected &opt opts]
  (default opts {})
  # goodness of fit; expected can be :uniform (calculated) or a
  # vector
  (def n (num/count observed))
  (when (< n 2) (error "stats: chisq-gof needs >= 2 cells"))
  (def exp (if (= expected :uniform)
             (seq [_ :range [0 n]] (/ (reduce + 0 observed) n))
             expected))
  (unless (= (num/count exp) n) (error "stats: chisq-gof expected/observed length mismatch"))
  (def df (- n 1))
  (if (= :num (num/resolve-impl opts))
    (let [k (num/kernels)
          of (num/stage observed k)
          ef (num/stage exp k)
          d ((k 'vmap2) :sub of ef)]
      (def x2 ((k 'vreduce) :sum ((k 'vmap2) :div ((k 'vmap2) :mul d d) ef)))
      {:statistic x2 :parameter df :p-value (- 1 (pchisq x2 df {}))})
    (do
      (num/check-impl-landed opts observed)
      (var x2 0)
      (for i 0 n
        (when (<= (exp i) 0) (errorf "stats: chisq-gof expected must be > 0, got %q" (exp i)))
        (+= x2 (/ (* (- (observed i) (exp i)) (- (observed i) (exp i))) (exp i))))
      {:statistic x2 :parameter df :p-value (- 1 (pchisq x2 df {}))})))

(defn chisq-test
  "Chi-squared contingency-table test, R chisq.test: the
  independence X^2 over row x column expected cells. Returns
  {:statistic :parameter :p-value :expected}. Composition
  (tables are small in the stats regime - flattening for the
  kernels would cost more than the cell loop). Opts: {:impl};
  throws on < 2 rows/cols or a ragged table."
  [table &opt opts]
  (default opts {})
  # contingency tables are small in the
  # stats regime - flattening for the kernels costs more than the
  # cell loop; the janet table glue IS the implementation
  (arrays-only opts table)
  (def nrows (length table))
  (when (< nrows 2) (error "stats: chisq-test needs >= 2 rows"))
  (def ncols (length (table 0)))
  (each r table (unless (= (length r) ncols) (error "stats: chisq-test ragged table")))
  (when (< ncols 2) (error "stats: chisq-test needs >= 2 cols"))
  (def row-sums (map (fn [r] (reduce + 0 r)) table))
  (var grand 0) (each s row-sums (+= grand s))
  (def col-sums (seq [c :range [0 ncols]] (var s 0) (each r table (+= s (r c))) s))
  (var x2 0)
  (for r 0 nrows
    (for c 0 ncols
      (def e (* (row-sums r) (/ (col-sums c) grand)))
      (def d (- (get (table r) c) e))
      (+= x2 (/ (* d d) e))))
  (def df (* (- nrows 1) (- ncols 1)))
  {:statistic x2 :parameter df :p-value (- 1 (pchisq x2 df {}))
   :expected (seq [r :range [0 nrows]]
                 (seq [c :range [0 ncols]]
                      (* (row-sums r) (/ (col-sums c) grand))))})

#==------------------------------------------------------------------------==#
# Fisher exact 2x2
#==------------------------------------------------------------------------==#

(defn/typed log-hyper-geo {:args [:number :number :number :number] :ret :number} [a b c d]
  "log P of the 2x2 table with these margins (Fisher's exact core);
  margins: row1=a+b, row2=c+d, col1=a+c, col2=b+d, grand.
  Internal to fisher-test-2x2."
  (- (+ (lgamma (+ a b 1)) (lgamma (+ c d 1))
        (lgamma (+ a c 1)) (lgamma (+ b d 1)))
     (+ (lgamma (+ a 1)) (lgamma (+ b 1))
        (lgamma (+ c 1)) (lgamma (+ d 1))
        (lgamma (+ a b c d 1)))))

(defn/typed fisher-c-mle-or {:args [:number :number :number :number] :ret :number} [a b c d]
  "Conditional MLE of the odds ratio - R fisher.test's $estimate:
  solve E_theta[X] = a under the noncentral hypergeometric with
  odds ratio theta (E is monotone in theta; solve-monotone on log
  theta). Internal to fisher-test-2x2."
  (def row1 (+ a b))
  (def col1 (+ a c))
  (def grand (+ a b c d))
  (def lo (max 0 (- (+ row1 col1) grand)))
  (def hi (min row1 col1))
  (defn/typed mean-at {:args [:number] :ret :number} [log-theta]
    (var wsum 0)
    (var xwsum 0)
    (for x lo (+ hi 1)
      (def lp (+ (- (+ (lgamma (+ row1 1)) (lgamma (+ (- grand row1) 1))
                       (lgamma (+ col1 1)) (lgamma (+ (- grand col1) 1)))
                    (+ (lgamma (+ x 1)) (lgamma (+ (- row1 x) 1))
                       (lgamma (+ (- col1 x) 1))
                       (lgamma (+ (- grand (- (+ row1 col1) x)) 1))))
                 (* x log-theta)))
      (def w (math/exp lp))
      (+= wsum w)
      (+= xwsum (* x w)))
    (/ xwsum wsum))
  (when (or (= a lo) (= a hi))
    (if (= a lo) 0 math/inf))
  (def log-or (solve-monotone mean-at -40 40 a))
  (math/exp log-or))

(defn fisher-test-2x2
  "Fisher's exact test on a 2x2 table [[a b] [c d]], R
  fisher.test: two-sided p = the sum of all tables with P <= the
  observed P. Returns {:p-value :odds-ratio :estimate} - the
  SAMPLE odds ratio kept separate from :estimate, the conditional
  MLE (R's $estimate). Opts: {:impl}."
  [table &opt opts]
  (default opts {})
  # table = [[a b] [c d]]; two-sided (sum of all tables
  # with P <= observed P, standard network-style two-sided).
  # Composition: 2x2 scalar + the small support loop
  (arrays-only opts table)
  (def a (get (table 0) 0))
  (def b (get (table 0) 1))
  (def c (get (table 1) 0))
  (def d (get (table 1) 1))
  (def row1 (+ a b))
  (def col1 (+ a c))
  (def grand (+ a b c d))
  # support of x (the [1,1] cell): max(0, row1+col1-grand) .. min(row1, col1)
  (def lo (max 0 (- (+ row1 col1) grand)))
  (def hi (min row1 col1))
  (def lp0 (log-hyper-geo a b c d))
  (var lsum (- math/inf))
  (def lps @[])
  (for x lo (+ hi 1)
    (def lp (log-hyper-geo x (- row1 x) (- col1 x) (- grand (+ row1 col1 (- x)))))
    (array/push lps lp)
    # two-sided: include tables with P <= observed (with small tolerance)
    (when (<= lp (+ lp0 1e-9))
      (set lsum (if (= lsum (- math/inf)) lp (logsumexp @[lsum lp])))))
  {:p-value (math/exp lsum)
   # R fisher.test: $estimate = conditional MLE; sample OR kept separate
   :odds-ratio (if (zero? (* b c)) math/inf (/ (* a d) (* b c)))
   :estimate (fisher-c-mle-or a b c d)})

#==------------------------------------------------------------------------==#
# rank tests
#==------------------------------------------------------------------------==#

(defn wilcoxon-sr
  "Wilcoxon signed-rank test, R wilcox.test(x, mu=): exact
  (psignrank DP) when no ties and n < 50, else R's normal
  approximation with tie correction + continuity (V scale,
  sigma /24-/48, cc = sign(z)*0.5 - R's formulas verbatim).
  Slot 2 is the null center (default 0; zero differences dropped).
  Returns {:statistic :n :p-value}. Opts: {:impl}; throws on a
  non-number mu or all-zero differences."
  [xs &opt mu opts]
  # R wilcox.test(x, mu=): exact (psignrank DP) when no ties and n < 50;
  # else normal approximation with tie correction + continuity - R's
  # formulas verbatim (V scale, sigma /24-/48, cc = sign(z)*0.5).
  # the |d| sort + rank ties + exact DPs stay janet (the
  # composition); the V sum goes through the kernels (indicator-
  # masked rank sum - adding 0s is exact, so parity is bitwise)
  (default mu 0)
  (default opts {})
  # the mu-slot trap, guarded
  (when (and mu (not (number? mu)))
    (error "stats: wilcoxon-sr expected a number (the null mean) in slot 2 - pass nil for the default 0, opts table trailing"))
  (arrays-only opts xs)
  (def diffs (map |(- $ mu) xs))
  (def nz (filter |(not= $ 0) diffs))
  (def n (length nz))
  (when (< n 1) (error "stats: wilcoxon-sr needs non-zero differences"))
  (def with-abs (seq [i :range [0 n]] [(math/abs (nz i)) i]))
  (def sorted-pts (sorted with-abs (fn [p q] (< (p 0) (q 0)))))
  (def ranks (rank-ties sorted-pts))
  (def v
    (if (= :num (num/resolve-impl opts))
      (let [k (num/kernels)
            pos (seq [i :range [0 n]] (if (pos? (nz i)) 1.0 0.0))]
        ((k 'vreduce) :sum ((k 'vmap2) :mul (num/stage ranks k) (num/stage pos k))))
      (do
        (var vv 0)
        (for i 0 n
          (when (pos? (nz i)) (+= vv (get ranks i))))
        vv)))
  (def ties? (< (length (rank-multiplicities (map first sorted-pts))) n))
  (def p-value
    (if (and (< n 50) (not ties?))
      (if (> v (* n (+ n 1) 0.25))
        (min 1 (* 2 (signrank-p (- v 1) n false)))
        (min 1 (* 2 (signrank-p v n true))))
      (let [ts (rank-multiplicities (map first sorted-pts))
            z0 (- v (* n (+ n 1) 0.25))
            sigma (math/sqrt (- (/ (* n (+ n 1) (+ (* 2 n) 1)) 24)
                                (/ (sum-t3-minus-t ts) 48)))
            z (/ (- z0 (* (if (>= z0 0) 1 (if (<= z0 0) -1 0)) 0.5)) sigma)
            pz (pnorm z 0 1 {})]
        (min 1 (* 2 (min pz (- 1 pz)))))))
  {:statistic v :n n :p-value p-value})

(defn mann-whitney
  "Mann-Whitney U / Wilcoxon rank-sum, R wilcox.test(x, y):
  exact (pwilcox DP) when no ties and n < 50, else R's normal
  approximation with tie correction + cc. Returns {:statistic
  :p-value :r-x} - W = R_x - nx(nx+1)/2 (R's statistic).
  Opts: {:impl}."
  [xs &opt ys opts]
  # R wilcox.test(x, y): exact (pwilcox DP) when no ties and n < 50;
  # else normal with tie correction + cc - R's formulas verbatim.
  # sorts/rank-ties/exact DPs stay janet; the R_x statistic sum
  # goes through the kernels (group-indicator masked rank sum)
  (default opts {})
  (arrays-only opts xs ys)
  (def nx (length xs))
  (def ny (length ys))
  (def n (+ nx ny))
  (def pts @[])
  (def group-of @[])
  (for i 0 nx (array/push pts [(xs i) (length pts)]) (array/push group-of 0))
  (for i 0 ny (array/push pts [(ys i) (length pts)]) (array/push group-of 1))
  (def sorted-pts (sorted pts (fn [p q] (< (p 0) (q 0)))))
  (def ranks (rank-ties sorted-pts))
  (def r-x
    (if (= :num (num/resolve-impl opts))
      (let [k (num/kernels)
            gi (seq [i :range [0 n]] (if (= (get group-of i) 0) 1.0 0.0))]
        ((k 'vreduce) :sum ((k 'vmap2) :mul (num/stage ranks k) (num/stage gi k))))
      (do
        (var rr 0)
        (for i 0 n
          (when (= (get group-of i) 0)
            (+= rr (get ranks i))))
        rr)))
  # W = R_x - nx(nx+1)/2 (R's statistic)
  (def w (- r-x (* 0.5 nx (+ nx 1))))
  (def ties? (< (length (rank-multiplicities (map first sorted-pts))) n))
  (def p-value
    (if (and (< n 50) (not ties?))
      (if (> w (* nx ny 0.5))
        (min 1 (* 2 (wilcox-p (- w 1) nx ny false)))
        (min 1 (* 2 (wilcox-p w nx ny true))))
      (let [ts (rank-multiplicities (map first sorted-pts))
            z0 (- w (* nx ny 0.5))
            sigma (math/sqrt (* (/ (* nx ny) 12)
                                (- (+ n 1)
                                   (/ (sum-t3-minus-t ts) (* n (- n 1))))))
            z (/ (- z0 (* (if (>= z0 0) 1 (if (<= z0 0) -1 0)) 0.5)) sigma)
            pz (pnorm z 0 1 {})]
        (min 1 (* 2 (min pz (- 1 pz)))))))
  {:statistic w :p-value p-value :r-x r-x})

(defn kruskal-wallis
  "Kruskal-Wallis rank-sum test, R kruskal.test: groups = array
  of sample arrays. Returns {:statistic :parameter :p-value
  :rank-sums} - H with R's tie correction (H / (1 -
  sum(t^3-t)/(n^3-n))). Composition (per-group rank sums are a
  scatter-add - no honest kernel shape). Opts: {:impl}; throws on
  < 2 groups."
  [groups &opt opts]
  (default opts {})
  # groups = array of sample arrays; R kruskal.test.
  # Composition: the per-group rank sums are a scatter-add - no
  # honest kernel shape (per-group masking would be O(kn)); the janet
  # accumulation IS the implementation
  (arrays-only opts groups)
  (def k (length groups))
  (when (< k 2) (error "stats: kruskal-wallis needs >= 2 groups"))
  (def pts @[])
  (def sizes @[])
  (each g groups
    (array/push sizes (length g))
    (each v g (array/push pts [v (length pts)])))
  (def n (length pts))
  (def sorted-pts (sorted pts (fn [p q] (< (p 0) (q 0)))))
  (def ranks (rank-ties sorted-pts))
  # rank sums per group (group id = the index of the sample the point
  # came from - track via a parallel array)
  (def group-ids @[])
  (eachp [gi g] groups
    (each v g (array/push group-ids gi)))
  (def rank-sums (array/new-filled k 0))
  (for i 0 n
    (def p (get sorted-pts i))
    (def orig-idx (get p 1))
    (+= (rank-sums (get group-ids orig-idx)) (get ranks orig-idx)))
  # H = 12/(n(n+1)) * sum(R_i^2/n_i) - 3(n+1), tie-corrected
  (var h 0)
  (for g 0 k
    (+= h (/ (* (rank-sums g) (rank-sums g)) (sizes g))))
  (set h (- (* 12 (/ h (* n (+ n 1)))) (* 3 (+ n 1))))
  # tie correction factor
  (var sum-t3 0)
  (var i 0)
  (while (< i n)
    (def v (get (get sorted-pts i) 0))
    (var j i)
    (while (and (< j n) (= (get (get sorted-pts j) 0) v)) (++ j))
    (def t (- j i))
    (+= sum-t3 (- (* t t t) t))
    (set i j))
  # R's tie correction: H / (1 - Σ(t³-t)/(n³-n))
  (when (pos? sum-t3)
    (def n3-n (- (* n n n) n))
    (set h (* h (/ n3-n (- n3-n sum-t3)))))
  {:statistic h :parameter (- k 1)
   :p-value (- 1 (pchisq h (- k 1) {}))
   :rank-sums rank-sums})

#==------------------------------------------------------------------------==#
# Kolmogorov-Smirnov. R ks.test conventions decoded from the
# R 4.5.0 source: ks.test.default (the walk + tie-block collapse,
# the exact switches), psmirnov/pkolmogorov (the R wrappers), and
# src/library/stats/src/ks.c (C_psmirnov_exact - the Viehmann DP
# with the tie grid; C_pkolmogorov_two_exact - the MTW matrix
# power; C_pkolmogorov_two_limit - K2l). Value parity with R
# requires the identical fp term sequence - every engine below
# is a verbatim port. Dual-impl: the two-sample WALK is the vscan
# composition (cumulative +-1/n over the union grid), D via
# vmap1 :abs + vreduce; sorts and the exact DPs stay janet (the
# documented composition). NO new kernels.
#==------------------------------------------------------------------------==#

(def PKOLM-INV-SQRT-2PI 0.398942280401432677939946059934)  # M_1_SQRT_2PI
(def PKOLM-PI2 1.570796326794896619231321691640)           # M_PI_2
(def PKOLM-PI4 0.785398163397448309615660845820)           # M_PI_4

(defn/typed ipow {:args [:number :number] :ret :number} [x e]
  "R_pow_di analogue: integer power by repeated multiplication
  (negative exponents invert). Internal."
  (var r 1)
  (var base x)
  (var n e)
  (when (< n 0) (set base (/ 1 base)) (set n (- n)))
  (for i 0 n (set r (* r base)))
  r)

(defn/typed k2l {:args [:number :boolean :number] :ret :number} [x lower tol]
  "ks.c K2l verbatim: the two-sided asymptotic Kolmogorov
  distribution (the alternative series for x < 1, the standard
  alternating series otherwise). Internal to ks-test."
  (cond
    (<= x 0) (if lower 0 1)
    (< x 1) (let [k-max (math/floor (math/sqrt (- 2 (math/log tol))))
                  w (math/log x)
                  z (- (/ (* PKOLM-PI2 PKOLM-PI4) (* x x)))]
              (var acc 0)
              (var k 1)
              (while (< k k-max)
                (+= acc (math/exp (- (* k k z) w)))
                (+= k 2))
              (if lower (/ acc PKOLM-INV-SQRT-2PI)
                (- 1 (/ acc PKOLM-INV-SQRT-2PI))))
    (let [z (* -2 x x)]
      (var s -1)
      (var k (if lower 1 2))
      (var old 0)
      (var new (if lower 1 (* 2 (math/exp z))))
      (while (> (math/abs (- old new)) tol)
        (set old new)
        (set new (+ new (* 2 s (math/exp (* z k k)))))
        (set s (* s -1))
        (++ k))
      new)))

(defn/typed m-multiply {:args [:any :any :any :number] :ret :any} [a b c m]
  "ks.c m-multiply verbatim (flat m x m matrices). Internal."
  (for i 0 m
    (for j 0 m
      (var s 0)
      (for kk 0 m (+= s (* (a (+ (* i m) kk)) (b (+ (* kk m) j)))))
      (set (c (+ (* i m) j)) s))))

(defn/typed m-power {:args [:any :number :any :number :number :number] :ret :any}
  [a e-a v e-v m n]
  "ks.c m_power verbatim (returns the updated [V eV] pair).
  Internal."
  (if (= n 1)
    [a e-a]
    (let [[v e-v2] (m-power a e-a v e-v m (math/floor (/ n 2)))
          b (array/new-filled (* m m) 0)]
      (m-multiply v v b m)
      (def e-b (* 2 e-v2))
      (def [v2 e-v3]
        (if (even? n)
          [b e-b]
          (let [v2 (array/new-filled (* m m) 0)]
            (m-multiply a b v2 m)
            [v2 (+ e-a e-b)])))
      (def e-v4
        (if (> (get v2 (+ (* (math/floor (/ m 2)) m) (math/floor (/ m 2)))) 1e140)
          (do (for i 0 (* m m) (set (v2 i) (* (v2 i) 1e-140)))
              (+ e-v3 140))
          e-v3))
      [v2 e-v4])))

(defn/typed k2x {:args [:number :number] :ret :number} [n d]
  "ks.c K2x verbatim: Marsaglia-Tsang-Wang 'Evaluating Kolmogorov's
  distribution' (JSS 8/i18) - the one-sample two-sided exact
  distribution via matrix power with 1e140 exponent scaling.
  Internal to ks-test."
  (def k (+ (math/floor (* n d)) 1))
  (def m (- (* 2 k) 1))
  (def h (- k (* n d)))
  (def hh (array/new-filled (* m m) 0))
  (for i 0 m
    (for j 0 m
      (set (hh (+ (* i m) j)) (if (< (+ 1 (- i j)) 0) 0 1))))
  (for i 0 m
    (-= (hh (+ (* i m) 0)) (ipow h (+ i 1)))
    (-= (hh (+ (* (- m 1) m) i)) (ipow h (- m i))))
  (+= (hh (* (- m 1) m))
      (if (> (- (* 2 h) 1) 0) (ipow (- (* 2 h) 1) m) 0))
  (for i 0 m
    (for j 0 m
      (when (> (+ 1 (- i j)) 0)
        (for g 1 (+ 2 (- i j))
          (set (hh (+ (* i m) j)) (/ (hh (+ (* i m) j)) g))))))
  (var [q e-q] (m-power hh 0 (array/new-filled (* m m) 0) 0 m n))
  (var s (get q (+ (* (- k 1) m) (- k 1))))
  (for i 1 (+ n 1)
    (set s (/ (* s i) n))
    (when (< s 1e-140)
      (set s (* s 1e140))
      (-= e-q 140)))
  (* s (ipow 10 e-q)))

(defn/typed pkolm-one-exact {:args [:number :number] :ret :number} [q n]
  "R pkolmogorov_one_exact verbatim (the upper-tail closed form).
  Internal to ks-test."
  (var p 0)
  (def j-max (math/floor (* n (- 1 q))))
  (for j 0 (+ j-max 1)
    (def lc (- (lgamma (+ n 1)) (lgamma (+ j 1)) (lgamma (+ (- n j) 1))))
    (def t1 (* (- n j) (math/log (- (- 1 q) (/ j n)))))
    (def t2 (* (- j 1) (math/log (+ q (/ j n)))))
    (+= p (math/exp (+ lc t1 t2))))
  (* q p))

(defn/typed psmirnov-test {:args [:number :number :number :boolean] :ret :boolean}
  [q r s two]
  "ks.c psmirnov test boundary predicate. Internal."
  (if two (>= (math/abs (- r s)) q) (>= (- r s) q)))

(defn/typed psmirnov-lower {:args [:number :number :number :any :boolean] :ret :number}
  [q m n z two]
  "ks.c psmirnov_exact_[uniq|ties]_lower verbatim: p = P(D < q)
  via the Viehmann weights; z non-nil gates the tie-aware boundary
  zeroing. Internal."
  (def u (array/new-filled (+ n 1) 0))
  (set (u 0) 1)
  (for j 1 (+ n 1)
    (set (u j) (if (and (psmirnov-test q 0 (/ j n) two)
                        (or (not z) (= 1 (get z j))))
                 0
                 (u (- j 1)))))
  (for i 1 (+ m 1)
    (def w (/ i (+ i n)))
    (set (u 0) (if (and (psmirnov-test q (/ i m) 0 two)
                        (or (not z) (= 1 (get z i))))
                 0
                 (* w (u 0))))
    (for j 1 (+ n 1)
      (set (u j) (if (and (psmirnov-test q (/ i m) (/ j n) two)
                          (or (not z) (= 1 (get z (+ i j)))))
                   0
                   (+ (* w (u j)) (u (- j 1)))))))
  (u n))

(defn/typed psmirnov-upper {:args [:number :number :number :any :boolean] :ret :number}
  [q m n z two]
  "ks.c psmirnov_exact_[uniq|ties]_upper verbatim: p = P(D >= q).
  Internal."
  (def u (array/new-filled (+ n 1) 0))
  (set (u 0) 0)
  (for j 1 (+ n 1)
    (set (u j) (if (and (psmirnov-test q 0 (/ j n) two)
                        (or (not z) (= 1 (get z j))))
                 1
                 (u (- j 1)))))
  (for i 1 (+ m 1)
    (when (and (psmirnov-test q (/ i m) 0 two)
               (or (not z) (= 1 (get z i))))
      (set (u 0) 1))
    (for j 1 (+ n 1)
      (set (u j) (if (and (psmirnov-test q (/ i m) (/ j n) two)
                          (or (not z) (= 1 (get z (+ i j)))))
                   1
                   (let [v (/ i (+ i j))
                         w (/ j (+ i j))]
                     (+ (* v (u j)) (* w (u (- j 1)))))))))
  (u n))

(defn/typed psmirnov-exact {:args [:number :number :number :any :boolean :boolean] :ret :number}
  [q m n z two lower]
  "ks.c psmirnov_exact verbatim incl. the q normalization.
  Internal."
  (def qq (/ (+ 0.5 (math/floor (- (* q m n) 1e-7))) (* m n)))
  (if lower (psmirnov-lower qq m n z two) (psmirnov-upper qq m n z two)))

(defn ks-cdf-fn
  "Resolve the one-sample reference: y keyword -> the landed CDF
  (:pnorm :pt :pchisq :pf :pbeta :pgamma, params via opts, R
  defaults); y function -> itself (must accept an array of points,
  R's vectorized-CDF convention). Internal to ks-test."
  [y &opt opts]
  (default opts {})
  (case y
    :pnorm (let [m (get opts :mean 0)
                 s (get opts :sd 1)]
             (fn [x] (pnorm x m s opts)))
    :pt (let [d (get opts :df)]
          (unless d (error "stats: ks-test needs :df for :pt"))
          (fn [x] (pt x d opts)))
    :pchisq (let [d (get opts :df)]
              (unless d (error "stats: ks-test needs :df for :pchisq"))
              (fn [x] (pchisq x d opts)))
    :pf (let [d1 (get opts :df1)
              d2 (get opts :df2)]
          (unless (and d1 d2) (error "stats: ks-test needs :df1/:df2 for :pf"))
          (fn [x] (pf x d1 d2 opts)))
    :pbeta (let [a (get opts :shape1)
                 b (get opts :shape2)]
             (unless (and a b)
               (error "stats: ks-test needs :shape1/:shape2 for :pbeta"))
             (fn [x] (pbeta x a b opts)))
    :pgamma (let [sh (get opts :shape)
                  r (get opts :rate)]
              (unless (and sh r)
                (error "stats: ks-test needs :shape/:rate for :pgamma"))
              (fn [x] (pgamma x sh r opts)))
    y))

(defn ks-check-finite
  "Assert all observations finite. Internal."
  [xs]
  (each v xs
    (unless (and (< v math/inf) (> v (- math/inf)))
      (error "stats: ks-test needs finite observations")))
  nil)

(defn ks-two-sample
  "The two-sample KS engine: the cumulative +-1/n walk over the
  union grid with tie-block collapse (D only measured at block
  ends), exact p via the Viehmann DP when nx*ny < 10000 (or
  :exact true), else the asymptotic K2l / exponential forms.
  Returns {:statistic :p-value :alternative :exact :ties}.
  Internal to ks-test."
  [x y alt &opt opts]
  (default opts {})
  (def nx (length x))
  (def ny (length y))
  (def n-tot (+ nx ny))
  (def sp (sorted (do (def pts @[])
                      (for i 0 nx (array/push pts [(x i) 1]))
                      (for i 0 ny (array/push pts [(y i) 0]))
                      pts)
                  (fn [p q] (< (p 0) (q 0)))))
  (def vals (seq [i :range [0 n-tot]]
               (if (= 1 ((sp i) 1)) (/ 1 nx) (- (/ 1 ny)))))
  # the walk + tie-block collapse: z counts ONLY at block ends (the
  # whole-block convention; intermediate positions discarded)
  (def block-end?
    (seq [i :range [0 n-tot]]
      (if (= i (- n-tot 1))
        true
        (not= ((sp i) 0) ((sp (+ i 1)) 0)))))
  (def kept
    (if (= :num (num/resolve-impl opts))
      (let [k (num/kernels)
            zc ((k 'unpack-f64) ((k 'vscan) :sum (num/stage vals k)))]
        (seq [i :range [0 n-tot] :when (block-end? i)] (zc i)))
      (do (var cur 0)
          (def out @[])
          (for i 0 n-tot
            (+= cur (vals i))
            (when (block-end? i) (array/push out cur)))
          out)))
  (def d
    (if (= :num (num/resolve-impl opts))
      (let [k (num/kernels)
            kh (num/stage kept k)]
        (case alt
          :two-sided ((k 'vreduce) :max ((k 'vmap1) :abs kh))
          :greater ((k 'vreduce) :max kh)
          :less (- ((k 'vreduce) :min kh))))
      (case alt
        :two-sided (do (var mx 0)
                       (each v kept (set mx (max mx (math/abs v))))
                       mx)
        :greater (do (var mx (- math/inf))
                     (each v kept (set mx (max mx v)))
                     mx)
        :less (do (var mn math/inf)
                  (each v kept (set mn (min mn v)))
                  (- mn)))))
  (def ties? (some (fn [i] (not (block-end? i))) (range 0 (- n-tot 1))))
  (def exact?
    (if (nil? (get opts :exact))
      (< (* nx ny) 10000)
      (get opts :exact)))
  # the tie grid for the exact DP (R: c(0L, diff(sort(w)) != 0, 1L))
  (def zgrid
    (when ties?
      (def zg (array/new-filled (+ n-tot 1) 0))
      (for i 1 n-tot
        (set (zg i) (if (= ((sp (- i 1)) 0) ((sp i) 0)) 0 1)))
      (set (zg n-tot) 1)
      zg))
  (def p
    (cond
      (<= d 0) 1
      (> d 1) 0
      exact?
      (let [[m sm-n] (if (= alt :less) [ny nx] [nx ny])]
        (psmirnov-exact d m sm-n zgrid (= alt :two-sided) false))
      (= alt :two-sided)
      (k2l (* (math/sqrt (/ (* nx ny) n-tot)) d) false 1e-6)
      (math/exp (* (* -2 (/ (* nx ny) n-tot)) (* d d)))))
  {:statistic d :p-value (min 1 (max 0 p)) :alternative alt
   :exact exact? :ties ties?})

(defn ks-one-sample
  "The one-sample KS engine against the resolved CDF (see
  ks-cdf-fn): D over the sorted points vs the CDF, exact p (K2x /
  pkolmogorov_one_exact) when n < 100 and no ties, else the
  asymptotic forms. Returns {:statistic :p-value :alternative
  :exact :ties}. Internal to ks-test."
  [x y alt &opt opts]
  (default opts {})
  (def n (length x))
  (def ties? (< (length (distinct x)) n))
  (def f0 (ks-cdf-fn y opts))
  (unless (or (function? f0) (cfunction? f0))
    (error "stats: ks-test y must be an array, a keyword naming a landed CDF, or a CDF function"))
  (def sx (sorted-copy x))
  (def f0v (f0 sx))
  (unless (indexed? f0v)
    (error "stats: ks-test's CDF must accept an array of points (R's vectorized-CDF convention)"))
  (def idx (seq [i :range [0 n]] (/ i n)))
  (def terms
    (if (= :num (num/resolve-impl opts))
      (let [k (num/kernels)]
        ((k 'unpack-f64) ((k 'vmap2) :sub (num/stage f0v k) (num/stage idx k))))
      (seq [i :range [0 n]] (- (f0v i) (idx i)))))
  (def d
    (case alt
      :two-sided
      (if (= :num (num/resolve-impl opts))
        (let [k (num/kernels)
              shift ((k 'unpack-f64) ((k 'vmap1-scalar) :sub (num/stage terms k) (/ 1 n)))
              lo ((k 'unpack-f64) ((k 'vmap1) :neg (num/stage shift k)))]
          (max ((k 'vreduce) :max (num/stage terms k))
               ((k 'vreduce) :max (num/stage lo k))))
        (do (var mx (- math/inf))
            (each t terms
              (set mx (max mx t (- (/ 1 n) t))))
            mx))
      :greater
      (if (= :num (num/resolve-impl opts))
        (let [k (num/kernels)]
          ((k 'vreduce) :max
           ((k 'vmap1) :neg
            ((k 'vmap1-scalar) :sub (num/stage terms k) (/ 1 n)))))
        (do (var mx (- math/inf))
            (each t terms (set mx (max mx (- (/ 1 n) t))))
            mx))
      :less
      (if (= :num (num/resolve-impl opts))
        (let [k (num/kernels)]
          ((k 'vreduce) :max (num/stage terms k)))
        (do (var mx (- math/inf))
            (each t terms (set mx (max mx t)))
            mx))))
  (def exact?
    (if (nil? (get opts :exact))
      (and (< n 100) (not ties?))
      (get opts :exact)))
  (def p
    (cond
      (<= d 0) 1
      (> d 1) 0
      (= alt :two-sided)
      (if exact? (- 1 (k2x n d)) (k2l (* (math/sqrt n) d) false 1e-6))
      (if exact? (pkolm-one-exact d n) (math/exp (* (* -2 n) (* d d))))))
  {:statistic d :p-value (min 1 (max 0 p)) :alternative alt
   :exact exact? :ties ties?})

(defn ks-test
  "Kolmogorov-Smirnov test, R ks.test(x, y, alternative, exact):
  two-sample (y an array) or one-sample against a CDF (y a
  keyword naming a landed distribution - :pnorm :pt :pchisq :pf
  :pbeta :pgamma, params via opts - or a CDF function).
  Returns {:statistic :p-value :alternative :exact :ties}
  (R-convention keys; :ties carries R's tie warning). Opts:
  {:alternative :two-sided (default)/:less/:greater :exact (nil =
  auto: two-sample nx*ny < 10000, one-sample n < 100 and no ties)
  :impl + the CDF params}. Throws on empty or non-finite input."
  [x y &opt opts]
  # R ks.test(x, y, alternative, exact)
  (default opts {})
  (arrays-only opts x)
  (when (zero? (length x))
    (error "stats: ks-test needs at least one observation"))
  (ks-check-finite x)
  (def alt (get opts :alternative :two-sided))
  (unless (or (= alt :two-sided) (= alt :less) (= alt :greater))
    (errorf "stats: ks-test :alternative must be :two-sided/:less/:greater, got %q" alt))
  (if (indexed? y)
    (do (arrays-only opts y)
        (when (zero? (length y))
          (error "stats: ks-test needs at least one observation in y"))
        (ks-check-finite y)
        (ks-two-sample x y alt opts))
    (ks-one-sample x y alt opts)))

#==------------------------------------------------------------------------==#
# p.adjust + the exact binomial test (the sign test's base) + the
# Poisson rate test + the sign test. R conventions decoded from
# the R 4.5.0 sources (p.adjust / binom.test / poisson.test). The
# exact tails use the accurate identities through the landed
# machinery (complements of near-1 recursions would lose the far
# tail; R's lower.tail=FALSE is STRICTLY greater):
#   binomial  P(K > k) = I_p(k+1, n-k)             [pbeta]
#   poisson   P(K >= k) = gamma-p(k, m)            [the LOWER regularized]
# NO new kernels; p-adjust is the sort-bound + vscan composition.
#==------------------------------------------------------------------------==#

(defn/typed stable-order {:args [:any :boolean] :ret :array} [ps desc?]
  "The (p, index) sort - ties keep ORIGINAL index order (R's
  stable order; explicit, independent of janet's sort stability).
  Internal to p-adjust."
  (sorted (seq [i :range [0 (length ps)]] i)
          (fn [a b]
            (if (= (ps a) (ps b))
              (if desc? (> a b) (< a b))
              (if desc? (> (ps a) (ps b)) (< (ps a) (ps b)))))))

(defn p-adjust
  "Adjust p-values for multiple comparisons, R p.adjust: methods
  :bonferroni :holm :hochberg :bh :by :none (+ :fdr alias -> BH).
  Returns the adjusted array (length preserved, R's stable-tie
  ordering). Opts: {:method :holm :n (default = length) :impl};
  p in [0,1] is validated (our convention - R is silent) and the
  method is validated before any short-circuit (R errors on an
  unknown method at any length)."
  [ps &opt opts]
  # R p.adjust; the n >= lp contract. Sort-bound composition: the
  # (p, index) sort is janet; the num lane stages the products +
  # the cummin/cummax via vscan.
  (default opts {})
  (def lp (length ps))
  (def n (get opts :n lp))
  (def method0 (get opts :method :holm))
  (def method (if (= method0 :fdr) :bh method0))
  # validate the method BEFORE any
  # short-circuit - R errors on an unknown method at any length
  (unless (or (= method :none) (= method :bonferroni) (= method :holm)
              (= method :hochberg) (= method :bh) (= method :by))
    (errorf "stats: p-adjust unknown method %q" method))
  (unless (>= n lp)
    (errorf "stats: p-adjust needs n >= the number of p-values, got %q vs %d" n lp))
  (each p ps
    (unless (and (>= p 0) (<= p 1))
      (errorf "stats: p-adjust needs p in [0,1], got %q" p)))
  (cond
    (or (= method :none) (<= n 1)) ps
    (= method :bonferroni) (map |(min 1 (* n $)) ps)
    (do
      (def asc? (= method :holm))
      (def o (stable-order ps (not asc?)))
      (def q (if (= method :by) (sum (seq [i :range [1 (+ n 1)]] (/ 1 i))) 1))
      (def mults
        (case method
          :holm (seq [k :range [0 lp]] (- n k))
          :hochberg (seq [k :range [0 lp]] (+ 1 (- n lp) k))
          :bh (seq [k :range [0 lp]] (/ n (- lp k)))
          :by (seq [k :range [0 lp]] (/ (* q n) (- lp k)))
          (errorf "stats: p-adjust unknown method %q" method)))
      (def s-ps (map |(ps $) o))
      (def cums
        (if (= :num (num/resolve-impl opts))
          (let [k (num/kernels)
                prods ((k 'unpack-f64)
                        ((k 'vmap2) :mul (num/stage mults k) (num/stage s-ps k)))
                scan-op (if asc? :max :min)
                scanned ((k 'unpack-f64) ((k 'vscan) scan-op (num/stage prods k)))
                ones (num/stage (array/new-filled lp 1) k)]
            ((k 'unpack-f64) ((k 'vmap2) :min ones (num/stage scanned k))))
          (do
            (var acc (if asc? (- math/inf) math/inf))
            (def out @[])
            (each pr (seq [kk :range [0 lp]] (* (mults kk) (s-ps kk)))
              (if asc?
                (set acc (max acc pr))
                (set acc (min acc pr)))
              (array/push out (min 1 acc)))
            out)))
      (def out (array/new-filled lp 0))
      (eachp [kk idx] o
        (put out idx (cums kk)))
      out)))

(defn/typed binom-upper {:args [:number :number :number] :ret :number}
  [k n p]
  "P(K > k) = I_p(k+1, n-k) - R's lower.tail=FALSE convention
  (STRICTLY greater); the I_p form is far-tail accurate. Internal."
  (cond
    (>= k n) 0
    (< k 0) 1
    (pbeta p (+ k 1) (- n k) {})))

(defn binom-test
  "Exact binomial test, R binom.test: x successes of n trials vs
  p (default 0.5). Two-sided p = the small-probability-tail
  convention (relErr = 1 + 1e-7); Clopper-Pearson CI via qbeta.
  Returns {:statistic :parameter :p-value :conf-interval
  :estimate :null-value :alternative}. Opts: {:alternative
  :two-sided :conf-level 0.95 :impl}; throws on non-integer x/n,
  x > n, p outside [0,1]."
  [x n &opt p opts]
  # R binom.test: the exact binomial test (landed here because the
  # sign test needs it)
  (default opts {})
  (default p 0.5)
  (def alt (get opts :alternative :two-sided))
  (def conf (get opts :conf-level 0.95))
  (def xr (math/round x))
  (when (or (< x 0) (> (math/abs (- x xr)) 1e-7))
    (errorf "stats: binom-test x must be a nonnegative integer, got %q" x))
  (def nr (math/round n))
  (when (or (< n 1) (> (math/abs (- n nr)) 1e-7) (> xr nr))
    (errorf "stats: binom-test n must be a positive integer >= x, got %q" n))
  (when (or (< p 0) (> p 1))
    # the actionable value summary, no whole-handle dumps
    # (%q on a struct prints its whole body)
    (errorf "stats: binom-test p must be a number in [0,1], got a %s" (type p)))
  (when (or (<= conf 0) (>= conf 1))
    (errorf "stats: binom-test conf-level must be in (0,1), got %q" conf))
  (def m (* n p))
  (def pv
    (case alt
      :less (pbinom xr nr p {})
      :greater (binom-upper (- xr 1) nr p)
      :two-sided
      (cond
        (= p 0) (if (= xr 0) 1 0)
        (= p 1) (if (= xr nr) 1 0)
        (do
          (def relErr (+ 1 1e-7))
          (def d (dbinom xr nr p {}))
          (cond
            (= xr m) 1
            (< xr m)
            (do
              (def y (count |(<= (dbinom $ nr p {}) (* d relErr))
                             (range (math/ceil m) (+ nr 1))))
              (+ (pbinom xr nr p {}) (binom-upper (- nr y) nr p)))
            (do
              (def y (count |(<= (dbinom $ nr p {}) (* d relErr))
                             (range 0 (+ (math/floor m) 1))))
              (+ (if (pos? y) (pbinom (- y 1) nr p {}) 0)
                 (binom-upper (- xr 1) nr p))))))
      (errorf "stats: binom-test :alternative must be :two-sided/:less/:greater, got %q" alt)))
  (defn p-lo [alpha]
    (if (= xr 0) 0 (qbeta alpha xr (+ (- nr xr) 1) {})))
  (defn p-hi [alpha]
    (if (= xr nr) 1 (qbeta (- 1 alpha) (+ xr 1) (- nr xr) {})))
  (def ci
    (case alt
      :less @[0 (p-hi (- 1 conf))]
      :greater @[(p-lo (- 1 conf)) 1]
      (let [alpha (/ (- 1 conf) 2)]
        @[(p-lo alpha) (p-hi alpha)])))
  {:statistic xr :parameter nr :p-value pv
   :conf-interval ci :estimate (/ xr nr) :null-value p
   :alternative alt})

(defn poisson-test
  "Exact Poisson rate test, R poisson.test: x a count (or a 2-array
  for the rate-ratio comparison), t the time base (scalar or
  2-array, default 1), r the null rate (default 1). The SMR lane:
  pass t = the expected count, r = 1 (estimate = the SMR).
  Returns the R htest shape; the 2-count lane delegates to
  binom-test (R's own route) with the CI transformed. Opts:
  {:alternative :two-sided :conf-level 0.95 :impl}; slot 2 must be
  a scalar or 2-array. Upper tails via the Q(k, m) gamma identity
  (far-tail accurate)."
  [x &opt t r opts]
  (default opts {})
  (default r 1)
  # the t-slot trap, guarded (t is a
  # scalar or 2-array - anything else is caller misuse)
  (when (and t (not (number? t)) (not (indexed? t)))
    (error "stats: poisson-test expected a number or 2-array (the time base) in slot 2 - pass nil for the default 1, opts table trailing"))
  (def alt (get opts :alternative :two-sided))
  (def conf (get opts :conf-level 0.95))
  (def k (if (indexed? x) (length x) 1))
  (when (> k 2) (error "stats: poisson-test supports 1 or 2 counts"))
  (defn check-count [c]
    (def cr (math/round c))
    (when (or (< c 0) (> (math/abs (- c cr)) 1e-7))
      (errorf "stats: poisson-test x must be a nonnegative integer, got %q" c))
    cr)
  (when (or (< r 0) (not (number? r)))
    (errorf "stats: poisson-test r must be a nonnegative number, got %q" r))
  (if (= k 2)
    # the rate-ratio comparison: R delegates to binom.test at
    # p = r*T1/(r*T1+T2), CI transformed pp/(1-pp) * T2/T1
    (let [x1 (check-count (x 0))
          x2 (check-count (x 1))
          t1 (if (indexed? t) (t 0) t)
          t2 (if (indexed? t) (t 1) t)
          pnull (/ (* r t1) (+ (* r t1) t2))
          base (binom-test x1 (+ x1 x2) pnull opts)
          pp (base :conf-interval)]
      (merge base
             {:statistic x1
              :estimate (/ (/ x1 t1) (/ x2 t2))
              :conf-interval @[(* (/ (pp 0) (- 1 (pp 0))) (/ t2 t1))
                              (* (/ (pp 1) (- 1 (pp 1))) (/ t2 t1))]}))
    (let [xc (check-count (if (indexed? x) (x 0) x))
          tc (if (indexed? t) (t 0) t)
          mm (* r tc)
          pv
          (case alt
            :less (ppois xc mm {})
            :greater (gamma-p xc mm)
            (do
              (def relErr (+ 1 1e-7))
              (def d (dpois xc mm {}))
              (cond
                (= mm 0) (if (= xc 0) 1 0)
                (= xc mm) 1
                (< xc mm)
                (do
                  (var nn (math/ceil (- (* 2 mm) xc)))
                  (while (> (dpois nn mm {}) d)
                    (set nn (* 2 nn)))
                  (def y (count |(<= (dpois $ mm {}) (* d relErr))
                                 (range (math/ceil mm) (+ nn 1))))
                  (+ (ppois xc mm {}) (gamma-p (+ 1 (- nn y)) mm)))
                (do
                  (def y (count |(<= (dpois $ mm {}) (* d relErr))
                                 (range 0 (+ (math/floor mm) 1))))
                  (+ (if (pos? y) (ppois (- y 1) mm {}) 0) (gamma-p xc mm))))))]
      (defn p-lo [alpha] (if (= xc 0) 0 (qgamma alpha xc 1 {})))
      (defn p-hi [alpha] (qgamma (- 1 alpha) (+ xc 1) 1 {}))
      (def ci
        (case alt
          :less @[0 (/ (p-hi (- 1 conf)) tc)]
          :greater @[(/ (p-lo (- 1 conf)) tc) math/inf]
          (let [alpha (/ (- 1 conf) 2)]
            @[(/ (p-lo alpha) tc) (/ (p-hi alpha) tc)])))
      {:statistic xc :parameter tc :p-value pv
       :conf-interval ci :estimate (/ xc tc) :null-value r
       :alternative alt})))

(defn sign-test
  "Sign test: signs of xs - :mu (default 0) reduced to the exact
  binomial test at p = 0.5; zeros DROPPED (R's
  binom.test-on-signs convention), :n-zero carries the drop
  count. Opts: {:mu 0 + binom-test's}."
  [xs &opt opts]
  # signs of xs - mu (default 0) reduced to the exact binomial test
  # at p = 0.5; zeros DROPPED
  (default opts {})
  (def mu (get opts :mu 0))
  (arrays-only opts xs)
  (def pos (count |(> (- $ mu) 0) xs))
  (def nz (count |(not= (- $ mu) 0) xs))
  (when (zero? nz)
    (error "stats: sign-test needs at least one nonzero difference"))
  (merge (binom-test pos nz 0.5 opts)
         {:n-zero (- (length xs) nz)}))

(defn/typed sgn {:args [:number] :ret :number} [v]
  "sign(v): 1 / -1 / 0. Internal."
  (cond (pos? v) 1 (neg? v) -1 0))

(defn prop-test
  "One/two-proportion test, R prop.test. x/n scalar or 1-2 array,
  p scalar (k=1) or 2-array or nil (k=2 pooled equality).
  k=1: Wilson score CI + continuity correction (the YATES cap
  makes the CI null-DEPENDENT when |x - n*p| < 0.5); k=2 pooled:
  cc Wald CI on DELTA (unpooled SE); k=2 with p given: chi-square
  only - df 2, forced two-sided, NO CI (R verbatim). Two-sided p
  = the gamma-q upper chisq tail; one-sided = sign*sqrt(X2)
  through the exact erfc normal tail. Returns {:statistic
  :parameter :p-value :conf-interval :estimate :null-value
  :alternative}. Opts: {:alternative :two-sided :conf-level 0.95
  :correct true :impl}; non-integer x allowed (R verbatim, unlike
  binom-test); p in (0,1) EXCLUSIVE."
  [x n &opt p opts]
  # R prop.test: one/two-proportion tests
  (default opts {})
  (defn as-counts [v what]
    (cond
      (number? v)
      (do (unless (and (< v math/inf) (> v (- math/inf)))
            (errorf "stats: prop-test %s must be finite, got %q" what v))
          @[v])
      (indexed? v)
      (do (each e v
           (unless (and (number? e) (< e math/inf) (> e (- math/inf)))
             (errorf "stats: prop-test %s must be finite, got %q" what e)))
          v)
      (errorf "stats: prop-test %s must be numbers or arrays; unpack packed handles first, got %q" what v)))
  (def xs (as-counts x "x"))
  (def ns (as-counts n "n"))
  (def k (length xs))
  (unless (= k (length ns))
    (errorf "stats: prop-test x and n must have the same length, got %d vs %d" k (length ns)))
  (unless (<= k 2)
    (errorf "stats: prop-test supports 1 or 2 proportions, got %d" k))
  (for i 0 k
    (when (<= (ns i) 0) (errorf "stats: prop-test n must be positive, got %q" (ns i)))
    (when (< (xs i) 0) (errorf "stats: prop-test x must be nonnegative, got %q" (xs i)))
    (when (> (xs i) (ns i))
      (errorf "stats: prop-test x must not be greater than n, got %q > %q" (xs i) (ns i))))
  (def p-given (not (nil? p)))
  (def pv
    (cond
      (not p-given) (if (= k 1) @[0.5] nil)
      (number? p)
      (do (unless (= k 1)
            (errorf "stats: prop-test p must match x and n (scalar for one proportion)"))
          (unless (and (> p 0) (< p 1))
            (errorf "stats: prop-test p must be in (0,1), got %q" p))
          @[p])
      (indexed? p)
      (do (unless (= (length p) k)
            (errorf "stats: prop-test p must have the same length as x and n, got %d" (length p)))
          (each e p
            (unless (and (> e 0) (< e 1))
              (errorf "stats: prop-test p must be in (0,1), got %q" e)))
          p)
      (errorf "stats: prop-test p must be a number or array, got %q" p)))
  (def alt0 (get opts :alternative :two-sided))
  (unless (or (= alt0 :two-sided) (= alt0 :less) (= alt0 :greater))
    (errorf "stats: prop-test :alternative must be :two-sided/:less/:greater, got %q" alt0))
  (def alt (if (and p-given (= k 2)) :two-sided alt0))
  (def conf (get opts :conf-level 0.95))
  (unless (and (> conf 0) (< conf 1))
    (errorf "stats: prop-test conf-level must be in (0,1), got %q" conf))
  (def yates0 (if (get opts :correct true) 0.5 0))
  (def phat (map (fn [a b] (/ a b)) xs ns))
  (def zc (qnorm (if (= alt :two-sided) (/ (+ 1 conf) 2) conf) 0 1 {}))
  (def [yates ci]
    (if (= k 1)
      (let [p0 (pv 0) n0 (ns 0) ph (phat 0)
            y (min yates0 (math/abs (- (xs 0) (* n0 p0))))
            z22n (/ (* zc zc) (* 2 n0))
            pcu (+ ph (/ y n0))
            pcl (- ph (/ y n0))
            sq (fn [pc] (math/sqrt (+ (/ (* pc (- 1 pc)) n0)
                                     (/ z22n (* 2 n0)))))
            pu (if (>= pcu 1) 1
                 (/ (+ pcu z22n (* zc (sq pcu))) (+ 1 (* 2 z22n))))
            pl (if (<= pcl 0) 0
                 (/ (- (+ pcl z22n) (* zc (sq pcl))) (+ 1 (* 2 z22n))))]
        (case alt
          :two-sided @[y @[(max pl 0) (min pu 1)]]
          :greater @[y @[(max pl 0) 1]]
          @[y @[0 (min pu 1)]]))
      (if p-given
        @[yates0 nil]
        (let [delta (- (phat 0) (phat 1))
              y (min yates0 (math/abs (/ delta (+ (/ 1 (ns 0)) (/ 1 (ns 1))))))
              se (math/sqrt (+ (/ (* (phat 0) (- 1 (phat 0))) (ns 0))
                              (/ (* (phat 1) (- 1 (phat 1))) (ns 1))))
              width (+ (* zc se) (* y (+ (/ 1 (ns 0)) (/ 1 (ns 1)))))]
          (case alt
            :two-sided @[y @[(max (- delta width) -1) (min (+ delta width) 1)]]
            :greater @[y @[(max (- delta width) -1) 1]]
            @[y @[-1 (min (+ delta width) 1)]])))))
  # chi-square over the 2k cells (p given or pooled; the k=1 default
  # 0.5); df k-1 pooled / k given (k=1 is ALWAYS df 1: R assigns the
  # 0.5 null before the parameter branch)
  (def p-use
    (cond p-given pv (= k 1) pv (array/new-filled k (/ (sum xs) (sum ns)))))
  (var x2 0)
  (for i 0 k
    (def pp (p-use i))
    (def e1 (* (ns i) pp))
    (def e2 (* (ns i) (- 1 pp)))
    (def dd (- (math/abs (- (xs i) e1)) yates))
    (set x2 (+ x2 (/ (* dd dd) e1) (/ (* dd dd) e2))))
  (def df (if (or p-given (= k 1)) k (- k 1)))
  (def pval
    (if (= alt :two-sided)
      (gamma-q (/ df 2) (/ x2 2))
      (let [z (* (if (= k 1)
                  (sgn (- (phat 0) (pv 0)))
                  (sgn (- (phat 0) (phat 1))))
                (math/sqrt x2))]
        (if (= alt :less)
          (* 0.5 (erfc (/ (- z) (math/sqrt 2))))
          (* 0.5 (erfc (/ z (math/sqrt 2))))))))
  {:statistic x2 :parameter df :p-value pval :conf-interval ci
   :estimate (if (= k 1) (phat 0) @[(phat 0) (phat 1)])
   :null-value (if (= k 2) (if p-given @[(pv 0) (pv 1)] nil) (pv 0))
   :alternative alt})

#==------------------------------------------------------------------------==#
# The two CRAN-oracle entries (epiR 2.0.99, effectsize 1.0.3).
#==------------------------------------------------------------------------==#

(defn/typed ncp-solve {:args [:number :number :number] :ret :number}
  [t df p]
  "The ncp where P(T_{df,ncp} <= t) = p - the .get_ncp_t equation
  (effectsize solves it heuristically via optim; ours is the exact
  monotone inversion: pnt is DECREASING in ncp, so solve -f on the
  t+-40 bracket). Internal."
  (solve-monotone (fn [ncp] (- (pt t df {:ncp ncp}))) (- t 40) (+ t 40)
                   (- p)))

(defn one-sample-stats
  "[m s n] for the one-sample/paired shapes; resident-safe.
  Internal."
  [data &opt opts]
  (default opts {})
  (def n (num/count data))
  (def m (mean data opts))
  [m (sd data m opts) n])

(defn cohens-d-core
  "The effect-size engine: estimate + se + df + the exact
  noncentral-t CI (no Hedges-Olkin normal approximation).
  Moments via the descriptive kernels (the t-test resident pattern:
  paired diffs are the vmap2 :sub of staged inputs - never
  unpacked). Internal."
  [xs ys opts]
  # the .effect_size_difference decode
  (def mu (get opts :mu 0))
  (def conf (get opts :conf-level 0.95))
  (def alt (get opts :alternative :two-sided))
  (unless (or (= alt :two-sided) (= alt :less) (= alt :greater))
    (errorf "stats: cohens-d :alternative must be :two-sided/:less/:greater, got %q" alt))
  (unless (and (> conf 0) (< conf 1))
    (errorf "stats: cohens-d conf-level must be in (0,1), got %q" conf))
  (def paired (truthy? (get opts :paired)))
  (when (and paired (nil? ys))
    (error "stats: cohens-d :paired needs a second sample"))
  (num/check-impl-landed opts xs)
  (def [m0 s n df1 se hn pooled-sd]
    (cond
      (and paired (not (nil? ys)))
      (do
        (when (not= (num/count xs) (num/count ys))
          (error "stats: cohens-d paired samples must have equal length"))
        (def d (if (and (= :num (num/resolve-impl opts))
                        (not (or (number? xs) (number? ys))))
                 (let [k (num/kernels)]
                   ((k 'vmap2) :sub (num/stage xs k) (num/stage ys k)))
                 (seq [i :range [0 (length xs)]] (- (xs i) (ys i)))))
        (let [[m s n] (one-sample-stats d opts)]
          [m s n (- n 1) (/ s (math/sqrt n)) (/ n) nil]))
      (nil? ys)
      (let [[m s n] (one-sample-stats xs opts)]
        [m s n (- n 1) (/ s (math/sqrt n)) (/ n) nil])
      (do
        (def n1 (num/count xs))
        (def n2 (num/count ys))
        (when (or (< n1 2) (< n2 2))
          (errorf "stats: cohens-d needs >= 2 observations per sample, got %d/%d" n1 n2))
        (def mx (mean xs opts))
        (def my (mean ys opts))
        (def s1 (sd xs mx opts))
        (def s2 (sd ys my opts))
        (def pooled (truthy? (get opts :pooled-sd true)))
        (if pooled
          (let [sp (math/sqrt (/ (+ (* (- n1 1) (* s1 s1)) (* (- n2 1) (* s2 s2)))
                                 (- (+ n1 n2) 2)))
                hn (+ (/ n1) (/ n2))]
            [(- mx my) sp (+ n1 n2) (- (+ n1 n2) 2)
             (* sp (math/sqrt hn)) hn true])
          (let [sv (math/sqrt (/ (+ (* s1 s1) (* s2 s2)) 2))
                se1 (/ s1 (math/sqrt n1))
                se2 (/ s2 (math/sqrt n2))
                se (math/sqrt (+ (* se1 se1) (* se2 se2)))
                hn (/ (* 2 (+ (* n2 (* s1 s1)) (* n1 (* s2 s2))))
                      (* n1 n2 (+ (* s1 s1) (* s2 s2))))
                dfw (/ (* se se se se)
                       (+ (/ (* se1 se1 se1 se1) (- n1 1))
                          (/ (* se2 se2 se2 se2) (- n2 1))))]
            [(- mx my) sv (+ n1 n2) dfw se hn false])))))
  (when (< n 2) (errorf "stats: cohens-d needs >= 2 observations, got %d" n))
  (when (= s 0)
    (error "stats: cohens-d needs a nonzero sample sd (constant input?)"))
  (def d (/ (- m0 mu) s))
  (def t1 (/ (- m0 mu) se))
  # the ncp CI: level two-sided conf; one-sided 2*conf-1 (.adjust_ci)
  (def level (if (= alt :two-sided) conf (- (* 2 conf) 1)))
  (def alpha (/ (- 1 level) 2))
  (def sq-hn (math/sqrt hn))
  (def lo (* (ncp-solve t1 df1 (- 1 alpha)) sq-hn))
  (def hi (* (ncp-solve t1 df1 alpha) sq-hn))
  (def ci (case alt
            :two-sided @[lo hi]
            :less @[(- math/inf) hi]
            @[lo math/inf]))
  {:estimate d :conf-interval ci :alternative alt :conf-level conf
   :mu mu :paired paired :df df1 :se se :pooled-sd pooled-sd
   :t t1 :hn hn})

(defn cohens-d
  "Standardized mean difference, effectsize::cohens_d: two-sample
  (:pooled-sd true default / false = the sqrt((s1^2+s2^2)/2)
  variant), one-sample (ys nil), paired (:paired true). Returns
  {:estimate :conf-interval :alternative :conf-level :mu :paired
  :pooled-sd :df :se} - the CI is the exact noncentral-t interval.
  Opts: {:mu 0 :paired false :pooled-sd true :alternative
  :two-sided :conf-level 0.95 :impl}; CONSTANT input (zero sd) ->
  VALIDATION (R silently yields +-Inf/NA - our convention)."
  [xs &opt ys opts]
  # effectsize::cohens_d: standardized mean difference
  (default opts {})
  (def r (cohens-d-core xs ys opts))
  {:estimate (r :estimate) :conf-interval (r :conf-interval)
   :alternative (r :alternative) :conf-level (r :conf-level)
   :mu (r :mu) :paired (r :paired) :pooled-sd (r :pooled-sd)
   :df (r :df) :se (r :se)})

(defn hedges-g
  "Bias-corrected standardized mean difference,
  effectsize::hedges_g: cohens-d with the EXACT J bias correction
  (J multiplies the estimate AND both CI bounds). Returns the
  cohens-d shape plus :J. Opts: cohens-d's."
  [xs &opt ys opts]
  # effectsize::hedges_g
  (default opts {})
  (def r (cohens-d-core xs ys opts))
  (def df1 (r :df))
  (def J (math/exp (- (lgamma (/ df1 2))
                      (math/log (math/sqrt (/ df1 2)))
                      (lgamma (/ (- df1 1) 2)))))
  (def ci (r :conf-interval))
  {:estimate (* J (r :estimate))
   :conf-interval @[(* J (ci 0)) (* J (ci 1))]
   :alternative (r :alternative) :conf-level (r :conf-level)
   :mu (r :mu) :paired (r :paired) :pooled-sd (r :pooled-sd)
   :df df1 :se (r :se) :J J})

(defn/typed zz2stat {:args [:number :number :number :number :number] :ret :number}
  [p1x nx p1y ny dif]
  "epiR's zz2stat verbatim: the Miettinen-Nurminen constrained-MLE
  score statistic for H0: p1 - p2 = dif (the cubic acos/cos
  root). pow mirrors R's ^ (b3^3 vs triple-multiply can differ 1
  ulp). Internal."
  (def diff (- p1x p1y dif))
  (if (= diff 0)
    0
    (let [tt (/ ny nx)
          a (+ 1 tt)
          b (- (+ 1 tt p1x (* tt p1y) (* dif (+ tt 2))))
          c (+ (* dif dif) (* dif (+ (* 2 p1x) tt 1)) p1x (* tt p1y))
          d (- (* p1x dif (+ 1 dif)))
          b3 (/ b (* 3 a))
          v (+ (math/pow b3 3) (- (/ (* b c) (* 6 a a))) (/ d (* 2 a)))
          s (math/sqrt (- (math/pow b3 2) (/ c (* 3 a))))
          u (if (> v 0) s (- s))
          w (/ (+ 3.141592654 (math/acos (/ v (math/pow u 3)))) 3)
          p1d (- (* 2 u (math/cos w)) b3)
          p2d (- p1d dif)
          var (+ (/ (* p1d (- 1 p1d)) nx) (/ (* p2d (- 1 p2d)) ny))]
      (/ (* diff diff) var))))

(defn epi-2x2
  "Epidemiologic 2x2 measures, epiR epi.2by2: rows = EXPOSURE
  (exposed/unexposed), cols = OUTCOME (yes/no) - the as.columns
  layout. cohort-count: {:rr :or :rd :rd-score}, each {:estimate
  :conf-interval} - RR log/Katz, OR Woolf, RD Wald + the
  Miettinen-Nurminen score interval (the halving bisection
  mirrored verbatim); zero-cell estimates NaN the RR/OR triples
  (epiR's guard; RD stays finite); estimates are PROPORTIONS
  (epiR's units=100 is display scale). :case-control -> :or ONLY
  (the cohort measures are undefined - keys absent). Opts:
  {:method :cohort-count :conf-level 0.95 :impl}."
  [table &opt opts]
  # epiR epi.2by2 - rows = EXPOSURE (exposed/unexposed),
  # cols = OUTCOME (yes/no) - the as.columns layout
  (default opts {})
  (def method (get opts :method :cohort-count))
  (unless (or (= method :cohort-count) (= method :case-control))
    (errorf "stats: epi-2x2 :method must be :cohort-count/:case-control, got %q" method))
  (def conf (get opts :conf-level 0.95))
  (unless (and (> conf 0) (< conf 1))
    (errorf "stats: epi-2x2 conf-level must be in (0,1), got %q" conf))
  (unless (and (indexed? table) (= 2 (length table))
               (all indexed? table) (all |(= 2 (length $)) table))
    (errorf "stats: epi-2x2 needs a 2x2 table @[@[a b] @[c d]], got %q" table))
  (def row1 (table 0))
  (def row2 (table 1))
  (def a (row1 0))
  (def b (row1 1))
  (def c (row2 0))
  (def d (row2 1))
  (each cell [a b c d]
    (unless (and (number? cell) (< cell math/inf) (> cell (- math/inf)))
      (errorf "stats: epi-2x2 cells must be finite numbers, got %q" cell))
    (when (< cell 0)
      (errorf "stats: epi-2x2 cells must be nonnegative, got %q" cell)))
  (def zc (qnorm (- 1 (/ (- 1 conf) 2)) 0 1 {}))
  (defn nan? [v] (or (not= v v) (= v math/inf) (= v (- math/inf))))
  (defn log-ci [est lnvar]
    # the RR/OR lane: CI on the log scale; epiR's NaN guard - a
    # 0/NaN/Inf estimate NaNs the whole triple
    (if (or (nan? est) (= est 0))
      {:estimate math/nan :conf-interval @[math/nan math/nan]}
      (let [se (math/sqrt lnvar)
            lne (math/log est)]
        {:estimate est
         :conf-interval @[(math/exp (- lne (* zc se)))
                          (math/exp (+ lne (* zc se)))]})))
  (def n1 (+ a b))
  (def n0 (+ c d))
  (def rr
    (log-ci (/ (/ a n1) (/ c n0))
             (+ (/ a) (- (/ n1)) (/ c) (- (/ n0)))))
  (def orv
    (log-ci (/ (/ a b) (/ c d))
             (+ (/ a) (/ b) (/ c) (/ d))))
  (def rd
    (let [p1 (/ a n1)
          p0 (/ c n0)
          est (- p1 p0)
          se (math/sqrt (+ (/ (* a (- n1 a)) (* n1 n1 n1))
                           (/ (* c (- n0 c)) (* n0 n0 n0))))]
      {:estimate est
       :conf-interval @[(- est (* zc se)) (+ est (* zc se))]}))
  (def rd-score
    # zARscore verbatim: the 50-step dp-halving bisection on the MN
    # score statistic (accept when score < qchisq(conf, 1); stop at
    # dp < 1e-7 or |z - score| < 1e-6; the LAST TESTED point is the
    # bound - epiR's convention, mirrored exactly)
    (let [px (/ a n1)
          py (/ c n0)
          z2 (qchisq conf 1 {})
          hi
          (do
            (var proot (- px py))
            (var dp (- 1 proot))
            (var out nil)
            (var it 0)
            (while (<= it 50)
              (++ it)
              (set dp (* 0.5 dp))
              (def up2 (+ proot dp))
              (when (< (zz2stat px n1 py n0 up2) z2)
                (set proot up2))
              (when (or (< dp 1e-7)
                        (< (math/abs (- z2 (zz2stat px n1 py n0 up2))) 1e-6))
                (set out up2)
                (break)))
            out)
          lo
          (do
            (var proot (- px py))
            (var dp (+ 1 proot))
            (var out nil)
            (var it 0)
            (while (<= it 50)
              (++ it)
              (set dp (* 0.5 dp))
              (def low2 (- proot dp))
              (when (< (zz2stat px n1 py n0 low2) z2)
                (set proot low2))
              (when (or (< dp 1e-7)
                        (< (math/abs (- z2 (zz2stat px n1 py n0 low2))) 1e-6))
                (set out low2)
                (break)))
            out)]
      {:estimate (- px py)
       :conf-interval @[lo hi]}))
  (if (= method :cohort-count)
    {:method method :conf-level conf
     :rr rr :or orv :rd rd :rd-score rd-score}
    {:method method :conf-level conf :or orv}))

(defn cor-test
  "Correlation test with a significance test, R cor.test:
  :pearson (t + Fisher-z CI, n > 3, the |r| = 1 pole gives t
  +-inf / p 0 / CI [+-1, +-1]), :kendall (exact DP / S-approx,
  auto-approx at n >= 50), :spearman (brute n <= 9 / AS-89
  Edgeworth / t-approx on ties). Returns {:estimate :statistic
  :p-value :null-value} (+ :parameter :conf-interval on pearson).
  CONSTANT input mirrors R's NA as NaN triples. Opts: {:method
  :pearson :alternative :two-sided :conf-level 0.95 :exact
  :continuity false :impl}."
  [xs ys &opt opts]
  # R cor.test: pearson (t + Fisher-z CI), kendall
  # (exact DP / S-approx), spearman (brute n<=9 / AS-89 Edgeworth /
  # t-approx). Estimates via the landed correlations (pearson
  # kernel-RESIDENT; spearman/kendall sort-bound composition - the
  # arrays-only guard delegated). Degenerate (constant) input
  # mirrors R's NA: NaN estimate/statistic/p.
  (default opts {})
  (def method (get opts :method :pearson))
  (unless (or (= method :pearson) (= method :spearman) (= method :kendall))
    (errorf "stats: cor-test :method must be :pearson/:spearman/:kendall, got %q" method))
  (def alt (get opts :alternative :two-sided))
  (unless (or (= alt :two-sided) (= alt :less) (= alt :greater))
    (errorf "stats: cor-test :alternative must be :two-sided/:less/:greater, got %q" alt))
  (def conf (get opts :conf-level 0.95))
  (unless (and (>= conf 0) (<= conf 1))
    (errorf "stats: cor-test conf-level must be in [0,1], got %q" conf))
  (def cont (truthy? (get opts :continuity)))
  (num/check-impl-landed opts xs)
  (def n (num/count xs))
  (unless (= n (num/count ys))
    (errorf "stats: cor-test x and y must have the same length, got %d vs %d"
        n (num/count ys)))
  # the non-finite VALIDATION only on janet-side arrays (packed
  # handles chain - the NaN-domain contract)
  (each v [xs ys]
    (when (indexed? v)
  (each x v
    (unless (and (< x math/inf) (> x (- math/inf)))
      (error "stats: cor-test needs finite observations")))))
  (defn const? [v]
    (and (indexed? v) (all |(= $ (v 0)) v)))
  (def degenerate (or (const? xs) (const? ys)))
  (case method
    :pearson
    (do
  (when (< n 3)
    (errorf "stats: cor-test pearson needs >= 3 observations, got %d" n))
  (if degenerate
    {:method method :alternative alt :estimate math/nan
     :statistic math/nan :p-value math/nan}
    (let [r (pearson xs ys opts)
          df (- n 2)
          t (if (= 1 (math/abs r))
              (* (sgn r) math/inf)
              (/ (* (math/sqrt df) r) (math/sqrt (- 1 (* r r)))))
          # the two tails via the landed pt + the direct
          # beta-inc upper (far-tail accurate both sides)
          lower (pt t df {})
          upper (if (pos? t)
                  (* 0.5 (beta-inc (/ df 2) 0.5 (/ df (+ df (* t t)))))
                  (pt (- t) df {}))
          pv (case alt
               :less lower
               :greater upper
               (* 2 (min lower upper)))
          ci (when (> n 3)
               (let [z (math/atanh r)
                     sigma (/ 1 (math/sqrt (- n 3)))
                     zc (qnorm (if (= alt :two-sided)
                                 (/ (+ 1 conf) 2)
                                 conf)
                               0 1 {})]
                 (case alt
                   :two-sided @[(math/tanh (- z (* zc sigma)))
                                (math/tanh (+ z (* zc sigma)))]
                   # R tanh's the WHOLE cint: tanh(-Inf) = -1
                   :less @[-1 (math/tanh (+ z (* zc sigma)))]
                   @[(math/tanh (- z (* zc sigma))) 1])))]
      {:method method :alternative alt :estimate r
       :statistic t :parameter df :p-value pv
       :conf-interval ci :null-value 0})))
    (do
  (when (< n 2)
    (errorf "stats: cor-test needs >= 2 observations, got %d" n))
  (if degenerate
    {:method method :alternative alt :estimate math/nan
     :statistic math/nan :p-value math/nan}
    (do
      (unless (and (indexed? xs) (indexed? ys))
        (error "stats: composition path consumes janet arrays; unpack packed handles first"))
      (let [ties? (< (min (length (distinct xs)) (length (distinct ys))) n)]
      (if (= method :kendall)
        (let [r (kendall xs ys opts)
              exact-flag (if (nil? (get opts :exact)) (< n 50)
                           (truthy? (get opts :exact)))
              half (/ (* n (- n 1)) 4)]
          (if (and exact-flag (not ties?))
            (let [q (math/round (* (+ r 1) half))
                  p-lo (kendall-p q n)
                  p-hi (- 1 (kendall-p (- q 1) n))
                  pv (case alt
                       :less p-lo
                       :greater p-hi
                       (let [p (if (> q half) p-hi p-lo)]
                         (min (* 2 p) 1)))]
              {:method method :alternative alt :estimate r
               :statistic q :p-value pv :null-value 0})
            # the S-approx lane (ties, n >= 50, or :exact false)
            (let [xt (rank-multiplicities (sorted-copy xs))
                  yt (rank-multiplicities (sorted-copy ys))
                  T0 (/ (* n (- n 1)) 2)
                  sum-tt (fn [ts] (var s 0) (each t ts (+= s (* t (- t 1)))) s)
                  T1 (* 0.5 (sum-tt xt))
                  T2 (* 0.5 (sum-tt yt))
                  S (* r (math/sqrt (* (- T0 T1) (- T0 T2))))
                  v0 (* n (- n 1) (+ (* 2 n) 5))
                  st (fn [ts f]
                       (var s 0)
                       (each t ts (+= s (f t)))
                       s)
                  vt (st xt |(* $ (- $ 1) (+ (* 2 $) 5)))
                  vu (st yt |(* $ (- $ 1) (+ (* 2 $) 5)))
                  v1 (* (st xt |(* $ (- $ 1)))
                        (st yt |(* $ (- $ 1))))
                  v2 (* (st xt |(* $ (- $ 1) (- $ 2)))
                        (st yt |(* $ (- $ 1) (- $ 2))))
                  var-S (+ (/ (- v0 vt vu) 18)
                           (/ v1 (* 2 n (- n 1)))
                           (if (> n 2)
                             (/ v2 (* 9 n (- n 1) (- n 2)))
                             math/nan))
                  Sc (if cont (* (sgn S) (- (math/abs S) 1)) S)
                  z (/ Sc (math/sqrt var-S))
                  pv (case alt
                       :less (* 0.5 (erfc (/ (- z) (math/sqrt 2))))
                       :greater (* 0.5 (erfc (/ z (math/sqrt 2))))
                       (min (* 2 (min (* 0.5 (erfc (/ (- z) (math/sqrt 2))))
                                      (* 0.5 (erfc (/ z (math/sqrt 2))))))
                            1))]
              {:method method :alternative alt :estimate r
               :statistic z :p-value pv :null-value 0})))
        # spearman
        (let [r (spearman xs ys opts)
              exact0 (if (nil? (get opts :exact)) true
                       (truthy? (get opts :exact)))
              exact (and exact0 (not ties?))
              q (/ (* (- (* n (* n n)) n) (- 1 r)) 6)
              half (/ (- (* n (* n n)) n) 6)]
          (if (and exact (<= n 1290))
            (let [ps (fn [lower]
                       (spearman-p (+ (math/round q)
                                      (if lower 2 0))
                                   n lower))
                  pv (case alt
                       :two-sided (let [p (if (> q half) (ps false) (ps true))]
                                    (min (* 2 p) 1))
                       :greater (ps true)
                       :less (ps false))]
              {:method method :alternative alt :estimate r
               :statistic q :p-value pv :null-value 0})
            # the t-approx (n > 1290, :exact false, or ties)
            (let [den (+ (/ (* n (- (* n n) 1)) 6)
                        (if cont 1 0))
                  r2 (- 1 (/ q den))
                  pv (case alt
                       :two-sided (let [t (/ r2 (math/sqrt (/ (- 1 (* r2 r2)) (- n 2))))
                                        p-lo (pt t (- n 2) {})]
                                    (min (* 2 (min p-lo (- 1 p-lo))) 1))
                       :greater (pt (- (/ r2 (math/sqrt (/ (- 1 (* r2 r2)) (- n 2)))))
                                    (- n 2) {})
                       (pt (/ r2 (math/sqrt (/ (- 1 (* r2 r2)) (- n 2))))
                           (- n 2) {}))]
              {:method method :alternative alt :estimate r
               :statistic q :p-value pv :null-value 0}))))))))))

#==------------------------------------------------------------------------==#
# McNemar + Mantel-Haenszel (base oracles + epiR's stratified
# extensions).
#==------------------------------------------------------------------------==#

(defn table-2x2
  "Shared 2x2 validation: @[a b] / @[c d], nonnegative finite
  cells. Internal."
  [table what]
  # shared validation: @[a b] / @[c d], nonnegative finite cells
  (unless (and (indexed? table) (= 2 (length table))
               (all indexed? table) (all |(= 2 (length $)) table))
    (errorf "stats: %s needs a 2x2 table @[@[a b] @[c d]], got %q" what table))
  (each row table
    (each cell row
      (unless (and (number? cell) (< cell math/inf) (> cell (- math/inf)))
        (errorf "stats: %s cells must be finite numbers, got %q" what cell))
      (when (< cell 0)
        (errorf "stats: %s cells must be nonnegative, got %q" what cell))))
  table)

(defn mcnemar-test
  "McNemar's chi-squared for the paired-proportions 2x2, R
  mcnemar.test: the correction applies ONLY when the discordant
  cells differ (R's any(x - t(x) != 0) gate); zero discordant
  pairs -> NaN (R's 0/0; janet division throws - the explicit
  guard). Returns {:statistic :parameter :p-value}. Opts:
  {:correct true :impl}."
  [table &opt opts]
  # R mcnemar.test: the paired-proportions 2x2
  (default opts {})
  (table-2x2 table "mcnemar-test")
  (def b ((table 0) 1))
  (def c ((table 1) 0))
  (def denom (+ b c))
  (def cc (if (and (truthy? (get opts :correct true)) (not= b c)) 1 0))
  (if (zero? denom)
    {:statistic math/nan :parameter 1 :p-value math/nan}
    (let [diff (- (math/abs (- b c)) cc)
          x2 (/ (* diff diff) denom)]
      {:statistic x2 :parameter 1 :p-value (gamma-q 0.5 (/ x2 2))})))

(defn mantelhaen-test
  "Stratified 2x2xK Mantel-Haenszel test, R mantelhaen.test:
  the MH chi-square (the 0.5 cc cap only when |DELTA| >= 0.5),
  common odds ratio + the Robins-Breslow-Greenland CI, one-sided
  CIs [0, hi]/[lo, inf] - PLUS the epiR-provenance keys: :rr
  (Greenland-Robins pooled RR + log CI) and :homogeneity
  ({:rr {:statistic :p-value} :or {...}} - the inverse-variance
  weighted log-scale Q tests, df K-1; the OR lane via the
  Haldane-Anscombe +0.5 correction; nil at K = 1). tables = K
  stratum tables, rows-major @[a b]/@[c d]. Opts: {:alternative
  :two-sided :conf-level 0.95 :correct true :impl}."
  [tables &opt opts]
  # R mantelhaen.test: the stratified 2x2xK approx lane - MH
  # chi-square (+ the |DELTA| >= 0.5 cc cap), common odds ratio +
  # the Robins-Breslow-Greenland CI - plus the epiR-provenance
  # extensions: the Greenland-Robins pooled RR and the
  # inverse-variance homogeneity tests (the RR Katz-weighted; the
  # OR via the Haldane-Anscombe +0.5)
  (default opts {})
  (unless (and (indexed? tables) (>= (length tables) 1))
    (errorf "stats: mantelhaen-test needs >= 1 stratum table, got %q" tables))
  (def alt (get opts :alternative :two-sided))
  (unless (or (= alt :two-sided) (= alt :less) (= alt :greater))
    (errorf "stats: mantelhaen-test :alternative must be :two-sided/:less/:greater, got %q" alt))
  (def conf (get opts :conf-level 0.95))
  (unless (and (> conf 0) (< conf 1))
    (errorf "stats: mantelhaen-test conf-level must be in (0,1), got %q" conf))
  (def K (length tables))
  (each t tables (table-2x2 t "mantelhaen-test"))
  (def cells
    (seq [i :range [0 K]]
      (def t (tables i))
      (def a ((t 0) 0))
      (def b ((t 0) 1))
      (def c ((t 1) 0))
      (def d ((t 1) 1))
      (def n (+ a b c d))
      (when (< n 2)
        (errorf "stats: mantelhaen-test sample size in each stratum must be > 1, got %d" n))
      [a b c d n]))
  (def as (map |($ 0) cells))
  (def bs (map |($ 1) cells))
  (def cs (map |($ 2) cells))
  (def ds (map |($ 3) cells))
  (def ns (map |($ 4) cells))
  (def sxs (map + as bs))
  (def sys (map + as cs))
  (def DELTA (sum (seq [i :range [0 K]]
                    (- (as i) (* (sxs i) (sys i) (/ (ns i)))))))
  (def yates (if (and (truthy? (get opts :correct true))
                      (>= (math/abs DELTA) 0.5))
               0.5 0))
  (var vv 0)
  (for i 0 K
    (def n (ns i))
    (+= vv (/ (* (sxs i) (sys i) (- n (sxs i)) (- n (sys i)))
              (* n n (- n 1)))))
  (def num (* (- (math/abs DELTA) yates) (- (math/abs DELTA) yates)))
  (def x2 (if (zero? vv)
            (if (zero? num) math/nan math/inf)
            (/ num vv)))
  (def pval
    (if (= alt :two-sided)
      (gamma-q 0.5 (/ x2 2))
      (let [z (* (sgn DELTA) (math/sqrt x2))]
        (if (= alt :less)
          (* 0.5 (erfc (/ (- z) (math/sqrt 2))))
          (* 0.5 (erfc (/ z (math/sqrt 2))))))))
  # the common OR + the RBG sd (uncorrected cell sums throughout)
  (def sad (sum (seq [i :range [0 K]] (* (as i) (ds i) (/ (ns i))))))
  (def sbc (sum (seq [i :range [0 K]] (* (bs i) (cs i) (/ (ns i))))))
  (def or-mh (if (zero? sbc)
               (if (zero? sad) math/nan math/inf)
               (/ sad sbc)))
  (defn s2 [f]
    (sum (seq [i :range [0 K]]
          (f (as i) (bs i) (cs i) (ds i) (ns i)))))
  (def sd
    (if (or (zero? sad) (zero? sbc))
      math/nan
      (math/sqrt
        (+ (/ (s2 (fn [a b c d n] (* (+ a d) a d (/ (* n n)))))
              (* 2 sad sad))
           (/ (s2 (fn [a b c d n]
                     (/ (+ (* (+ a d) b c) (* (+ b c) a d)) (* n n))))
              (* 2 sad sbc))
           (/ (s2 (fn [a b c d n] (* (+ b c) b c (/ (* n n)))))
              (* 2 sbc sbc))))))
  (defn zc [pp] (qnorm pp 0 1 {}))
  (def ci
    (if (not= sd sd)
      @[math/nan math/nan]
      (case alt
        :two-sided (let [q (* (zc (/ (- 1 conf) 2)) sd)]
                     @[(math/exp (+ (math/log or-mh) q))
                      (math/exp (- (math/log or-mh) q))])
        :less @[0 (math/exp (+ (math/log or-mh) (* (zc conf) sd)))]
        @[(math/exp (- (math/log or-mh) (* (zc conf) sd))) math/inf])))
  # the Greenland-Robins pooled RR (epiR's sRR)
  (def m1s (map + as cs))
  (def sar (sum (seq [i :range [0 K]] (* (as i) (+ (cs i) (ds i)) (/ (ns i))))))
  (def scr (sum (seq [i :range [0 K]] (* (cs i) (+ (as i) (bs i)) (/ (ns i))))))
  (def rr-est (if (zero? scr)
                (if (zero? sar) math/nan math/inf)
                (/ sar scr)))
  (def rr-var
    (if (or (zero? sar) (zero? scr))
      math/nan
      (/ (sum (seq [i :range [0 K]]
                 (- (* (m1s i) (+ (as i) (bs i)) (+ (cs i) (ds i))
                         (/ (* (ns i) (ns i))))
                    (* (as i) (cs i) (/ (ns i))))))
         (* sar scr))))
  (def rr-ci
    (if (or (not= rr-var rr-var) (not= rr-est rr-est))
      @[math/nan math/nan]
      (let [w (* (zc (/ (+ 1 conf) 2)) (math/sqrt rr-var))]
        @[(math/exp (- (math/log rr-est) w))
          (math/exp (+ (math/log rr-est) w))])))
  # the homogeneity (heterogeneity) tests - nil at K = 1 (epiR
  # computes them only in the stratified branch)
  (def homog
    (when (>= K 2)
      (defn q-stat [est-f var-f]
        (def lns (map est-f cells))
        (def ws (map (fn [cl] (/ 1 (var-f cl))) cells))
        (def sw (sum ws))
        (def lns-w (/ (sum (map * ws lns)) sw))
        (var q 0)
        (for i 0 K
          (+= q (* (ws i) (* (- (lns i) lns-w) (- (lns i) lns-w)))))
        {:statistic q :p-value (gamma-q (/ (- K 1) 2) (/ q 2))})
      {:rr (q-stat (fn [cl]
                     (math/log (/ (/ (cl 0) (+ (cl 0) (cl 1)))
                                  (/ (cl 2) (+ (cl 2) (cl 3))))))
                   (fn [cl]
                     (+ (/ (cl 1) (* (cl 0) (+ (cl 0) (cl 1))))
                        (/ (cl 3) (* (cl 2) (+ (cl 2) (cl 3)))))))
       :or (q-stat (fn [cl]
                     (math/log (/ (* (+ (cl 0) 0.5) (+ (cl 3) 0.5))
                                  (* (+ (cl 1) 0.5) (+ (cl 2) 0.5)))))
                   (fn [cl]
                     (+ (/ 1 (+ (cl 0) 0.5)) (/ 1 (+ (cl 1) 0.5))
                        (/ 1 (+ (cl 2) 0.5)) (/ 1 (+ (cl 3) 0.5)))))}))
  {:statistic x2 :parameter 1 :p-value pval
   :estimate or-mh :conf-interval ci :alternative alt :null-value 1
   :rr {:estimate rr-est :conf-interval rr-ci}
   :homogeneity homog})

(defn oneway-test
  "One-way ANOVA test, R oneway.test verbatim: {:var-equal
  false} (the DEFAULT) = the WELCH lane (R's own default);
  {:var-equal true} = the classic F. Returns {:statistic
  :parameter [num-df den-df] :p-value :var-equal}; the classic
  lane's F equals the aov table's. Opts: {:var-equal false
  :impl}; VALIDATION before the branch, both lanes: k >= 2,
  n_i >= 2, finite observations (R's messages)."
  [groups &opt opts]
  # R oneway.test verbatim - default :var-equal false = the WELCH
  # lane (R's own default). VALIDATION before the branch, both
  # lanes (R's order): k >= 2, n_i >= 2. Moments via the
  # descriptive kernels (the t-test class); small-k glue is janet.
  (default opts {})
  (num/check-impl-landed opts)
  (def k (length groups))
  (when (< k 2) (error "stats: oneway-test needs >= 2 groups (R: not enough groups)"))
  (each g groups
   (unless (indexed? g)
     (error "stats: oneway-test consumes janet arrays per group; unpack packed handles first"))
   (when (< (length g) 2)
     (error "stats: oneway-test needs >= 2 observations per group (R: not enough observations)"))
   (each v g
     (unless (and (number? v) (< (math/abs v) math/inf))
       (error "stats: oneway-test needs finite observations"))))
  (def ns (seq [g :in groups] (length g)))
  (def ms (seq [g :in groups] (mean g opts)))
  (def vs (seq [[g m] :in (map (fn [a b] [a b]) groups ms)] (variance g m opts)))
  (if (get opts :var-equal false)
   (do
     # classic: F = [sum n_i (m_i - grand)^2 / (k-1)] /
     # [sum (n_i-1) v_i / (n-k)], grand = the PLAIN mean of all y
     (def ys (array ;(mapcat identity groups)))
     (def ybar (mean ys opts))
     (def n (reduce + 0 ns))
     (var ssb 0)
     (var ssw 0)
     (for i 0 k
       (def d (- (ms i) ybar))
       (+= ssb (* (ns i) (* d d)))
       (+= ssw (* (- (ns i) 1) (vs i))))
     # ssw = 0 (all groups constant): R's Inf/NaN - the explicit
     # guards (janet division throws)
     (def f (cond
              (zero? ssw) (if (zero? ssb) math/nan math/inf)
              (zero? ssb) 0
              (/ (/ ssb (- k 1)) (/ ssw (- n k)))))
     (def df-g (- k 1))
     (def df-r (- n k))
     {:statistic f :parameter @[df-g df-r]
      :p-value (- 1 (pf f df-g df-r {})) :var-equal true})
   (do
     # Welch: w_i = n_i/v_i (v_i = 0 -> inf, R's value), tmp from
     # the (1 - w_i/W)^2 terms, the WEIGHTED grand mean m
     (def ws (seq [i :range [0 k]]
               (if (zero? (vs i)) math/inf (/ (ns i) (vs i)))))
     (def wsum (reduce + 0 ws))
     (var tmp 0)
     (for i 0 k
       (def r (- 1 (/ (ws i) wsum)))
       (+= tmp (/ (* r r) (- (ns i) 1))))
     (set tmp (/ tmp (- (* k k) 1)))
     (var wm 0)
     (for i 0 k (+= wm (* (ws i) (ms i))))
     # the weighted grand mean m = sum(w m)/W - R sums first, then
     # the squared deviations use the SAME m
     (def m (/ wm wsum))
     (var num 0)
     (for i 0 k
       (def d (- (ms i) m))
       (+= num (* (ws i) (* d d))))
     (def den (+ 1 (* 2 (- k 2) tmp)))
     (def f (/ (/ num (- k 1)) den))
     (def df2 (/ 1 (* 3 tmp)))
     {:statistic f :parameter @[(- k 1) df2]
      :p-value (- 1 (pf f (- k 1) df2 {})) :var-equal false})))


