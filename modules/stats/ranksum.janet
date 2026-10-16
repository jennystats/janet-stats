# stats/ranksum.janet - exact rank-sum DPs (R psignrank /
# pwilcox conventions) + shared tie helpers + the correlation exact
# DPs (kendall.c / prho.c mirrors). Pure janet.
(import typed :prefix "")
(import ./utils :prefix "")

(defn/typed rank-multiplicities {:args [:array] :ret :array} [sorted-vals]
  "Sizes of tie groups in an ALREADY-SORTED value array."
  (def ts @[])
  (var i 0)
  (def n (length sorted-vals))
  (while (< i n)
    (def v (get sorted-vals i))
    (var j i)
    (while (and (< j n) (= (get sorted-vals j) v)) (++ j))
    (array/push ts (- j i))
    (set i j))
  ts)

(defn/typed sum-t3-minus-t {:args [:array] :ret :number} [ts]
  "Sum of (t^3 - t) over tie-group sizes - the shared tie
  correction for rank statistics."
  (var s 0)
  (each t ts (+= s (- (* t t t) t)))
  s)

(defn/typed signrank-p {:args [:number :number :boolean] :ret :number} [q n lower]
  "Exact P(V <= q) (lower) / P(V > q) (upper) for the signed-rank
  sum V over ranks {1..n} (R psignrank): DP over the 2^n sign
  assignments, uniform."
  (def total (math/pow 2 n))
  (var counts @{})
  (put counts 0 1)
  (for r 1 (+ n 1)
    (def nc @{})
    (eachp [s c] counts
      (put nc s (+ (get nc s 0) c))          # rank r assigned a minus
      (put nc (+ s r) (+ (get nc (+ s r) 0) c))) # rank r assigned a plus
    (set counts nc))
  (var acc 0)
  (eachp [s c] counts
    (when (if lower (<= s q) (> s q)) (+= acc c)))
  (/ acc total))

(defn/typed wilcox-p {:args [:number :number :number :boolean] :ret :number} [q nx ny lower]
  "Exact P(W <= q) (lower) / P(W > q) (upper) for the Mann-Whitney
  W = rank-sum of x - nx(nx+1)/2 (R pwilcox): DP over the equally
  likely size-nx subsets of {1..n}, C(n, nx) total."
  (def n (+ nx ny))
  # C(n, nx) - equally likely size-nx subsets
  (var prod 1)
  (for k (+ (- n nx) 1) (+ n 1) (*= prod k))
  (var nx-fact 1)
  (for k 1 (+ nx 1) (*= nx-fact k))
  (def total (/ prod nx-fact))
  # cnt[s][k] = subsets of {1..i} of size k summing to s - rolling
  (var cnt @{}) # cnt[k] -> table of sum -> count
  (put cnt 0 @{0 1})
  (for i 1 (+ n 1)
    (def nc @{})
    (eachp [k subs] cnt
      (eachp [s c] subs
        (put-in nc [k s] (+ (get-in nc [k s] 0) c))      # skip rank i
        (when (< k nx)
          (put-in nc [(+ k 1) (+ s i)] (+ (get-in nc [(+ k 1) (+ s i)] 0) c)))))
    (set cnt nc))
  (var acc 0)
  (eachp [s c] (get cnt nx)
    (def w (- s (* 0.5 nx (+ nx 1))))
    (when (if lower (<= w q) (> w q)) (+= acc c)))
  (/ acc total))

#==------------------------------------------------------------------------==#
# The correlation exact DPs - R's kendall.c and
# prho.c (AS 89), decoded verbatim from the r-source mirror.
#==------------------------------------------------------------------------==#

(def ck-table @{})

(defn/typed ckendall {:args [:number :number] :ret :number} [k n]
  "R's ckendall: number of permutations of n with exactly k
  concordant pairs; c(k,n) = sum_{i<n} c(k-i, n-1), c(k,1) =
  [k=0]; memoized over (n,k); counts exact in doubles to 2^53
  (R identical)."
  (def u (/ (* n (- n 1)) 2))
  (cond
    (or (< k 0) (> k u)) 0
    (= n 1) (if (= k 0) 1 0)
    (do
      (def row (or (get ck-table n)
                   (let [r (array/new-filled (+ u 1) -1)]
                     (put ck-table n r)
                     r)))
      (def known (row k))
      (if (>= known 0)
        known
        (do
          (var v 0)
          (for i 0 n (+= v (ckendall (- k i) (- n 1))))
          (put row k v)
          v)))))

(defn/typed kendall-p {:args [:number :number] :ret :number} [q n]
  "R's pkendall: P(T <= q), T the concordant-pair count; q floors
  at +1e-7 (R's round); n! via the integer product (doubles past
  22! - the same rounding class as R's gammafn)."
  (def qq (math/floor (+ q 1e-7)))
  (def u (/ (* n (- n 1)) 2))
  (cond
    (< qq 0) 0
    (> qq u) 1
    (do
      (var s 0)
      (for j 0 (+ qq 1) (+= s (ckendall j n)))
      (var f 1)
      (for i 1 (+ n 1) (set f (* f i)))
      (/ s f))))

(def spearman-hist @{})

(defn/typed spearman-histogram {:args [:number] :ret :table} [n]
  "ONE pass over all n! permutations (R's rotation enumeration,
  mirrored verbatim) building ise = sum(i+1 - l[i])^2 -> count;
  serves every q at this n. Memoized per n."
  (or (get spearman-hist n)
      (let [l @[]
            hist @{}]
        (for i 1 (+ n 1) (array/push l i))
        (def perms (do (var f 1)
                       (for i 1 (+ n 1) (set f (* f i)))
                       f))
        (for m 0 perms
          (var ise 0)
          (for i 0 n
            (def n1 (- (+ i 1) (l i)))
            (+= ise (* n1 n1)))
          (put hist ise (+ 1 (get hist ise 0)))
          # the next permutation by cyclic rotation (prho.c verbatim)
          (var n1 n)
          (var mt 0)
          (var go true)
          (while go
            (set mt (l 0))
            (for i 1 n1 (put l (- i 1) (l i)))
            (-- n1)
            (put l n1 mt)
            (set go (and (= mt (+ n1 1)) (> n1 1)))))
        (put spearman-hist n hist)
        hist)))

(defn/typed spearman-p {:args [:number :number :boolean] :ret :number}
  [is n lower]
  "R's prho (AS 89) verbatim: P(S < is) (lower) / P(S >= is);
  brute force for n <= 9, the Edgeworth series for
  10 <= n <= 1290 (R's 'exact' IS the Edgeworth series there).
  prho's is == n3 shortcut reproduces naturally via the histogram."
  (when (<= n 1) (errorf "stats: spearman-p needs n >= 2, got %q" n))
  (def n3 (/ (* n (- (* n n) 1)) 3))
  (defn clamp [v] (min 1 (max 0 v)))
  (cond
    (<= is 0) (if lower 0 1)
    (> is n3) (if lower 1 0)
    (<= n 9)
    (let [h (spearman-histogram n)
          nfac (do (var f 1)
                   (for i 1 (+ n 1) (set f (* f i)))
                   f)]
      (var ifr 0)
      (eachp [ise c] h
        (when (<= is ise) (+= ifr c)))
      (if lower (/ (- nfac ifr) nfac) (/ ifr nfac)))
    # the AS 89 Edgeworth series (10 <= n; R's 'exact' is THIS)
    (let [b (/ 1 n)
          x (* (- (/ (* 6 (- is 1) b) (- (* n n) 1)) 1)
               (math/sqrt (- n 1)))
          y (* x x)
          dd (+ 0.0072 (- (* 0.0831 b))
                (* y b (- 0.0131 (* 0.00046 y))))
          cc (- (+ 0.0879 (* 0.0151 b)) (* y dd))
          bb (- (+ -0.0758 (* b (+ 0.1033 (* 0.3932 b)))) (* y b cc))
          aa (+ 0.2274 (* b (+ 0.2531 (* 0.1745 b))))
          u (* x b (+ aa (* y bb)))
          w (/ u (math/exp (* 0.5 y)))
          phi (if lower
                (* 0.5 (erfc (/ (- x) (math/sqrt 2))))
                (* 0.5 (erfc (/ x (math/sqrt 2)))))]
      (clamp (if lower (- phi w) (+ phi w))))))

