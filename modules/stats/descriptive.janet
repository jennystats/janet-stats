# stats/descriptive.janet - moments, quantiles (R type-7),
# correlations, robust. {:impl} plumbing throughout; the moment
# family (mean/var/sd/skew/kurt) is kernel-backed (resident-packed:
# array in stages, handle in consumes - stats/num owns the
# resolution contract). pearson is kernel-backed (the five raw sums
# via vdot/vreduce, both input modes); quantile/median/iqr/
# spearman/kendall/mad/trimmed-mean are SORT-BOUND: no kernel path
# exists for a janet sort or a pairwise sign loop - :num runs the
# SAME documented janet composition (identical-formula policy);
# packed handles are rejected with an actionable unpack-first error.
(import typed :prefix "")
(import ./utils :prefix "")
(import ./num)

(defn arrays-only
  "The sort-bound-composition guard: resolves/validates the impl,
  then asserts janet-array inputs - a packed handle cannot feed a
  pure loop or a janet sort; it needs the num path (where a kernel
  path exists) or an unpack first."
  [opts & args]
  (num/check-impl-landed opts)
  (each x args
    (unless (indexed? x)
      (error "stats: this path consumes janet arrays; a packed handle needs {:impl :num} (where a kernel path exists) or an unpack first")))
  nil)

(defn mean
  "Arithmetic mean. Accepts a janet array or a packed handle under
  :num. Throws on an empty input. Opts: {:impl}."
  [xs &opt opts]
  (default opts {})
  (def n (num/count xs))
  (when (zero? n) (error "stats: mean of empty"))
  (if (= :num (num/resolve-impl opts))
    (let [k (num/kernels)]
      (/ ((k 'vreduce) :sum (num/stage xs k)) n))
    (do
      (num/check-impl-landed opts xs)
      (var s 0)
      (each x xs (+= s x))
      (/ s n))))

(defn variance
  "Sample variance (n-1 denominator, R var default). Slot 2 is an
  optional precomputed mean (R's var shape: pass nil to compute it);
  the opts table is TRAILING. num: centered pass + fused vdot - the
  same sequential accumulation as the pure loop (bitwise). Throws on
  fewer than 2 observations or a non-number mean."
  [xs &opt m opts]
  (default opts {})
  # the m-slot trap - an opts table bound
  # into m fails deep in the kernel; guard it actionably
  (when (and m (not (number? m)))
    (error "stats: variance expected a number (precomputed mean) in slot 2 - pass nil to compute it, opts table trailing"))
  (default m (mean xs opts))
  (def n (num/count xs))
  (when (< n 2) (error "stats: variance needs >= 2 observations"))
  (if (= :num (num/resolve-impl opts))
    (let [k (num/kernels)
          c ((k 'vmap1-scalar) :sub (num/stage xs k) m)]
      (/ ((k 'vdot) c c) (- n 1)))
    (do
      (num/check-impl-landed opts xs)
      (var s 0)
      (each x xs (+= s (* (- x m) (- x m))))
      (/ s (- n 1)))))

(defn sd
  "Sample standard deviation (sqrt of variance); slot 2 is an
  optional precomputed mean, opts table trailing."
  [xs &opt m opts]
  (default opts {})
  (when (and m (not (number? m)))
    (error "stats: sd expected a number (precomputed mean) in slot 2 - pass nil to compute it, opts table trailing"))
  (default m (mean xs opts))
  (math/sqrt (variance xs m opts)))

(defn skewness
  "Sample skewness G1, R e1071 type=2: (n/((n-1)(n-2))) *
  sum((x-m)^3)/s^3. Slots 2/3 are optional precomputed mean and sd;
  the opts table is trailing. num: the m3 dot product precomputes
  c-squared elementwise (association differs from the pure (d*d)*d;
  parity bar 1e-12). Throws on fewer than 3 observations."
  [xs &opt m s opts]
  (default opts {})
  # the m/s slot traps, guarded
  (when (or (and m (not (number? m))) (and s (not (number? s))))
    (error "stats: skewness expected numbers (precomputed moments) in slots 2/3 - pass nil to compute them, opts table trailing"))
  (default m (mean xs opts))
  (default s (sd xs m opts))
  (def n (num/count xs))
  (when (< n 3) (error "stats: skewness needs >= 3 observations"))
  (if (= :num (num/resolve-impl opts))
    (let [k (num/kernels)
          c ((k 'vmap1-scalar) :sub (num/stage xs k) m)
          c2 ((k 'vmap2) :mul c c)
          m3 ((k 'vdot) c c2)]
      (* (/ n (* (- n 1) (- n 2))) (/ m3 (* s s s))))
    (do
      (num/check-impl-landed opts xs)
      (var m3 0)
      (var m2 0)
      (each x xs
        (def d (- x m))
        (+= m3 (* d d d))
        (+= m2 (* d d)))
      (* (/ n (* (- n 1) (- n 2))) (/ m3 (* s s s))))))

(defn kurtosis
  "Excess kurtosis G2, R e1071 type=2 (default):
  ((n+1)n/((n-1)(n-2)(n-3)))*m4/s^4 - 3(n-1)^2/((n-2)(n-3)).
  Slots 2/3 are optional precomputed mean and sd; the opts table
  is trailing. Throws on fewer than 4 observations."
  [xs &opt m s opts]
  (default opts {})
  # the m/s slot traps, guarded
  (when (or (and m (not (number? m))) (and s (not (number? s))))
    (error "stats: kurtosis expected numbers (precomputed moments) in slots 2/3 - pass nil to compute them, opts table trailing"))
  (default m (mean xs opts))
  (default s (sd xs m opts))
  (def n (num/count xs))
  (when (< n 4) (error "stats: kurtosis needs >= 4 observations"))
  (if (= :num (num/resolve-impl opts))
    (let [k (num/kernels)
          c ((k 'vmap1-scalar) :sub (num/stage xs k) m)
          c2 ((k 'vmap2) :mul c c)
          m4 ((k 'vdot) c2 c2)]
      (- (* (/ (* (+ n 1) n) (* (- n 1) (- n 2) (- n 3)))
            (/ m4 (* s s s s)))
         (* 3 (/ (* (- n 1) (- n 1)) (* (- n 2) (- n 3))))))
    (do
      (num/check-impl-landed opts xs)
      (var m4 0)
      (each x xs
        (def d (- x m))
        (+= m4 (* d d d d)))
      (- (* (/ (* (+ n 1) n) (* (- n 1) (- n 2) (- n 3)))
            (/ m4 (* s s s s)))
         (* 3 (/ (* (- n 1) (- n 1)) (* (- n 2) (- n 3))))))))

#==------------------------------------------------------------------------==#
# quantiles - R type-7 (the default: h = (n-1)p + 1, linear interpolation)
#==------------------------------------------------------------------------==#

(defn/typed sorted-copy {:args [:any] :ret :array} [xs]
  "Sorted copy of an array."
  (sorted (array/slice xs)))

(defn quantile-t7
  "Type-7 quantile (R default) of a SORTED array: h = (n-1)p,
  linear interpolation between adjacent order statistics. Throws on
  p outside [0,1] or an empty input. Sort-bound composition; opts:
  {:impl}."
  [sorted-xs p &opt opts]
  (default opts {})
  # sorted input; sort-bound composition:
  # :num runs the identical janet interpolation
  (arrays-only opts sorted-xs)
  (def n (length sorted-xs))
  (when (or (< p 0) (> p 1)) (errorf "stats: quantile p in [0,1], got %q" p))
  (when (zero? n) (error "stats: quantile of empty"))
  (def h (* p (- n 1)))
  (def lo (math/floor h))
  (def hi (min (+ lo 1) (- n 1)))
  (+ (get sorted-xs lo) (* (- (get sorted-xs hi) (get sorted-xs lo)) (- h lo))))

(defn median
  "Median (the type-7 quantile at 0.5). Opts: {:impl}."
  [xs &opt opts]
  (default opts {})
  (arrays-only opts xs)
  (quantile-t7 (sorted-copy xs) 0.5 opts))

(defn iqr
  "Interquartile range, Q3 - Q1 by the type-7 convention. Opts:
  {:impl}."
  [xs &opt opts]
  (default opts {})
  (arrays-only opts xs)
  (def s (sorted-copy xs))
  (- (quantile-t7 s 0.75 opts) (quantile-t7 s 0.25 opts)))

(defn/typed sig {:args [:number] :ret :number} [x]
  "sign(x): 1 / -1 / 0."
  (cond (pos? x) 1 (neg? x) -1 0))

#==------------------------------------------------------------------------==#
# correlations
#==------------------------------------------------------------------------==#

(defn/typed check-pairs {:args [:any :any] :ret :number} [xs ys]
  "Correlation-pair guard: equal lengths, >= 2 pairs; returns n."
  # num/count: arrays and packed handles both report their length
  (def n (num/count xs))
  (unless (= n (num/count ys)) (errorf "stats: length mismatch %d vs %d" n (num/count ys)))
  (when (< n 2) (error "stats: correlations need >= 2 pairs"))
  n)

(defn pearson
  "Pearson product-moment correlation, computed from the five raw
  sums in one pass (sx/sy/sxx/syy/sxy; the R-verified closed form).
  Resident input: arrays stage in, packed handles chain through.
  Throws on a length mismatch or fewer than 2 pairs. Opts: {:impl}."
  [xs ys &opt opts]
  (default opts {})
  # the five raw sums - sx/sy via vreduce :sum,
  # sxx/syy/sxy via fused vdot - same products in the same order as
  # the pure single pass (bitwise by construction)
  (def n (check-pairs xs ys))
  (if (= :num (num/resolve-impl opts))
    (let [k (num/kernels)
          xf (num/stage xs k)
          yf (num/stage ys k)
          sx ((k 'vreduce) :sum xf)
          sy ((k 'vreduce) :sum yf)
          sxx ((k 'vdot) xf xf)
          syy ((k 'vdot) yf yf)
          sxy ((k 'vdot) xf yf)]
      (def cov (- sxy (/ (* sx sy) n)))
      (def vx (- sxx (/ (* sx sx) n)))
      (def vy (- syy (/ (* sy sy) n)))
      (/ cov (math/sqrt (* vx vy))))
    (do
      (arrays-only opts xs ys)
      (var sx 0) (var sy 0) (var sxx 0) (var syy 0) (var sxy 0)
      (for i 0 n
        (def x (xs i))
        (def y (ys i))
        (+= sx x) (+= sy y) (+= sxx (* x x)) (+= syy (* y y)) (+= sxy (* x y)))
      (def cov (- sxy (/ (* sx sy) n)))
      (def vx (- sxx (/ (* sx sx) n)))
      (def vy (- syy (/ (* sy sy) n)))
      (/ cov (math/sqrt (* vx vy))))))

(defn/typed rank-ties {:args [:array] :ret :array} [sorted-with-idx]
  "Average ranks for ties (R rank default); input: array of [val
  idx] pairs already sorted by val; output: ranks indexed by
  original position (1-based)."
  (def n (length sorted-with-idx))
  (def ranks (array/new-filled n 0))
  (var i 0)
  (while (< i n)
    (def v (get (get sorted-with-idx i) 0))
    (var j i)
    (while (and (< j n) (= (get (get sorted-with-idx j) 0) v)) (++ j))
    # positions i..j-1 share rank (i + j - 1)/2 + 1 (1-based average)
    (def avg-rank (/ (+ i j 1) 2))
    (for k i j (put ranks (get (get sorted-with-idx k) 1) avg-rank))
    (set i j))
  ranks)

(defn/typed ranks {:args [:any] :ret :array} [xs]
  "Ranks with average ties (R rank default)."
  (def with-idx (seq [i :range [0 (length xs)]] [(xs i) i]))
  (def sorted-pts (sorted with-idx (fn [a b] (< (a 0) (b 0)))))
  (rank-ties sorted-pts))

(defn spearman
  "Spearman rank correlation: pearson of the ranks (average ranks
  for ties, R rank default). Sort-bound composition (the janet rank
  sort has no kernel path); the correlation itself may finish on the
  kernel path under :num. Opts: {:impl}."
  [xs ys &opt opts]
  (default opts {})
  # janet ranks (the sort has no kernel
  # path), then pearson at its own resolution - :num finishes the
  # correlation on the kernel path
  (arrays-only opts xs ys)
  (pearson (ranks xs) (ranks ys) opts))

(defn/typed kendall-core {:args [:any :any] :ret :number} [xs ys]
  "Kendall tau-b (R cor method='kendall'): pairs tied in BOTH
  variables excluded; denominator sqrt((C+D+Tx)(C+D+Ty))."
  (def n (check-pairs xs ys))
  (var c 0)
  (var d 0)
  (var tx 0)
  (var ty 0)
  (for i 0 (- n 1)
    (for j (+ i 1) n
      (def dx (sig (- (xs i) (xs j))))
      (def dy (sig (- (ys i) (ys j))))
      (cond
        (= dx 0) (if (= dy 0) nil (++ tx))
        (= dy 0) (++ ty)
        (pos? (* dx dy)) (++ c)
        (++ d))))
  (def denom (math/sqrt (* (+ c d tx) (+ c d ty))))
  (if (zero? denom) 0 (/ (- c d) denom)))

(defn kendall
  "Kendall tau-b rank correlation (see kendall-core). Sort-bound
  composition (branchy pairwise sign loop - no kernel applies).
  Opts: {:impl}."
  [xs ys &opt opts]
  (default opts {})
  # the pairwise sign loop is branchy
  # scalar work - no vector kernel applies; :num runs it identically
  (arrays-only opts xs ys)
  (kendall-core xs ys))

(defn/typed mad-core {:args [:any :any :any] :ret :number} [xs m c]
  "Median absolute deviation x c around the center m (R mad
  passes the median as m)."
  (def s (sorted-copy (map |(math/abs (- $ m)) xs)))
  (* c (quantile-t7 s 0.5 {})))

(defn mad
  "Median absolute deviation, R mad: median(|x - center|) x
  {:constant 1.4826}. Slot 2 is the CENTER (R's second positional,
  default the median); :constant carries R's third positional
  ((mad @[1 2 3 4] 1 {:constant 2}) = 3.0, R's own value); the opts
  table is trailing. Sort-bound composition. Opts: {:impl}."
  [xs &opt m opts]
  (default opts {})
  (def c (get opts :constant 1.4826))
  # the m-slot trap, guarded
  (when (and m (not (number? m)))
    (error "stats: mad expected a number (the center) in slot 2 - pass nil for the median, opts table trailing"))
  # the center + the absolute-deviation
  # quantile are janet sorts - :num runs them identically
  (arrays-only opts xs)
  (default m (median xs opts))
  (mad-core xs m c))

(defn trimmed-mean
  "Symmetrically trimmed mean, R mean(trim=): a fraction trim of n
  observations (floor) is dropped from EACH end before averaging.
  Sort-bound composition; throws when the trim removes everything.
  Opts: {:impl}."
  [xs trim &opt opts]
  (default opts {})
  # (fraction of n trimmed each side, floor)
  (arrays-only opts xs)
  (def s (sorted-copy xs))
  (def n (length s))
  (def k (math/floor (* trim n)))
  (when (>= (* 2 k) n) (error "stats: trim removes everything"))
  (var sum 0)
  (var count 0)
  (for i k (- n k)
    (+= sum (s i))
    (++ count))
  (/ sum count))

(defn/typed ecdf-count {:args [:any :number] :ret :number} [s q]
  "#{sorted[i] <= q} by bisect-right - the right-continuous step:
  F(x) counts every observation <= x. Internal to ecdf."
  (def n (length s))
  (var lo 0)
  (var hi n)
  (while (< lo hi)
    (def mid (math/floor (/ (+ lo hi) 2)))
    (if (<= (s mid) q) (set lo (+ mid 1)) (set hi mid)))
  lo)

(defn ecdf
  "Empirical CDF as a CALLABLE (R ecdf): F(x) = #{x_i <= x}/n,
  right-continuous, ties share one step value, F(max) = 1. F
  accepts a scalar or an indexed grid (scalar/array out, R's
  vectorized evaluation). Throws on an empty or non-finite input;
  packed handles are rejected (evaluation is a binary search, not
  a kernel op - SPEC ecdf-ks.qnt). Opts: {:impl}."
  [xs &opt opts]
  (default opts {})
  # R ecdf; sort-bound composition - the
  # sort runs identical janet code under both impls.
  (arrays-only opts xs)
  (def n (length xs))
  (when (zero? n) (error "stats: ecdf needs at least one observation"))
  (each x xs
    (unless (and (< x math/inf) (> x (- math/inf)))
      (error "stats: ecdf needs finite observations")))
  (def s (sorted-copy xs))
  (fn ecdf-fn [q]
    (if (indexed? q)
      (map |(/ (ecdf-count s $) n) q)
      (/ (ecdf-count s q) n))))
