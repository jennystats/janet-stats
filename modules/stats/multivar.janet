# multivar.janet - PCA (jacobi eigensolver) + k-means (Lloyd),
# R-conformant with base R prcomp + kmeans (algorithm="Lloyd").
# pca's covariance + scores go through the GEMM/matvec kernels
# (bitwise: same products in the same order); the eigensolver stays
# janet (small-p, branchy - composition); k-means assignment = the
# dist2-rows kernel (the O(nkp) hot loop) + the masked centroid
# accumulation (vmap2+vreduce per cluster off the transposed columns).
(import typed :prefix "")
(import ./utils :prefix "")
(import ./num)
(import ./linalg)
(import ./regsurv)

(defn/typed m-rows {:args [:any] :ret :number} [m] "Row-major matrix row count." (length m))
(defn/typed m-cols {:args [:any] :ret :number} [m] "Row-major matrix column count." (length (get m 0)))

(defn/typed zeros-m {:args [:number :number] :ret :array} [rows cols]
  "Zero matrix with distinct row arrays - array/new-filled shares
  the fill value between rows."
  (def out @[])
  (for i 0 rows (array/push out (array/new-filled cols 0)))
  out)

(defn/typed copy-m {:args [:any] :ret :array} [m]
  "Deep-copy a row-major matrix."
  (map |(array ;$) m))

(defn/typed mean-of {:args [:any] :ret :number} [xs]
  "Plain mean (the 1.0-force keeps int data honest). Internal."
  (var s 0)
  (each v xs (+= s v))
  (/ s (* 1.0 (length xs))))

#==------------------------------------------------------------------------==#
# jacobi-eig - cyclic jacobi rotations on a SYMMETRIC matrix.
# Off-diagonal mass -> 0; eigenvalues on the diagonal, eigenvectors
# accumulated in the columns of :vectors. Sorted DESCENDING (prcomp
# PC order). Rotation formulas: the classic tau/t/c/s (Golub & Van
# Loan, Symmetric Schur decomposition).
#==------------------------------------------------------------------------==#

(defn/typed off-diag {:args [:any] :ret :number} [a]
  "Frobenius norm of the off-diagonal (the jacobi sweep measure)."
  (var s 0)
  (def n (length a))
  (for i 0 n (for j (+ i 1) n (def v (get (get a i) j)) (+= s (* v v))))
  (math/sqrt s))

(defn jacobi-eig
  "Eigen-decomposition of a SYMMETRIC matrix by cyclic jacobi
  rotations (Golub & Van Loan symmetric Schur). Returns
  {:values :vectors :sweeps} with eigenvalues sorted DESCENDING
  (prcomp PC order) and eigenvectors as the columns of :vectors.
  Opts: {:tol 1e-14 :max-sweeps 100}; throws on non-square,
  non-symmetric, p < 2, or non-convergence. Pure composition (no
  vector kernel - branchy small-p scalar work)."
  [a &opt opts]
  (default opts {})
  # cyclic rotations are branchy
  # small-p scalar work - no vector kernel applies; :num runs it
  # identically
  (num/check-impl-landed opts)
  (def tol (get opts :tol 1e-14))
  (def max-sweeps (get opts :max-sweeps 100))
  (def n (length a))
  (when (or (< n 1) (not= n (length (get a 0))))
    (error "stats: jacobi-eig needs a square symmetric matrix"))
  (when (= n 1) (error "stats: jacobi-eig needs p >= 2"))
  # symmetry check
  (for i 0 n
    (for j (+ i 1) n
      (when (> (math/abs (- (get (get a i) j) (get (get a j) i))) 1e-9)
        (error "stats: jacobi-eig needs a symmetric matrix"))))
  (def work (copy-m a))
  (def v (let [out (zeros-m n n)] (for i 0 n (put (get out i) i 1)) out))
  (var sweeps 0)
  (var converged false)
  (while (and (not converged) (< sweeps max-sweeps))
    (++ sweeps)
    (def o0 (off-diag work))
    (when (<= o0 tol) (set converged true) (break))
    (for i 0 (dec n)
      (for j (inc i) n
        (def aij (get (get work i) j))
        (when (> (math/abs aij) 1e-300)
          (def tau (/ (- (get (get work j) j) (get (get work i) i)) (* 2 aij)))
          # t = sign(tau) / (|tau| + sqrt(1 + tau^2))
          (def t
            (let [ta (math/abs tau)
                  rt (math/sqrt (+ 1 (* tau tau)))]
              (if (>= tau 0)
                (/ 1 (+ ta rt))
                (/ -1 (+ ta rt)))))
          (def c (/ 1 (math/sqrt (+ 1 (* t t)))))
          (def s (* t c))
          # A <- J^T A J, applied only to rows/cols i and j
          (for k 0 n
            (def aik (get (get work i) k))
            (def ajk (get (get work j) k))
            (put (get work i) k (- (* c aik) (* s ajk)))
            (put (get work j) k (+ (* s aik) (* c ajk))))
          (for k 0 n
            (def aki (get (get work k) i))
            (def akj (get (get work k) j))
            (put (get work k) i (- (* c aki) (* s akj)))
            (put (get work k) j (+ (* s aki) (* c akj))))
          # V <- V J (accumulate eigenvectors as columns)
          (for k 0 n
            (def vki (get (get v k) i))
            (def vkj (get (get v k) j))
            (put (get v k) i (- (* c vki) (* s vkj)))
            (put (get v k) j (+ (* s vki) (* c vkj))))))))
  (unless converged
    (errorf "stats: jacobi-eig did not converge in %d sweeps" max-sweeps))
  # extract + sort descending
  (def order (sorted (seq [i :range [0 n]] i)
                     (fn [x y] (> (get (get work x) x) (get (get work y) y)))))
  (def values (map |(get (get work $) $) order))
  (def vectors
    (let [out (zeros-m n n)]
      (for c 0 n
        (def src (get order c))
        (for r 0 n (put (get out r) c (get (get v r) src))))
      out))
  {:values values :vectors vectors :sweeps sweeps})

#==------------------------------------------------------------------------==#
# pca - prcomp(center=TRUE, scale.=FALSE) convention: sdev = sqrt of
# the covariance eigenvalues (n-1 denominator), PCs by DECREASING
# sdev, column signs arbitrary (fixtures compare |loadings| and
# scores up to sign).
#==------------------------------------------------------------------------==#

(defn pca
  "Principal component analysis, R prcomp(center=TRUE, scale.=FALSE)
  convention. x = n x p row-major matrix (array of rows); returns
  {:sdev :rotation :scores :var-explained :center :scale :values} -
  sdev from the covariance eigenvalues (n-1 denominator), PCs by
  DECREASING sdev, rotation column signs arbitrary. Opts:
  {:center true :scale false :impl}. Throws when n <= p."
  [x &opt opts]
  (default opts {})
  (num/check-impl-landed opts)
  (def center? (get opts :center true))
  (def scale? (get opts :scale false))
  (def n (m-rows x))
  (def p (m-cols x))
  (when (< n (+ p 1)) (error "stats: pca needs n > p"))
  # column means / sds
  (def means (seq [j :range [0 p]] (mean-of (map |(get $ j) x))))
  (def sds
    (when scale?
      (seq [j :range [0 p]]
        (var s 0)
        (each row x (def d (- (get row j) (get means j))) (+= s (* d d)))
        (math/sqrt (/ s (* 1.0 (- n 1)))))))
  (def xc
    (map (fn [row]
           (seq [j :range [0 p]]
             (def v (- (get row j) (get means j)))
             (if sds (/ v (get sds j)) v)))
         x))
  # covariance matrix, n-1 denominator - cov = X'X/(n-1) via
  # the GEMM (division per element, matching the pure rounding);
  # the eigensolver + tail stay janet
  (def cov
    (if (= :num (num/resolve-impl opts))
      (let [kk (num/kernels)]
        (map (fn [row] (map |(/ $ (* 1.0 (- n 1))) row))
             (linalg/rebuild-rows
              ((kk 'unpack-f64)
               ((kk 'matmul) (num/stage (regsurv/xt-flat xc n p) kk)
                             (num/stage (linalg/flat-row-major xc) kk) p n p))
              p p)))
      (do
        (def covp (zeros-m p p))
        (for i 0 p
         (for j 0 p
            (var s 0)
            (each row xc (+= s (* (get row i) (get row j))))
            (put (get covp i) j (/ s (* 1.0 (- n 1))))))
        covp)))
  (def eig (jacobi-eig cov {}))
  (def values (eig :values))
  (def rotation (eig :vectors))
  # scores = xc . rotation - num: one GEMM
  (def scores
    (if (= :num (num/resolve-impl opts))
      (let [kk (num/kernels)]
        (linalg/rebuild-rows
         ((kk 'unpack-f64)
         ((kk 'matmul) (num/stage (linalg/flat-row-major xc) kk)
                        (num/stage (linalg/flat-row-major rotation) kk) n p p))
         n p))
      (map (fn [row]
             (seq [j :range [0 p]]
               (var s 0)
               (for k2 0 p (+= s (* (get row k2) (get (get rotation k2) j))))
               s))
           xc)))
  (def total (reduce + 0 values))
  {:sdev (map (fn [v] (math/sqrt v)) values)
   :rotation rotation
   :scores scores
   :var-explained (map (fn [v] (if (zero? total) 0 (/ v total))) values)
   :center means :scale sds :values values})

#==------------------------------------------------------------------------==#
# k-means - Lloyd's algorithm (R kmeans algorithm="Lloyd"): assignment
# by nearest centroid (euclidean), centroids = member means;
# converged when no point moves. Seeded Forgy init (k distinct rows,
# math/rng deterministic per seed) or explicit :centers (the R
# same-start conformance path - RNG paths are never cross-compared).
#==------------------------------------------------------------------------==#

(defn/typed dist2 {:args [:any :any] :ret :number} [a b]
  "Squared euclidean distance."
  (var s 0)
  (for j 0 (length a) (def d (- (get a j) (get b j))) (+= s (* d d)))
  s)

(defn k-means
  "Lloyd k-means (R kmeans, algorithm='Lloyd'). x = n x p row-major
  matrix; returns {:centers :labels :withinss :tot-withinss :iters
  :converged} - withinss measured against the FINAL centroids. Opts:
  {:iters 10 :seed \"stats\" :centers nil} - :centers gives k explicit
  starting rows (the R same-start path); otherwise seeded Forgy init
  (deterministic per :seed). Empty clusters keep their centroids.
  Throws unless 1 <= k <= n (integer k)."
  [x k &opt opts]
  (default opts {})
  (num/check-impl-landed opts)
  (def n (m-rows x))
  (def p (m-cols x))
  (unless (and (int? k) (>= k 1) (<= k n))
    (error "stats: k-means needs 1 <= k <= n"))
  (def max-iters (get opts :iters 10))
  (def explicit (get opts :centers))
  (when (and explicit (not= (length explicit) k))
    (error "stats: k-means :centers must have k rows"))
  # initial centroids
  (def centers
    (cond
      explicit (map |(array ;$) explicit)
      (do
        # seeded Forgy: k distinct row indices via math/rng
        (def rng (math/rng (string (get opts :seed "stats"))))
        (def chosen @{})
        (while (< (length chosen) k)
          (def i (math/rng-int rng n))
          (unless (get chosen i) (put chosen i (get x i))))
        (map |(array ;$) (values chosen)))))
  (def labels (array/new-filled n -1))
  (var iters 0)
  (var converged false)
  (var withinss (array/new-filled k 0))
  # x (and its transposed columns) stage ONCE - x never changes;
  # centers repack per iteration (k x p - trivial). Assignment =
  # dist2-rows (same per-element products/order as the pure dist2 -
  # identical labels, identical trajectory)
  (def num? (= :num (num/resolve-impl opts)))
  (def kacc (when num? (num/kernels)))
  (def xf (when num? (num/stage (linalg/flat-row-major x) kacc)))
  (def colf (when num?
              (seq [c :range [0 p]]
                (num/stage (seq [i :range [0 n]] (get (get x i) c)) kacc))))
  (while (and (not converged) (< iters max-iters))
    (++ iters)
    (var moved false)
    # assignment
    (if num?
      (let [D ((kacc 'unpack-f64)
               ((kacc 'dist2-rows) xf
                                  (num/stage (linalg/flat-row-major centers) kacc)
                                  n p k))]
        (for i 0 n
          (def base (* i k))
          (var best 0)
          (var bd (get D base))
          (for c 1 k
            (def d (get D (+ base c)))
            (when (< d bd) (set best c) (set bd d)))
          (unless (= (get labels i) best)
            (set moved true)
            (put labels i best))))
      (for i 0 n
        (var best 0)
        (var bd (dist2 (get x i) (get centers 0)))
        (for c 1 k
          (def d (dist2 (get x i) (get centers c)))
          (when (< d bd) (set best c) (set bd d)))
        (unless (= (get labels i) best)
          (set moved true)
          (put labels i best))))
    # centroid update (member means); empty cluster keeps its centroid.
    # num: the masked accumulation - per
    # cluster, an indicator mask (O(n) glue) then vmap2+vreduce per
    # column off the staged transpose; adding 0s is exact - bitwise
    # with the scatter
    (if num?
      (for c 0 k
        (def mask (num/stage (seq [i :range [0 n]] (if (= (get labels i) c) 1.0 0.0)) kacc))
        (def cnt (math/round ((kacc 'vreduce) :sum mask)))
        (when (pos? cnt)
          (for j 0 p
            (put (get centers c) j
                 (/ ((kacc 'vreduce) :sum ((kacc 'vmap2) :mul (get colf j) mask))
                    (* 1.0 cnt))))))
      (do
        (def sums (zeros-m k p))
        (def counts (array/new-filled k 0))
        (for i 0 n
          (def c (get labels i))
          (put counts c (inc (get counts c)))
          (for j 0 p (put (get sums c) j (+ (get (get sums c) j) (get (get x i) j)))))
        (for c 0 k
          (when (pos? (get counts c))
            (for j 0 p
              (put (get centers c) j (/ (get (get sums c) j) (* 1.0 (get counts c)))))))))
    (unless moved (set converged true)))
  # within-cluster SS against the FINAL centroids - num: one final
  # dist2-rows, masked janet accumulation (same order as pure)
  (set withinss (array/new-filled k 0))
  (if num?
    (let [D ((kacc 'unpack-f64)
             ((kacc 'dist2-rows) xf
                                (num/stage (linalg/flat-row-major centers) kacc)
                                n p k))]
      (for i 0 n
        (def c (get labels i))
        (put withinss c (+ (get withinss c) (get D (+ (* i k) c))))))
    (for i 0 n
      (def c (get labels i))
      (put withinss c (+ (get withinss c) (dist2 (get x i) (get centers c))))))
  (var tot 0)
  (each w withinss (+= tot w))
  {:centers centers :labels labels :withinss withinss :tot-withinss tot
   :iters iters :converged converged})
