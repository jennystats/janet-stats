# linalg.janet - small-dim linear algebra with R conventions
# (pure-janet reference lane; janet-num kernels behind {:impl
# :num}). Row-major arrays-of-arrays throughout.
# The landed kernels resolve through num/resolve-impl ({:impl
# :impl} routes dynamically); qr is documented composition (branchy
# Householder loops - no kernel applies)
(import typed :prefix "")

(import ./num)

(defn/typed m-rows {:args [:any] :ret :number} [m] "Row-major matrix row count." (length m))
(defn/typed m-cols {:args [:any] :ret :number} [m] "Row-major matrix column count." (length (get m 0)))

(defn/typed check-square {:args [:any] :ret :nil} [a]
  "Assert a square matrix. Internal."
  (def n (length a))
  (when (or (< n 1) (not= n (length (get a 0))))
    (error "stats: linalg needs a square matrix")))

(defn/typed zeros-m {:args [:number :number] :ret :array} [rows cols]
  "Zero matrix with distinct row arrays - array/new-filled shares
  the fill value between rows."
  (def out @[])
  (for i 0 rows (array/push out (array/new-filled cols 0)))
  out)

(defn/typed copy-m {:args [:any] :ret :array} [m]
  "Deep-copy a row-major matrix."
  (map |(array ;$) m))

(defn/typed k* {:args [] :ret :any} []
  "LAZY kernel accessor - the require runs ONLY inside a :num call
  site (the pure path never loads janet-num). defn/typed bindings
  are ref-tables {:value ...} - unwrap on access."
  (def m (require "janet-num"))
  (fn [nm]
    (def w (get m nm))
    (if (table? w) (get w :value) w)))

(defn/typed flat-row-major {:args [:any] :ret :any} [a]
  "Flatten a row-major matrix into a flat array for pack-f64."
  (def out @[])
  (each row a (each v row (array/push out v)))
  out)

(defn/typed rebuild-rows {:args [:any :number :number] :ret :array} [flat n p]
  "Rebuild n x p row-major rows from a flat array."
  (seq [i :range [0 n]] (array ;(slice flat (* i p) (* (inc i) p)))))

#==------------------------------------------------------------------------==#
# mmul - row-major C = A(m x k) . B(k x n)
#==------------------------------------------------------------------------==#

(defn mmul
  "Row-major matrix product C = A(m x k) . B(k x n). Throws on
  dimension mismatch. Opts: {:impl}."
  [a b &opt opts]
  (default opts {})
  (num/check-impl-landed opts)
  (def m (m-rows a))
  (def k (m-cols a))
  (def n (m-cols b))
  (when (not= k (m-rows b)) (error "stats: mmul dimension mismatch"))
  (if (= :num (num/resolve-impl opts))
    (let [kn (k*)
          fm ((kn 'pack-f64) (flat-row-major a))
          fb ((kn 'pack-f64) (flat-row-major b))
          out ((kn 'unpack-f64) ((kn 'matmul) fm fb m k n))]
      (rebuild-rows out m n))
    (let [out (zeros-m m n)]
      (for i 0 m
        (for j 0 n
          (var s 0)
          (for t 0 k (+= s (* (get (get a i) t) (get (get b t) j))))
          (put (get out i) j s)))
      out)))

#==------------------------------------------------------------------------==#
# LU solve with partial pivoting (R solve(a, b) semantics).
#==------------------------------------------------------------------------==#

(defn lu-solve-pure
  "Pure-lane LU solve with partial pivoting (the reference for
  solve's num lane). Internal."
  [a b]
  (def n (length a))
  (check-square a)
  (when (not= (length b) n) (error "stats: solve needs b of length n"))
  # augmented LU with partial pivoting (Gauss-Jordan-free: standard
  # Doolittle + row swaps; mirrors LAPACK dgesv structure)
  (def lu (map |(array ;$) a))
  (def piv (seq [i :range [0 n]] i))
  (def x (array ;b))
  (for col 0 n
    # pivot: largest |a[i][col]| in the remaining submatrix
    (var piv-row col)
    (var best (math/abs (get (get lu col) col)))
    (for r (inc col) n
      (def v (math/abs (get (get lu r) col)))
      (when (> v best) (set best v) (set piv-row r)))
    (when (< best 1e-300)
      (error "stats: solve: singular (or near-singular) matrix"))
    (unless (= piv-row col)
      (def tmp (get lu col))
      (put lu col (get lu piv-row))
      (put lu piv-row tmp)
      (def tx (get x col))
      (put x col (get x piv-row))
      (put x piv-row tx))
    # eliminate below
    (def diag (get (get lu col) col))
    (for r (inc col) n
      (def factor (/ (get (get lu r) col) diag))
      (when (not= factor 0)
        (for c col n
          (put (get lu r) c (- (get (get lu r) c) (* factor (get (get lu col) c)))))
      (put x r (- (get x r) (* factor (get x col)))))))
  # back substitution - janet `for` CANNOT iterate downward (zero
  # iterations even with a negative step) and loop :down is EXCLUSIVE
  # of its stop (unlike :range - so the stop must be -1 to reach 0)
  (loop [i :down [(dec n) -1]]
    (var s (get x i))
    (for j (inc i) n (-= s (* (get (get lu i) j) (get x j))))
    (put x i (/ s (get (get lu i) i))))
  x)

(defn solve
  "Solve a x = b (R solve(a, b) semantics) by LU with partial
  pivoting; b is a vector of length n. Throws on singular or
  near-singular matrices. Opts: {:impl}."
  [a b &opt opts]
  (default opts {})
  (num/check-impl-landed opts)
  (if (= :num (num/resolve-impl opts))
    (let [kn (k*)
          fa ((kn 'pack-f64) (flat-row-major a))
          fb ((kn 'pack-f64) b)
          n (length a)
          status ((kn 'lu-solve) fa fb n)]
      (unless status (error "stats: solve: singular (or near-singular) matrix"))
      ((kn 'unpack-f64) fb))
    (lu-solve-pure a b)))

(defn inverse
  "Matrix inverse via columnwise solves against the identity.
  Opts: {:impl}; throws on singular input."
  [a &opt opts]
  (default opts {})
  (num/check-impl-landed opts)
  (def n (length a))
  (check-square a)
  # solve A x = e_i columnwise
  (def cols (seq [i :range [0 n]]
              (solve a (seq [j :range [0 n]] (if (= i j) 1.0 0.0)) opts)))
  (rebuild-rows (mapcat identity cols) n n))

#==------------------------------------------------------------------------==#
# Cholesky - UPPER factor U with A = U'U (R chol convention).
# The kernel produces LOWER L (A = L L'); pure mirrors it and
# transposes to U so both impls agree on the returned factor.
#==------------------------------------------------------------------------==#

(defn chol-lower-pure
  "Pure-lane lower Cholesky factor (the reference for chol's num
  lane). Internal."
  [a]
  (def n (length a))
  (check-square a)
  (def l (zeros-m n n))
  (for i 0 n
    (for j 0 (inc i)
      (var s (get (get a i) j))
      (for k 0 j (-= s (* (get (get l i) k) (get (get l j) k))))
      (cond
        (= i j)
        (do
          (when (<= s 0)
            (error "stats: chol: matrix not positive-definite"))
          (put (get l i) j (math/sqrt s)))
        (do
          (put (get l i) j (/ s (get (get l j) j)))))))
  l)

(defn chol
  "Cholesky factor, returning the UPPER factor U with A = U'U (R
  chol convention); both impls agree on the returned factor.
  Throws on non-positive-definite input. Opts: {:impl}."
  [a &opt opts]
  (default opts {})
  (num/check-impl-landed opts)
  (def l
    (if (= :num (num/resolve-impl opts))
      (let [kn (k*)
            fa ((kn 'pack-f64) (flat-row-major a))
            n (length a)
            [status lh] ((kn 'cholesky-l) fa n)]
        (unless status (error "stats: chol: matrix not positive-definite"))
        (rebuild-rows ((kn 'unpack-f64) lh) n n))
      (chol-lower-pure a)))
  # return the UPPER factor (R chol: A = U'U); U = L' - the UPPER
  # region (j >= i) of U carries L[j][i], the lower region is zero
  (def n (length l))
  (seq [i :range [0 n]]
    (seq [j :range [0 n]]
      (if (>= j i) (get (get l j) i) 0.0))))

(defn chol-solve
  "Solve U'U x = b given the UPPER factor u from chol; b is a
  vector of length n. Opts: {:impl}; throws on a failed solve."
  [u b &opt opts]
  # u: the UPPER factor from chol; solves U'U x = b
  (default opts {})
  (num/check-impl-landed opts)
  (def n (length u))
  (when (not= (length b) n) (error "stats: chol-solve needs b of length n"))
  (if (= :num (num/resolve-impl opts))
    (let [kn (k*)
          # kernel wants the FULL n x n LOWER L; U' = L (pack the transpose)
          lflat @[]
          _ (do (for i 0 n (for j 0 n
                             (array/push lflat (if (<= j i) (get (get u j) i) 0.0)))))
          fl ((kn 'pack-f64) lflat)
          fb ((kn 'pack-f64) b)
          status ((kn 'chol-solve) fl fb n)]
      (unless status (error "stats: chol-solve failed"))
      ((kn 'unpack-f64) fb))
    (do
      # pure: forward-solve U' y = b (U' is lower), then U x = y
      (def y (array ;b))
      (for i 0 n
        (var s (get y i))
        (for j 0 i (-= s (* (get (get u j) i) (get y j))))
        (put y i (/ s (get (get u i) i))))
      (def x (array ;y))
      (loop [i :down [(dec n) -1]]
        (var s (get x i))
        (for j (inc i) n (-= s (* (get (get u i) j) (get x j))))
        (put x i (/ s (get (get u i) i))))
      x)))

#==------------------------------------------------------------------------==#
# QR - Householder, NO pivoting. {:Q :R}, A = Q R, R upper with
# raw Householder diag signs (the signs follow the algorithm);
# conformance via reconstruction + solve.
#==------------------------------------------------------------------------==#

(defn qr
  "QR decomposition by Householder reflections, no pivoting.
  Returns {:Q :R} with A = Q R; R carries the raw Householder
  diagonal signs. Throws when n < p. Composition (branchy scalar
  loops - no kernel applies). Opts: {:impl}."
  [a &opt opts]
  (default opts {})
  # the Householder loop is branchy
  # scalar work - no vector kernel applies; :num runs it identically
  (num/check-impl-landed opts)
  (def n (m-rows a))
  (def p (m-cols a))
  (when (< n p) (error "stats: qr needs n >= p"))
  (def r (map |(array ;$) a))
  # accumulate Householder reflectors, apply to a growing Q
  (def q (zeros-m n n))
  (for i 0 n (put (get q i) i 1.0))
  (for k 0 p
    # x = R[k:n, k]
    (def x (seq [i :range [k n]] (get (get r i) k)))
    (def xnorm (math/sqrt (reduce (fn [acc v] (+ acc (* v v))) 0 x)))
    (when (> xnorm 1e-300)
      # sign choice avoids cancellation
      (def a0 (get x 0))
      (def alpha (- (if (>= a0 0) xnorm (- xnorm))))
      (put x 0 (- (get x 0) alpha))
      (def vnorm2 (reduce (fn [acc v] (+ acc (* v v))) 0 x))
      (when (> vnorm2 1e-300)
        # R <- R - 2 v v' R / vnorm2  (rows/cols k..)
        # apply to R columns k..p-1
        (for c k p
          (var dot 0)
          (for i k n (+= dot (* (get x (- i k)) (get (get r i) c))))
          (def s2 (/ (* 2 dot) vnorm2))
          (for i k n
            (put (get r i) c (- (get (get r i) c) (* s2 (get x (- i k)))))))
        # Q <- Q - 2 (Q v) v' / vnorm2  (columns of Q, rows k..n)
        (for c 0 n
          (var dot 0)
          (for i k n (+= dot (* (get (get q c) i) (get x (- i k)))))
          (def s2 (/ (* 2 dot) vnorm2))
          (for i k n
            (put (get q c) i (- (get (get q c) i) (* s2 (get x (- i k))))))))))
  # q accumulated H1 H2 ... Hp by RIGHT-multiplication (row-wise
  # updates) = Q itself (each reflector symmetric) - NO transpose
  {:Q q
   :R r})