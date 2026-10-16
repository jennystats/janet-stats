# stats/utils.janet - numeric-stability plumbing.
# lgamma (Lanczos), erf/erfc via the incomplete gamma, continued
# fractions for P(a,x) / I_x(a,b) - each R-conformant at extreme
# args. Pure janet, zero deps; num-independent (the pure path is the
# reference lane).
(import typed :prefix "")

(def SQRT-2PI (math/sqrt (* 2 math/pi)))
(def LN-SQRT-2PI 0.91893853320467274178032973640562)

(defn/typed lgamma-stirling {:args [:number] :ret :number} [z]
  "Stirling asymptotic (Bernoulli) series for log Gamma; error < 1e-15
  for z >= 10."
  (def z2 (* z z))
  (def z3 (* z2 z))
  (def z5 (* z3 z2))
  (def z7 (* z5 z2))
  (def z9 (* z7 z2))
  (def z11 (* z9 z2))
  (def z13 (* z11 z2))
  (+ (* (- z 0.5) (math/log z))
     (- z)
     LN-SQRT-2PI
     (/ 1 (* 12 z))
     (- (/ 1 (* 360 z3)))
     (/ 1 (* 1260 z5))
     (- (/ 1 (* 1680 z7)))
     (/ 1 (* 1188 z9))
     (- (/ 691 (* 360360 z11)))
     (/ 1 (* 156 z13))))

(defn/typed lgamma-pos {:args [:number] :ret :number} [x]
  "log Gamma for x >= 0.5: shift up to >= 10 via the recurrence,
  Stirling there."
  (var n 0)
  (var z x)
  (while (< z 10) (++ n) (++ z))
  (var acc (lgamma-stirling z))
  (for i 0 n (-= acc (math/log (+ x i))))
  acc)

(defn/typed lgamma {:args [:number] :ret :number} [x]
  "log |Gamma(x)|; reflection (|Gamma(x)Gamma(1-x)| = pi/|sin(pi x)|)
  for x < 0.5. Errors at the poles (non-positive integers)."
  (when (and (= x (math/floor x)) (<= x 0))
    (errorf "stats: lgamma pole at %q" x))
  (if (< x 0.5)
    (do # reflection: |Γ(x)Γ(1-x)| = π/|sin(πx)|
      (- (math/log (/ math/pi (math/abs (math/sin (* math/pi x)))))
         (lgamma-pos (- 1 x))))
    (lgamma-pos x)))

(defn/typed log-beta {:args [:number :number] :ret :number} [a b]
  "log Beta(a,b) = lgamma(a) + lgamma(b) - lgamma(a+b)."
  (- (+ (lgamma a) (lgamma b)) (lgamma (+ a b))))

(def GAMMA-MAXIT 1000)
(def GAMMA-EPS 1e-15)
(def FPMIN 1e-300)

(defn/typed gamma-p-series {:args [:number :number] :ret :number} [a x]
  "Series evaluation of the regularized lower incomplete gamma
  P(a,x); the convergent lane for x < a+1. Internal to gamma-p."
  (var ap a)
  (var sum (/ 1 a))
  (var del sum)
  (for i 0 GAMMA-MAXIT
    (+= ap 1)
    (set del (* del (/ x ap)))
    (+= sum del)
    (when (< (math/abs del) (* (math/abs sum) GAMMA-EPS)) (break)))
  (* sum (math/exp (- (+ (* a (math/log x)) (- x)) (lgamma a)))))

(defn/typed gamma-q-cf {:args [:number :number] :ret :number} [a x]
  "Continued-fraction evaluation of the regularized upper
  incomplete gamma Q(a,x); the convergent lane for x >= a+1
  (modified Lentz). Internal to gamma-q."
  (def tiny 1e-300)
  (var b (+ x 1 (- a)))
  (var c (/ 1 tiny))
  (var d (/ 1 b))
  (var h d)
  (var done false)
  (for i 1 GAMMA-MAXIT
    (when (not done)
      (def an (- (* i (- i a))))
      (+= b 2)
      (set d (+ (* an d) b))
      (when (< (math/abs d) tiny) (set d tiny))
      (set c (+ b (/ an c)))
      (when (< (math/abs c) tiny) (set c tiny))
      (set d (/ 1 d))
      (def del (* d c))
      (*= h del)
      (when (< (math/abs (- del 1)) GAMMA-EPS) (set done true))))
  (* h (math/exp (- (+ (* a (math/log x)) (- x)) (lgamma a)))))

(defn/typed gamma-p {:args [:number :number] :ret :number} [a x]
  "Regularized lower incomplete gamma P(a,x). Requires a > 0,
  x >= 0."
  (when (<= a 0) (errorf "stats: gamma-p needs a > 0, got %q" a))
  (when (< x 0) (errorf "stats: gamma-p needs x >= 0, got %q" x))
  (cond
    (= x 0) 0.0
    (< x (+ a 1)) (gamma-p-series a x)
    (- 1 (gamma-q-cf a x))))

(defn/typed gamma-q {:args [:number :number] :ret :number} [a x]
  "Regularized upper incomplete gamma Q(a,x) = 1 - P(a,x).
  Requires a > 0, x >= 0."
  (when (<= a 0) (errorf "stats: gamma-q needs a > 0, got %q" a))
  (when (< x 0) (errorf "stats: gamma-q needs x >= 0, got %q" x))
  (cond
    (= x 0) 1.0
    (< x (+ a 1)) (- 1 (gamma-p-series a x))
    (gamma-q-cf a x)))

(defn/typed beta-cf {:args [:number :number :number] :ret :number} [a b x]
  "Continued fraction for the regularized incomplete beta
  (betacf). Internal to beta-inc."
  (def tiny 1e-300)
  (def qab (+ a b))
  (def qap (+ a 1))
  (def qam (- a 1))
  (var c 1.0)
  (var d (- 1 (/ (* qab x) qap)))
  (when (< (math/abs d) tiny) (set d tiny))
  (set d (/ 1 d))
  (var h d)
  (var done false)
  (var m 1)
  (while (and (<= m 5000) (not done))
    (def m2 (* 2 m))
    (def aa (/ (* m (- b m) x) (* (+ qam m2) (+ a m2))))
    (set d (+ 1 (* aa d)))
    (when (< (math/abs d) tiny) (set d tiny))
    (set c (+ 1 (/ aa c)))
    (when (< (math/abs c) tiny) (set c tiny))
    (set d (/ 1 d))
    (*= h (* d c))
    (def aa2 (/ (* (- (+ a m)) (+ qab m) x) (* (+ a m2) (+ qap m2))))
    (set d (+ 1 (* aa2 d)))
    (when (< (math/abs d) tiny) (set d tiny))
    (set c (+ 1 (/ aa2 c)))
    (when (< (math/abs c) tiny) (set c tiny))
    (set d (/ 1 d))
    (def del (* d c))
    (*= h del)
    (when (< (math/abs (- del 1)) (* 3 1e-15)) (set done true))
    (++ m))
  h)

(defn/typed packed-handle? {:args [:any] :ret :boolean} [x]
  "True for a full packed-handle marker: :n and :bytes numbers plus
  a :ptr buffer of at least :bytes bytes - buffer sufficiency is
  part of shape, because the kernel validates :bytes == :n x 8
  (element width) and then reads :n elements, so a short :ptr is an
  out-of-bounds read. Partial markers and tables are caller misuse
  (see api-contract input dispatch). Shape is gated, not content:
  caller-supplied :ptr bytes unpack as data. Real pack-f64 handles
  are exact, so there are no false positives; :n 0 needs zero
  bytes."
  (and (struct? x)
       (number? (get x :n))
       (number? (get x :bytes))
       (buffer? (get x :ptr))
       (>= (length (get x :ptr)) (get x :bytes))))

(defn/typed beta-inc {:args [:number :number :number] :ret :number} [a b x]
  "Regularized incomplete beta I_x(a,b) for 0 <= x <= 1; the
  symmetric swap keeps the continued fraction converging. Requires
  a, b > 0."
  (when (or (<= a 0) (<= b 0)) (errorf "stats: beta-inc needs a,b > 0"))
  (when (or (< x 0) (> x 1)) (errorf "stats: beta-inc needs 0<=x<=1, got %q" x))
  (cond
    (= x 0) 0.0
    (= x 1) 1.0
    (> x (/ (+ a 1) (+ a b 2))) (- 1 (beta-inc b a (- 1 x)))
    (do
      # log-space combination: bt and h can each over/underflow at
      # extreme params; only their product is well-scaled. h > 0 for
      # all valid (a,b,x).
      (def lbt (+ (- (log-beta a b))
                 (* a (math/log x))
                 (* b (math/log1p (- x)))))
      (* (math/exp (+ lbt (math/log (beta-cf a b x)))) (/ 1 a)))))

(defn/typed erf {:args [:number] :ret :number} [x]
  "Error function via P(1/2, x^2) - the same machinery as the
  chi-squared and Poisson CDFs. The +inf pole returns the limit 1
  (the continued fraction alone returns nan; the libm kernel and R
  give the limit). Odd by reflection."
  (if (neg? x)
    (- (erf (- x)))
    (if (= x math/inf) 1.0 (gamma-p 0.5 (* x x)))))

(defn/typed erfc {:args [:number] :ret :number} [x]
  "Complementary error function via Q(1/2, x^2) - no cancellation
  in the deep tail. The +inf pole returns the limit 0; erfc(-inf) = 2
  by reflection."
  (if (neg? x) (- 2 (erfc (- x)))
    (if (= x math/inf) 0.0 (gamma-q 0.5 (* x x)))))

(defn/typed logsumexp {:args [:array] :ret :number} [xs]
  "log(sum(exp(xs))) computed stably; an empty array gives -inf."
  (if (empty? xs) math/-inf
    (do
      (def m (reduce (fn [acc x] (if (> x acc) x acc)) (- math/inf) xs))
      (if (= m (- math/inf)) (- math/inf)
        (+ m (math/log (reduce (fn [acc x] (+ acc (math/exp (- x m)))) 0 xs)))))))

(defn/typed solve-monotone {:args [:function :number :number :number] :ret :number} [f lo hi target]
  "Solve f(x) = target for f monotone increasing on [lo, hi] -
  hybrid bisection/secant. The secant step is only accepted inside
  the middle half of the bracket: a stale far anchor against a
  saturated tail (f == 1.0 for all large x) otherwise makes the
  secant crawl by ~5% per iteration instead of bisecting. Shared by
  the q-family inversions."
  (var a lo)
  (var b hi)
  (var fa (- (f a) target))
  (when (> fa 0) (errorf "stats: solve-monotone not bracketed at lo (f=%q)" (+ fa target)))
  (var fb (- (f b) target))
  (var x (* 0.5 (+ a b)))
  (for _ 0 200
    (def fx (- (f x) target))
    (if (> fx 0) (do (set b x) (set fb fx)) (do (set a x) (set fa fx)))
    (when (< (- b a) (* 1e-15 (+ 1 (math/abs x)))) (break))
    (def x-next (- b (/ (* fb (- b a)) (- fb fa))))
    (def q1 (+ a (* 0.25 (- b a))))
    (def q2 (- b (* 0.25 (- b a))))
    (set x (if (and (>= x-next q1) (<= x-next q2)) x-next (* 0.5 (+ a b)))))
  x)
