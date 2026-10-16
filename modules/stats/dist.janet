# stats/dist.janet - the d/p/q family: normal, t, chi², F,
# binomial, poisson. R-convention semantics (lower.tail, log.p where
# meaningful). Scalar-or-array inputs; {:impl} plumbing. Complete:
# densities (dnorm/dt/dchisq/df/dbinom/
# dpois/dbeta/dgamma) + pnorm/pchisq/ppois/pbeta/pgamma are
# kernel-backed vector paths (vmap ops incl. :exp/:log 0.5.0 + the
# specfun kernels; arrays stage, handles chain); pbinom/pt/pf + the
# q-family are documented composition - the R-exact recurrence,
# per-element sign branches, and monotone inversions run identically
# under :num. VECTOR+num does NO validation pass (the NaN domain
# contract - out-of-domain → NaN, kernels trust) EXCEPT dbinom/
# dpois, whose integer prepass is required (lgamma of a bad k is a
# silent value, not a NaN).
(import typed :prefix "")
(import ./utils :prefix "")
(import ./num)

(defn/typed vapply {:args [:function :any] :ret :any} [f x]
  "Apply f elementwise over an indexed x, or directly to a scalar.
  Internal."
  (if (indexed? x) (map f x) (f x)))

(defn/typed k* {:args [] :ret :any} []
  "The janet-num kernel accessor (lazy require; defn/typed bindings
  are ref-tables - unwrap on access). Internal."
  (num/kernels))

(defn comp-only
  "Documented-composition guard (identical-formula policy): the
  R-exact scalar algorithms run identically under :num; a packed
  handle unpacks first (actionable)."
  [opts x]
  # documented composition (identical-formula policy): R-exact scalar
  # algorithms - monotone inversions, the pbinom recurrence, pt/pf
  # per-element sign branches - run identically under :num (a kernel
  # path would break R-exactness or need reversed-div/branching ops
  # a kernel
  # path would need. Scalars or janet arrays only - a
  # packed handle unpacks first (actionable)
  (num/check-impl-landed opts)
  (unless (or (number? x) (indexed? x))
    (error "stats: composition path takes scalars or janet arrays; unpack packed handles first"))
  nil)

#==------------------------------------------------------------------------==#
# normal
#==------------------------------------------------------------------------==#

(defn dnorm
  "Normal density (log with {:log true}). Kernel-backed on vectors
  (arrays stage, packed handles chain); scalars run the identical
  pure formula. Vector+num does no validation pass - out-of-domain
  elements give NaN by contract (see api-contract input dispatch).
  Opts: {:log false :impl}."
  [x mean sd &opt opts]
  (default opts {})
  # kernel-backed: log-density = -0.5*z² - log(sd)
  # - 0.5*log(2π) composed from the vmap op set; resident input
  # (array stages, handle chains); vector+num does NO validation pass
  # (NaN domain contract - support is all-reals anyway); scalar+num
  # runs the identical pure formula
  (def log? (get opts :log))
  # the FULL-marker gate (the dlnorm marker
  # shape): a partial-marker struct is actionable misuse,
  # never a num/stage staging crash
  (when (and (or (struct? x) (table? x)) (not (packed-handle? x)))
    (error "stats: dnorm expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (if (and (= :num (num/resolve-impl opts)) (not (number? x)))
    (let [k (k*)
          z ((k 'vmap1-scalar) :div
             ((k 'vmap1-scalar) :sub (num/stage x k) mean) sd)
          z2 ((k 'vmap2) :mul z z)
          ld ((k 'vmap1-scalar) :add
              ((k 'vmap1-scalar) :mul z2 -0.5)
              (- (+ (math/log sd) (* 0.5 (math/log (* 2 math/pi))))))]
      (if log?
        ((k 'unpack-f64) ld)
        ((k 'unpack-f64) ((k 'vmap1) :exp ld))))
    (do
      (num/check-impl-landed opts x)
      (vapply
        (fn [v]
          (def z (/ (- v mean) sd))
          (def ld (- (+ (* 0.5 z z) (math/log sd) (* 0.5 (math/log (* 2 math/pi))))))
          (if log? ld (math/exp ld)))
        x))))

(defn pnorm
  "Normal CDF, erfc route (no cancellation in either tail); log
  with {:log true}. Kernel-backed on vectors; scalars run the
  identical pure formula. Opts: {:log false :impl}."
  [x mean sd &opt opts]
  (default opts {})
  (def log? (get opts :log))
  # the FULL-marker gate (the dlnorm marker
  # shape): a partial-marker struct is actionable misuse,
  # never a num/stage staging crash
  (when (and (or (struct? x) (table? x)) (not (packed-handle? x)))
    (error "stats: pnorm expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (if (and (= :num (num/resolve-impl opts)) (not (number? x)))
    (let [k (k*)
          out ((k 'unpack-f64)
               ((k 'pnorm-v) (num/stage x k) mean sd))]
      (if log? (map math/log out) out))
    (do
      (num/check-impl-landed opts x)
      (vapply
        (fn [v]
          (def z (/ (- v mean) sd))
          (def s (if (neg? z) -1 1))
          # erfc route: no cancellation in either tail
          (def p (* 0.5 (erfc (* s z (/ 1 (math/sqrt 2))))))
          (def v (if (neg? z) p (- 1 p)))
          (if log? (math/log v) v))
        x))))

(defn horner
  "Horner evaluation, coefficients high-to-low. Internal (AS 241
  aid)."
  [r & cs]
  # Horner evaluation, coefficients high-to-low (AS 241 port aid)
  (var acc (first cs))
  (each c (slice cs 1)
    (set acc (+ (* acc r) c)))
  acc)

(defn q-lower-tail
  "Read the q-family :lower-tail (default true) and reject
  :log-p loudly (no v1 lane - a raw p in [0,1] is the only
  supported form). Internal."
  [name opts]
  # :log-p was silently ignored by every q family - the p stayed
  # raw, so a legitimate R log.p call NaN-poisoned with no error.
  # Rejected loud instead (same class as the sibling-option rejects)
  (when (get opts :log-p)
    (errorf "stats: %s: :log-p is not supported in v1 - pass the raw p in [0,1]" name))
  (get opts :lower-tail true))

(defn qnorm
  "Normal quantile: the AS 241 port (Wichura 1988), decoded
  verbatim from R 4.5.0 nmath/qnorm.c - rational approximations
  accurate to ~1e-16 over the whole open interval. :lower-tail
  false computes the UPPER-tail quantile (dpq.h's 0.5-p+0.5
  complement - subnormal p survives). Returns NaN elementwise for p
  outside [0,1] (R's convention). Opts: {:lower-tail true :impl}."
  [p mean sd &opt opts]
  # AS 241 port (Wichura 1988), decoded
  # verbatim from R 4.5.0 nmath/qnorm.c (coefficients copied from
  # the C source) - rational approximations accurate to ~1e-16 over
  # the whole open interval, replacing the erfc-bracket solver
  # whose absolute tolerance could not see the far tails.
  (default opts {})
  (comp-only opts p)
(def lower-tail (q-lower-tail "qnorm" opts))
  (vapply
    (fn [v]
      (def pp (if lower-tail v (- 1 v)))
      # R_Q_P01_boundaries on the RAW
      # p with the flag (dpq.h) - the old (= pp 1) hit the materialized
      # 1-v == 1.0 for subnormal v and returned Inf
      (if (or (not= pp pp) (< pp 0) (> pp 1)) math/nan (do
      (cond
        (= v 0) (if lower-tail (- math/inf) math/inf)
        (= v 1) (if lower-tail math/inf (- math/inf))
        (do
          (def q (- pp 0.5))
          (def val
            (if (<= (math/abs q) 0.425)
              # central: |p - 0.5| <= .425 (r = .425^2 - q^2)
              (do
                (def r (- 0.180625 (* q q)))
                (* q
                   (/ (horner r 2509.0809287301226727 33430.575583588128105
                              67265.770927008700853 45921.953931549871457
                              13731.693765509461125 1971.5909503065514427
                              133.14166789178437745 3.387132872796366608)
                      (horner r 5226.495278852854561 28729.085735721942674
                              39307.89580009271061 21213.794301586595867
                              5394.1960214247511077 687.1870074920579083
                              42.313330701600911252 1))))
              # near a boundary: r = sqrt(-log(min(p,1-p)))
              (do
                # R qnorm.c:98 - lp = log of
            # the small side via R_DT_CIv/R_DT_Iv on the RAW p (the
            # 0.5-p+0.5 complement), never the materialized 1-pp:
            # subnormal upper survives (qnorm(5e-324,up) = 38.467,
            # R-exact; was +Inf)
            (def lp (math/log (if (> q 0)
                             (if lower-tail (+ (- 0.5 v) 0.5) v)
                             (if lower-tail v (+ (- 0.5 v) 0.5)))))
                (def r (math/sqrt (- lp)))
                (cond
                  # r <= 5 <==> min(p,1-p) >= exp(-25) ~ 1.39e-11
                  (<= r 5)
                  (do
                    (def r2 (- r 1.6))
                    (/ (horner r2 7.7454501427834140764e-4
                               0.0227238449892691845833 0.24178072517745061177
                               1.27045825245236838258 3.64784832476320460504
                               5.7694972214606914055 4.6303378461565452959
                               1.42343711074968357734)
                       (horner r2 1.05075007164441684324e-9
                               5.475938084995344946e-4 0.0151986665636164571966
                               0.14810397642748007459 0.68976733498510000455
                               1.6763848301838038494 2.05319162663775882187 1)))
                  # 5 < r <= 27 <==> min(p,1-p) >= exp(-729) ~ 2.5e-317
                  (<= r 27)
                  (do
                    (def r2 (- r 5))
                    (/ (horner r2 2.01033439929228813265e-7
                               2.71155556874348757815e-5 0.0012426609473880784386
                               0.026532189526576123093 0.29656057182850489123
                               1.7848265399172913358 5.4637849111641143699
                               6.6579046435011037772)
                       (horner r2 2.04426310338993978564e-15
                               1.4215117583164458887e-7 1.8463183175100546818e-5
                               7.868691311456132591e-4 0.0148753612908506148525
                               0.13692988092273580531 0.59983220655588793769 1)))
                  # r > 27: practically only log-scale p; the
                  # asymptotic ladder (Maechler 2022, xs_1..xs_5)
                  (>= r 6.4e8)
                  (* r (math/sqrt 2))
                  (do
                    (def s2 (* -2 lp))
                    (var x2 (- s2 (math/log (* 2 math/pi s2))))
                    (when (< r 36000)
                      (set x2 (- s2 (math/log (* 2 math/pi x2)) (/ 2 (+ 2 x2))))
                      (when (< r 840)
                        (set x2 (- s2 (math/log (* 2 math/pi x2))
                                    (* 2 (math/log1p (- (/ 1 (/ (- 4 x2) (+ 2 x2))))))))
                        (when (< r 109)
                          (set x2 (- s2 (math/log (* 2 math/pi x2))
                                      (* 2 (math/log1p (- (/ (- 1 (/ 1 (+ 4 x2)))
                                                             (+ 2 x2)))))))
                          (when (< r 55)
                            (set x2 (- s2 (math/log (* 2 math/pi x2))
                                        (* 2 (math/log1p
                                               (/ (- 1 (/ (- 5 (/ 9 (+ 8 x2))) (+ 6 x2)))
                                                  (+ 4 x2))))))))))
                    (math/sqrt x2))))))
          # R's sign flip applies ONLY in the boundary branches
          # (the central branch's q already carries the sign);
          # then the location-scale transform
          (def z (if (<= (math/abs q) 0.425) val (if (< q 0) (- val) val)))
          (+ mean (* sd z))))
      )))
    p))

#==------------------------------------------------------------------------==#
# t
#==------------------------------------------------------------------------==#

(defn dt
  "Student-t density (log with {:log true}). Kernel-backed on
  vectors; Opts: {:log false :impl}."
  [x df &opt opts]
  (default opts {})
  # kernel-backed: log-density = lgamma((ν+1)/2) - lgamma(ν/2)
  # - 0.5*log(νπ) - ((ν+1)/2)*log(1+t²/ν) - the per-element term via
  # vmap + :log (log(1+x) vs pure log1p: ~1e-16 association diff,
  # inside every conformance bar); vector+num: no validation pass
  (def log? (get opts :log))
  # the FULL-marker gate (the dlnorm marker
  # shape): a partial-marker struct is actionable misuse,
  # never a num/stage staging crash
  (when (and (or (struct? x) (table? x)) (not (packed-handle? x)))
    (error "stats: dt expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (if (and (= :num (num/resolve-impl opts)) (not (number? x)))
    (let [k (k*)
          xf (num/stage x k)
          t2 ((k 'vmap2) :mul xf xf)
          q ((k 'vmap1-scalar) :add ((k 'vmap1-scalar) :div t2 df) 1)
          lq ((k 'vmap1) :log q)
          ld ((k 'vmap1-scalar) :add
               ((k 'vmap1-scalar) :mul lq (- (/ (+ df 1) 2)))
               (- (lgamma (/ (+ df 1) 2))
                  (lgamma (/ df 2))
                  (* 0.5 (math/log (* df math/pi)))))]
      (if log?
        ((k 'unpack-f64) ld)
        ((k 'unpack-f64) ((k 'vmap1) :exp ld))))
    (do
      (num/check-impl-landed opts x)
      (vapply
        (fn [v]
          (def ld (+ (- (lgamma (/ (+ df 1) 2))
                        (lgamma (/ df 2))
                        (* 0.5 (math/log (* df math/pi))))
                    (* (/ (+ df 1) 2) (- (math/log1p (/ (* v v) df))))))
          (if log? ld (math/exp ld)))
        x))))

(defn/typed pnt-core {:args [:number :number :number] :ret :number}
  [t df ncp]
  "The noncentral t CDF - R's pnt.c verbatim (AS 243, Lenth 1989;
  the Guenther 1978 twin series over doubly-truncated incomplete
  betas), lower tail only. Internal to pt."
  # betas), decoded from the r-source mirror. Lower tail only.
  # C's LDOUBLE accumulators -> doubles: the errbd < 1e-12 gate
  # keeps both within ~1e-12 (the residual absorbs into the row tol).
  (when (<= df 0)
    (errorf "stats: pnt needs df > 0, got %q" df))
  (def t-fin? (and (< t math/inf) (> t (- math/inf))))
  (cond
    (not t-fin?) (if (neg? t) 0 1)
    (and (< t 0) (> ncp 40)) 0
    (let [del (if (>= t 0) ncp (- ncp))
          negdel (< t 0)
          tt (math/abs t)]
      (if (or (> df 4e5) (> (* del del) 1416.0))
        # A&S 26.7.10 escape (|del| > 37.63 or huge df): the shifted-
        # scaled normal approx, upper tail when the reflection ran
        (let [s (/ 1 (* 4 df))
              phi (pnorm (* tt (- 1 s)) del (math/sqrt (+ 1 (* tt tt 2 s))) {})]
          (if negdel (- 1 phi) phi))
        (let [x2 (* t t)
              x (/ x2 (+ x2 df))
              tnc
              (if (<= x 0)
                0
                (do
                  (def lambda (* del del))
                  (var p (* 0.5 (math/exp (* -0.5 lambda))))
                  # underflow unreachable post-escape (del <= 37.63)
                  (var q (* 0.7978845608028654 p del))
                  (var s (- 0.5 p))
                  (when (< s 1e-7)
                    (set s (* -0.5 (math/expm1 (* -0.5 lambda)))))
                  (var a 0.5)
                  (def b (* 0.5 df))
                  (def rxb (math/pow (/ df (+ x2 df)) b))
                  (def albeta (+ 0.5723649429247001
                                 (- (lgamma b) (lgamma (+ b 0.5)))))
                  (var xodd (beta-inc a b x))
                  (var godd (* 2 rxb (math/exp (- (* a (math/log x)) albeta))))
                  (def tnc0 (* b x))
                  (var xeven (if (< tnc0 2.2204460492503131e-16)
                               tnc0
                               (- 1 rxb)))
                  (var geven (* tnc0 rxb))
                  (var acc (+ (* p xodd) (* q xeven)))
                  (var it 0)
                  (while (< it 1000)
                    (++ it)
                    (set a (+ a 1))
                    (set xodd (- xodd godd))
                    (set xeven (- xeven geven))
                    (set godd (* godd x (/ (+ a b -1) a)))
                    (set geven (* geven x (/ (- (+ a b) 0.5) (+ a 0.5))))
                    (set p (* p (/ lambda (* 2 it))))
                    (set q (* q (/ lambda (+ (* 2 it) 1))))
                    (set acc (+ acc (* p xodd) (* q xeven)))
                    (set s (- s p))
                    (when (or (< s -1e-10)
                              (and (<= s 0) (> it 1))
                              (< (math/abs (* 2 s (- xodd godd))) 1e-12))
                      (break)))
                  acc))]
          (def fin (min (+ tnc (pnorm (- del) 0 1 {})) 1))
          (if negdel (- 1 fin) fin))))))

(defn pt
  "Student-t CDF: the central lane via the regularized incomplete
  beta, the NONCENTRAL lane (R's pnt.c verbatim, AS 243 - the
  Guenther twin series) with {:ncp}. Log with {:log true}.
  Opts: {:ncp nil :log false :impl}."
  [x df &opt opts]
  (default opts {})
  (comp-only opts x)
  (def log? (get opts :log))
  (def ncp (get opts :ncp))
  (vapply
    (fn [v]
      (def val
        (if (and ncp (not= ncp 0))
          (pnt-core v df ncp)
          (do
            # P(T > x) = 0.5 * I_{df/(df+x^2)}(df/2, 1/2) for x > 0 (upper tail)
            (def up (* 0.5 (beta-inc (/ df 2) 0.5 (/ df (+ df (* v v))))))
            # v < 0: P(T<=v) = P(T>|v|) = up; v > 0: 1 - up
            (if (neg? v) up (- 1 up)))))
      (if log? (math/log val) val))
    x))


#==------------------------------------------------------------------------==#
# chi-square
#==------------------------------------------------------------------------==#

(defn dchisq
  "Chi-squared density (log with {:log true}). Kernel-backed on
  vectors. Opts: {:log false :impl}."
  [x df &opt opts]
  (default opts {})
  # kernel-backed: log-density = (ν/2−1)·log x − x/2 −
  # (ν/2)·log2 − lgamma(ν/2). NaN domain contract: x < 0 → log → NaN
  # (the pure path throws; the vector path trusts - documented)
  (def log? (get opts :log))
  # the FULL-marker gate (the dlnorm marker
  # shape): a partial-marker struct is actionable misuse,
  # never a num/stage staging crash
  (when (and (or (struct? x) (table? x)) (not (packed-handle? x)))
    (error "stats: dchisq expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (if (and (= :num (num/resolve-impl opts)) (not (number? x)))
    (let [k (k*)
          xf (num/stage x k)
          lv ((k 'vmap1) :log xf)
          a ((k 'vmap1-scalar) :mul lv (- (/ df 2) 1))
          b ((k 'vmap1-scalar) :mul xf -0.5)
          ld ((k 'vmap1-scalar) :add ((k 'vmap2) :add a b)
              (- (+ (* (/ df 2) (math/log 2)) (lgamma (/ df 2)))))]
      (if log?
        ((k 'unpack-f64) ld)
        ((k 'unpack-f64) ((k 'vmap1) :exp ld))))
    (do
      (num/check-impl-landed opts x)
      (vapply
        (fn [v]
          # density poles to R - out of
          # domain or inf -> 0 (log: -inf)
          (if (or (< v 0) (= v math/inf) (= v (- math/inf)))
            (if log? (- math/inf) 0.0)
            (let [ld (- (* (- (/ df 2) 1) (math/log v))
                       (/ v 2)
                       (* (/ df 2) (math/log 2))
                       (lgamma (/ df 2)))]
              (if log? ld (math/exp ld)))))
        x))))

(defn pchisq
  "Chi-squared CDF = gamma-p(df/2, x/2) (log with {:log true});
  poles follow R (x=+inf -> 1, x<0 -> 0). Kernel-backed on
  vectors. Opts: {:log false :impl}."
  [x df &opt opts]
  (default opts {})
  (def log? (get opts :log))
  # the FULL-marker gate (the dlnorm marker
  # shape): a partial-marker struct is actionable misuse,
  # never a num/stage staging crash
  (when (and (or (struct? x) (table? x)) (not (packed-handle? x)))
    (error "stats: pchisq expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (if (and (= :num (num/resolve-impl opts)) (not (number? x)))
    (let [k (k*)
          n (num/count x)
          half-x ((k 'vmap1-scalar) :div (num/stage x k) 2)
          out ((k 'unpack-f64)
               ((k 'gamma-p-v) ((k 'pack-f64) (array/new-filled n (/ df 2))) half-x))]
      (if log? (map math/log out) out))
    (do
      (num/check-impl-landed opts x)
      (vapply
        (fn [v]
          # R poles: x=+inf -> 1, x<0 (incl -inf) -> 0
          (def p (cond (= v math/inf) 1.0
                       (< v 0) 0.0
                       (gamma-p (/ df 2) (/ v 2))))
          (if log? (math/log p) p))
        x))))


#==------------------------------------------------------------------------==#
# F
#==------------------------------------------------------------------------==#

(defn df-f
  "F density (R's df; exported as df-f because df reads as
  degrees-of-freedom elsewhere). Log with {:log true}.
  Kernel-backed on vectors. Opts: {:log false :impl}."
  [x df1 df2 &opt opts]
  (default opts {})
  # kernel-backed: log-density = −logB(d1/2,d2/2) +
  # (d1/2)·log(d1/d2) + (d1/2−1)·log x − ((d1+d2)/2)·log(1+x·d1/d2).
  # NaN domain contract: x < 0 → log → NaN (pure throws; vector trusts)
  (def log? (get opts :log))
  # the FULL-marker gate (the dlnorm marker
  # shape): a partial-marker struct is actionable misuse,
  # never a num/stage staging crash
  (when (and (or (struct? x) (table? x)) (not (packed-handle? x)))
    (error "stats: df-f expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (if (and (= :num (num/resolve-impl opts)) (not (number? x)))
    (let [k (k*)
          xf (num/stage x k)
          lv ((k 'vmap1) :log xf)
          a ((k 'vmap1-scalar) :mul lv (- (/ df1 2) 1))
          u ((k 'vmap1-scalar) :add ((k 'vmap1-scalar) :mul xf (/ df1 df2)) 1)
          lu ((k 'vmap1) :log u)
          b ((k 'vmap1-scalar) :mul lu (- (/ (+ df1 df2) 2)))
          ld ((k 'vmap1-scalar) :add ((k 'vmap2) :add a b)
              (+ (- (log-beta (/ df1 2) (/ df2 2)))
                 (* (/ df1 2) (math/log (/ df1 df2)))))]
      (if log?
        ((k 'unpack-f64) ld)
        ((k 'unpack-f64) ((k 'vmap1) :exp ld))))
    (do
      (num/check-impl-landed opts x)
      (vapply
        (fn [v]
          # poles to R - x < 0 or inf
          # -> 0 (log: -inf)
          (if (or (< v 0) (= v math/inf) (= v (- math/inf)))
            (if log? (- math/inf) 0.0)
            (let [ld (+ (- (log-beta (/ df1 2) (/ df2 2)))
                       (* (/ df1 2) (math/log (/ df1 df2)))
                       (* (- (/ df1 2) 1) (math/log v))
                       (* (- (/ (+ df1 df2) 2)) (math/log1p (* v (/ df1 df2)))))]
              (if log? ld (math/exp ld)))))
        x))))

(defn pf
  "F CDF = beta-inc(df1/2, df2/2, ...); poles follow R (x<0 -> 0,
  +inf -> 1); log with {:log true}. Opts: {:log false :impl}."
  [x df1 df2 &opt opts]
  (default opts {})
  (comp-only opts x)
  (def log? (get opts :log))
  (vapply
    (fn [v]
      # x < 0 (incl -inf) -> 0
      # the +inf pole: R's pbeta(1) limit is 1; the beta-inc route
      # yields inf/inf -> NaN (measured)
      (if (< v 0)
        (if log? (- math/inf) 0.0)
        (let [p (if (= v math/inf)
                 1
                 (beta-inc (/ df1 2) (/ df2 2)
                           (/ (* v df1) (+ (* v df1) df2))))]
          (if log? (math/log p) p))))
    x))


#==------------------------------------------------------------------------==#
# binomial
#==------------------------------------------------------------------------==#

(defn r-int-tol
  "R's non-integer tolerance (nmath R_NON_INT_TOL): round-half-up
  within the relative 1e-9 band, else nil. Internal."
  [v]
  [v]
  # R's non-integer tolerance (nmath R_NON_INT_TOL): m = round-half-
  # up (floor(x+0.5)), tolerance |x - m| <= 1e-9*max(1,|x|) -> m;
  # beyond -> nil (non-integer). Boundary-measured vs R 4.5.0 from
  # BOTH SIDES (R rounds - an earlier one-sided decode floored, and
  # the below band NaNed):
  # x=2: 2±2e-9, x=7: 7±6.9e-9, x=200: 200±1.9e-7. Kernel-ok gates
  # stay EXACT-integer: a tolerance-accepted non-exact element falls
  # to the elementwise path, which rounds here - no hot-path scan.
  (def m (math/floor (+ v 0.5)))
  (if (<= (math/abs (- v m)) (* 1e-9 (max 1 (math/abs v)))) m nil))

(defn dbinom
  "Binomial PMF (log with {:log true}). SUPPORT CONVENTION:
  near-integers within R's relative 1e-9 tolerance round to the
  integer; beyond the tolerance or out of support -> 0 (log:
  -inf; R's own value + warning). Kernel-backed on exact-integer
  arrays; any other element falls to the elementwise pure path.
  Opts: {:log false :impl}."
  [k n prob &opt opts]
  (default opts {})
  # kernel-backed. SUPPORT CONVENTION: near-integers within R's
  # relative 1e-9 tolerance ROUND
  # to the integer (r-int-tol); beyond the tolerance -> 0
  # (log: -inf) - R's own value + warning (an early draft NaNed
  # here on a mis-decode, reversed on live evidence);
  # out-of-support k -> 0 (log: -inf). Arrays with bad
  # elements fall to the elementwise pure path (fitness check)
  (def log? (get opts :log))
  (def kernel-ok
    (and (= :num (num/resolve-impl opts)) (indexed? k)
         (not= prob 0) (not= prob 1)
         (all |(and (= (math/round $) $)
                    (>= (math/round $) 0)
                    (<= (math/round $) n))
              k)))
  (if kernel-ok
    (do
      (let [acc (k*)
            kf (num/stage k acc)
            lg1 ((acc 'lgamma-v) ((acc 'vmap1-scalar) :add kf 1))
            nm ((acc 'vmap1) :neg ((acc 'vmap1-scalar) :sub kf n))
            lg2 ((acc 'lgamma-v) ((acc 'vmap1-scalar) :add nm 1))
            t1 ((acc 'vmap1-scalar) :mul kf (math/log prob))
            t2 ((acc 'vmap1-scalar) :mul nm (math/log1p (- prob)))
            g-lg1 ((acc 'vmap1) :neg
                   ((acc 'vmap1-scalar) :sub lg1 (lgamma (+ n 1))))
            s1 ((acc 'vmap2) :sub g-lg1 lg2)
            s2 ((acc 'vmap2) :add s1 t1)
            ld ((acc 'vmap2) :add s2 t2)]
        (if log?
          ((acc 'unpack-f64) ld)
          ((acc 'unpack-f64) ((acc 'vmap1) :exp ld)))))
    (do
      (comp-only opts k)
      (vapply
        (fn [v]
          # the relative tolerance decode
          (def kk (r-int-tol v))
          (cond
            # poles BEFORE the tolerance decode: ±inf/out-of-support
            # -> 0 (r-int-tol of inf is nil - not the non-integer class)
            (or (< v 0) (= v math/inf)) (if log? (- math/inf) 0)
            (nil? kk) (if log? (- math/inf) 0)
            (> kk n) (if log? (- math/inf) 0)
            (or (= kk 0) (= prob 0))
            (if (= kk 0) (if log? (* n (math/log1p (- prob))) (math/pow (- 1 prob) n)) (if log? (- math/inf) 0))
            (= kk n)
            (if (= prob 1) (if log? 0 1) (if log? (* n (math/log prob)) (math/pow prob n)))
            (do
              (def ld (+ (lgamma (+ n 1))
                         (- (lgamma (+ kk 1)))
                         (- (lgamma (+ (- n kk) 1)))
                         (* kk (math/log prob))
                         (* (- n kk) (math/log1p (- prob)))))
              (if log? ld (math/exp ld)))))
        k))))

(defn pbinom
  "Binomial CDF: the exact upward recurrence (R-exact on
  knife-edge dyadic p), the backward-tail recompute in the
  saturated upper region (R's beta-complement precision); R's
  floor(q + 1e-7) coercion, the v<0 pole before it. Log with
  {:log true}; poles follow R. Opts: {:log false :impl}."
  [k n prob &opt opts]
  (default opts {})
  # documented composition: the exact upward recurrence IS the
  # R-exact algorithm on knife-edge dyadic p (a beta-inc kernel path
  # would trade that exactness for speed - recorded decision)
  (comp-only opts k)
  (def log? (get opts :log))
  (vapply
    (fn [v]
      # R coerces with floor(q + 1e-7)
      # (nmath semantics, verified live: pbinom(2.5)=pbinom(2),
      # pbinom(2.9999999)=pbinom(3)) - math/round was an early draft's
      # class of error. The raw v<0 pole runs BEFORE the
      # coercion - the tiny-negative zone (-1e-10..-1.5e-7) must be
      # 0 (R), not P(X<=0) via floor(v+1e-7)=0
      (cond
        (< v 0) (if log? (- math/inf) 0.0)
        (do
          (def kk (math/floor (+ v 1e-7)))
          # R poles: k<0 -> 0; k>=n saturates at 1
          (def kmax (min kk n))
      (def praw
        (cond
          (< kk 0) 0.0
          (>= kk n) 1.0
          (= prob 0) (if (>= kk 0) 1.0 0.0)
          (= prob 1) (if (>= kk n) 1.0 0.0)
          # exact upward recurrence while terms are representable -
          # matches R on knife-edge p (dyadic probs sum exactly);
          # log-space only for the extreme regime
          (let [t0 (math/pow (- 1 prob) n)]
            (if (and (> t0 0) (< n 1500))
              (do
                (var sum t0)
                (var term t0)
                (def ratio (/ prob (- 1 prob)))
                (for i 0 kmax
                  (*= term (/ (* (- n i) ratio) (+ i 1)))
                  (+= sum term))
                sum)
              (do
                (def lterms (array/new (+ kmax 1)))
                (for i 0 (+ kmax 1)
                  (array/push lterms (dbinom i n prob {:log true})))
                (math/exp (logsumexp lterms)))))))
      # BOTH summation branches
      # round ABOVE 1 in the saturated upper region
      # (pbinom(960,1000,0.9) = 1.0000000000001048, an impossible CDF
      # value that derailed the upper walks by up to 8 counts).
      # Recompute from the BACKWARD tail (P(K=n) down to P(K=kmax+1))
      # and return 1 - tail - R's beta-complement precision
      (def p (if (> praw 1)
             (do
               (def ratio (/ prob (- 1 prob)))
               (var t2 (math/pow prob n))
               (var tail t2)
               (var j n)
               (while (> j (+ kmax 1))
                 # term_{j-1} = term_j * j / ((n-j+1)*ratio)
                 (set t2 (/ (* t2 j) (* (- (+ n 1) j) ratio)))
                 (-- j)
                 (+= tail t2))
               (- 1 tail))
             praw))
          (if log? (math/log p) p))))
    k))

(defn qbinom
  "Binomial quantile: R's q_DISCRETE_BODY search decoded
  (qDiscrete_search.h) - the 8eps fuzz on the raw p, the
  Cornish-Fisher seed on the flag-aware qnorm, do_search on the
  flag-threaded CDF. Returns NaN elementwise for p outside
  [0,1]. Opts: {:lower-tail true :impl}."
  [p n prob &opt opts]
  (default opts {})
  (comp-only opts p)
  # R's
  # q_DISCRETE_BODY search decoded (qDiscrete_search.h) - the 8eps
  # left-continuity fuzz on the RAW p (lower shrinks p, upper grows
  # p only when 1-p > 32eps = _xf_*8eps, the source-verified
  # constant; the old lane used 16eps and fuzzed cells R does not),
  # the Cornish-Fisher seed on the flag-aware qnorm (the old
  # materialized target fed qnorm(1.0) = inf at subnormal p), and
  # do_search on the FLAG-THREADED CDF: lower walks P(K<=k) =
  # pbinom; upper walks P(K>k) = I_prob(k+1, n-k) DIRECTLY (a
  # materialized 1-pbinom saturates in the far upper region where
  # pbinom's recurrence rounds above 1 - the off-by-8
  # class)
(def lower-tail (q-lower-tail "qbinom" opts))
  (vapply
    (fn [v]
      (if (or (not= v v) (< v 0) (> v 1)) math/nan (do
      (def e 1.7763568394002505e-15)
      (def pt (if lower-tail (* v (- 1 e))
                (if (> (- 1 v) (* 4 e)) (* v (+ 1 e)) v)))
      (defn pdist [ky]
        (if lower-tail
          (pbinom ky n prob {})
          # P(K>k) = I_prob(k+1, n-k); k >= n -> 0 (beta-inc needs
          # b > 0)
          (if (>= ky n) 0.0 (beta-inc (+ ky 1) (- n ky) prob))))
      (var k0 (math/round (+ (* n prob)
                (* (qnorm v 0 1 (if lower-tail {} {:lower-tail false}))
                   (math/sqrt (* n prob (- 1 prob)))))))
      (set k0 (max 0 (min n k0)))
      (def zc (pdist k0))
      (if (if lower-tail (>= zc pt) (< zc pt))
        # search left: step down while P(k-1) still clears the target
        (while (and (> k0 0)
                    (if lower-tail
                      (>= (pdist (- k0 1)) pt)
                      (< (pdist (- k0 1)) pt)))
          (-- k0))
        # search right: step up until P(k) clears the target
        (while (if lower-tail
                 (< (pdist k0) pt)
                 (>= (pdist k0) pt))
          (++ k0)))
      k0)))
    p))

#==------------------------------------------------------------------------==#
# poisson
#==------------------------------------------------------------------------==#

(defn dpois
  "Poisson PMF; the same support convention as dbinom, plus
  lambda=0 -> the point mass at 0 (log with {:log true}).
  Kernel-backed on exact-integer arrays. Opts: {:log false
  :impl}."
  [k lambda &opt opts]
  (default opts {})
  # kernel-backed. SUPPORT CONVENTION: same as
  # dbinom - non-integer k -> NaN, k < 0 -> 0 (log: -inf); bad
  # element arrays fall to the elementwise pure path
  (def log? (get opts :log))
  # the FULL-marker gate (the dlnorm marker
  # shape) - a partial-marker struct is actionable misuse
  # (the reproduced cell: dpois({:n 2}) died in the int prepass with
  # "can not convert nil nil to 64 bit signed integer")
  (when (and (or (struct? k) (table? k)) (not (packed-handle? k)))
    (error "stats: dpois expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  # lambda=0 -> the point mass at
  # 0 (dpois(0,0)=1, log 0; dpois(k>0,0)=0); +inf k -> 0. Bad
  # element arrays fall to the elementwise pure path
  (def kernel-ok
    (and (= :num (num/resolve-impl opts)) (not (number? k)) (> lambda 0)
         (or (struct? k)
             (and (indexed? k)
                  (all |(and (= (math/round $) $) (>= (math/round $) 0)
                             (not= $ math/inf))
                       k)))))
  (if kernel-ok
    (do
      (let [acc (k*)
            kf (num/stage k acc)
            lg ((acc 'lgamma-v) ((acc 'vmap1-scalar) :add kf 1))
            t ((acc 'vmap1-scalar) :mul kf (math/log lambda))
            ld ((acc 'vmap1-scalar) :add ((acc 'vmap2) :sub t lg) (- lambda))]
        (if log?
          ((acc 'unpack-f64) ld)
          ((acc 'unpack-f64) ((acc 'vmap1) :exp ld)))))
    (do
      (comp-only opts k)
      (vapply
        (fn [v]
          # the relative tolerance decode
          # (r-int-tol); beyond tolerance = NaN (OUR convention)
          (def kk (r-int-tol v))
          (cond
            # poles before the tolerance decode (inf is not the
            # non-integer class); beyond-tolerance -> 0
            (or (< v 0) (= v math/inf)) (if log? (- math/inf) 0)
            (nil? kk) (if log? (- math/inf) 0)
            (= lambda 0) (if (= kk 0) (if log? 0 1) (if log? (- math/inf) 0))
            (do
              (def ld (- (* kk (math/log lambda)) lambda (lgamma (+ kk 1))))
              (if log? ld (math/exp ld)))))
        k))))

(defn ppois
  "Poisson CDF = Q(k+1, lambda) (the regularized upper gamma);
  R's floor(q + 1e-7) coercion, poles follow R. Log with
  {:log true}. Opts: {:log false :impl}."
  [k lambda &opt opts]
  (default opts {})
  # R coerces with floor(q + 1e-7) -
  # the old math/round rounded 2.5 UP (P(<=3)); R floors (P(<=2)).
  # Non-integer/inf arrays fall to the elementwise pure path
  (def log? (get opts :log))
  # the FULL-marker gate (the dlnorm marker
  # shape) - a partial-marker struct is actionable misuse,
  # never a num/stage staging crash
  (when (and (or (struct? k) (table? k)) (not (packed-handle? k)))
    (error "stats: ppois expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (def kernel-ok
    (and (= :num (num/resolve-impl opts)) (not (number? k))
         (or (struct? k)
             (and (indexed? k)
                  (all |(and (= (math/round $) $) (>= (math/round $) 0)
                             (not= $ math/inf))
                       k)))))
  (if kernel-ok
    (let [kn (k*)
          n (num/count k)
          kk (num/stage k kn)
          a1 ((kn 'vmap1-scalar) :add kk 1)
          out ((kn 'unpack-f64)
               ((kn 'gamma-q-v) a1 ((kn 'pack-f64) (array/new-filled n lambda))))]
      (if log? (map math/log out) out))
    (do
      (num/check-impl-landed opts k)
      (vapply
        (fn [v]
          # raw v<0 BEFORE the coercion
          # (the tiny-negative zone is 0 per R, not P(X<=0))
          (cond
            (< v 0) (if log? (- math/inf) 0.0)
            (do
              (def kk (math/floor (+ v 1e-7)))
              # R poles: k<0 -> 0, +inf -> 1
              (def p (cond (= v math/inf) 1.0
                           # P(K <= k) = Q(k+1, lambda) - regularized UPPER gamma
                           (gamma-q (+ kk 1) lambda)))
              (if log? (math/log p) p))))
        k))))

(defn qpois
  "Poisson quantile: R's q_DISCRETE_BODY search (qDiscrete_search.h
  + qpois.c) - the 8eps fuzz, the Cornish-Fisher seed, the upper
  lane walking gamma-p directly (no materialized 1-p). Returns
  NaN elementwise for p outside [0,1]. Opts: {:lower-tail true
  :impl}."
  [p lambda &opt opts]
  (default opts {})
  (comp-only opts p)
  # R's q_DISCRETE_BODY
  # search decoded (qDiscrete_search.h + qpois.c) - the 8eps fuzz
  # on the RAW p (lower shrinks, upper grows only when 1-p > 32eps),
  # the Cornish-Fisher seed (z + (z^2-1)/(6*sigma)) on the flag-aware
  # qnorm, and do_search on the FLAG-THREADED CDF: lower walks
  # P(K<=y) = ppois vs p; upper walks P(K>y) = gamma-p(y+1, lambda)
  # vs the raw p DIRECTLY. The old complement walk's threshold
  # 1-v+8eps saturated ABOVE 1 at v < ~2e-15 - the walk never
  # advanced and qpois(1e-16,100) returned the normal-approx seed
  # (18 vs R 30); the upper lane's materialized target rounded to
  # 1.0 at subnormal p and walked to inf (qpois(1e-300,5,up) vs
  # R 238)
(def lower-tail (q-lower-tail "qpois" opts))
  (vapply
    (fn [v]
      (if (or (not= v v) (< v 0) (> v 1)) math/nan (do
      # R's boundaries (qpois.c:61-62): p==0 -> 0, p==1 -> +Inf,
      # flag-composed
      (cond
        (= v 0) (if lower-tail 0 math/inf)
        (= v 1) (if lower-tail math/inf 0)
        (= lambda 0) 0
        (do
          (def e 1.7763568394002505e-15)
          (def pt (if lower-tail (* v (- 1 e))
                    (if (> (- 1 v) (* 4 e)) (* v (+ 1 e)) v)))
          (def sg (math/sqrt lambda))
          (def z (qnorm v 0 1 (if lower-tail {} {:lower-tail false})))
          (var y (math/round (+ lambda
                               (* sg (+ z (/ (- (* z z) 1) (* 6 sg)))))))
          (when (< y 0) (set y 0))
          (defn pdist [ky]
            (if lower-tail (ppois ky lambda {}) (gamma-p (+ ky 1) lambda)))
          (def zc (pdist y))
          (if (if lower-tail (>= zc pt) (< zc pt))
            (while (and (> y 0)
                        (if lower-tail
                          (>= (pdist (- y 1)) pt)
                          (< (pdist (- y 1)) pt)))
              (-- y))
            (while (if lower-tail
                     (< (pdist y) pt)
                     (>= (pdist y) pt))
              (++ y)))
          y)))))
    p))

#==------------------------------------------------------------------------==#
# beta - plumbing; regularized incomplete beta (utils beta-inc =
# R pbeta lower tail, verified 15-digit). gamma below in RATE
# parameterisation (R default dgamma(x, shape, rate)).
#==------------------------------------------------------------------------==#

(defn dbeta
  "Beta density on (0,1) (log with {:log true}). Kernel-backed on
  vectors. Opts: {:log false :impl}."
  [x a b &opt opts]
  (default opts {})
  # kernel-backed: log-density = (a−1)·log x + (b−1)·log(1−x)
  # − logB(a,b). NaN domain contract: x outside (0,1) → NaN via log
  # (the pure path throws; the vector path trusts - documented)
  (def log? (get opts :log))
  # the FULL-marker gate (the dlnorm marker
  # shape): a partial-marker struct is actionable misuse,
  # never a num/stage staging crash
  (when (and (or (struct? x) (table? x)) (not (packed-handle? x)))
    (error "stats: dbeta expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (if (and (= :num (num/resolve-impl opts)) (not (number? x)))
    (let [k (k*)
          xf (num/stage x k)
          lv ((k 'vmap1) :log xf)
          t1 ((k 'vmap1-scalar) :mul lv (- a 1))
          lo ((k 'vmap1) :log ((k 'vmap1-scalar) :add ((k 'vmap1) :neg xf) 1))
          t2 ((k 'vmap1-scalar) :mul lo (- b 1))
          ld ((k 'vmap1-scalar) :add ((k 'vmap2) :add t1 t2) (- (log-beta a b)))]
      (if log?
        ((k 'unpack-f64) ld)
        ((k 'unpack-f64) ((k 'vmap1) :exp ld))))
    (do
      (num/check-impl-landed opts x)
      (vapply
        (fn [v]
          # poles to R - out of range
          # or inf -> 0 (log: -inf)
          (if (or (< v 0) (> v 1))
            (if log? (- math/inf) 0.0)
            (let [ld (- (+ (* (- a 1) (math/log v))
                          (* (- b 1) (math/log (- 1 v))))
                      (log-beta a b))]
              (if log? ld (math/exp ld)))))
        x))))

(defn pbeta
  "Beta CDF = the regularized incomplete beta I_x(a,b) (log with
  {:log true}); poles follow R. Kernel-backed on vectors.
  Opts: {:log false :impl}."
  [x a b &opt opts]
  (default opts {})
  # kernel-backed; resident input - arrays stage, handles chain
  (def log? (get opts :log))
  # the FULL-marker gate (the dlnorm marker
  # shape): a partial-marker struct is actionable misuse,
  # never a num/stage staging crash
  (when (and (or (struct? x) (table? x)) (not (packed-handle? x)))
    (error "stats: pbeta expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (if (and (= :num (num/resolve-impl opts)) (not (number? x)))
    (let [k (k*)
          n (num/count x)
          out ((k 'unpack-f64)
               ((k 'beta-inc-v) ((k 'pack-f64) (array/new-filled n a))
                                ((k 'pack-f64) (array/new-filled n b))
                                (num/stage x k)))]
      (if log? (map math/log out) out))
    (do
      (num/check-impl-landed opts x)
      (vapply
        (fn [v]
          # R poles: x<=0 -> 0 (incl -inf), x>=1 -> 1
          (def p (cond (<= v 0) 0.0 (>= v 1) 1.0 (beta-inc a b v)))
          (if log? (math/log p) p))
        x))))

(defn qbeta
  "Beta quantile: the mirrored solve + Newton polish (R qbeta.c's
  swap_tail in miniature) above pp 0.5, the DIRECT log-x Newton
  lane (R's use_log_x block) for far-tail targets - the
  documented extreme-(a,b) floor lives in api-contract. Returns
  NaN elementwise for p outside [0,1]; 0/1 endpoints exact.
  Opts: {:lower-tail true :impl}."
  [p a b &opt opts]
  (default opts {})
  (comp-only opts p)
(def lower-tail (q-lower-tail "qbeta" opts))
  (vapply
    (fn [v]
      (def pp (if lower-tail v (- 1 v)))
      (if (or (not= pp pp) (< pp 0) (> pp 1)) math/nan (do
      # R endpoints: exact 0/1, no solver noise
      (if (or (= pp 0) (= pp 1)) pp
        # the mirrored solve + Newton
        # polish, R qbeta.c's swap_tail (lines 190-212) in miniature.
        # At pp near 1 the CDF's double spacing (1.1e-16) EQUALS the
        # complement's own scale, so solve-monotone's comparisons are
        # staircase noise (a +-50% band: qbeta(1e-16,3.5,2,up) was 25%
        # off). Mirrored to the small side the CDF has full relative
        # resolution. tgt per dpq.h: the LOWER lane's complement via
        # the 0.5-pp+0.5 trick (R's CIv); the UPPER lane targets the
        # RAW v (R's CIv(p) = p when !lower_tail) - never the
        # double-rounded 1-(1-v). The polish is R's swap-side Newton
        # - bracket tol 1e-15 absolute is ~3e-7 RELATIVE at the far
        # tail without it
        (if (> pp 0.5)
          (let [tgt (if lower-tail (+ (- 0.5 pp) 0.5) v)]
            (var xc (solve-monotone (fn [t] (beta-inc b a t)) 0 1 tgt))
            (for _ 0 4
              (def fx (- (beta-inc b a xc) tgt))
              (def dfx (dbeta xc b a {}))
              (if (or (= dfx 0) (not= dfx dfx)) (break))
              (def nx (- xc (/ fx dfx)))
              (if (or (< nx 0) (> nx 1) (not= nx nx)) (break))
              (set xc nx))
            (- 1 xc))
          # the DIRECT lane's log-x Newton, R qbeta.c's
          # use_log_x block (lines 444-509) - solve-monotone's 1e-15
          # ABSOLUTE bracket width strands every root below ~1e-15
          # (qbeta(1e-12,0.5,0.5) was 8.88e-16 vs R 2.47e-24). u =
          # log(x); y = log I_x(a,b) computed DIRECTLY from
          # beta-inc's own log-space combination (a linear-scale I
          # underflows at the far tail); R's step w = (y-la)*exp(y-u
          # + logB + (1-a)u + (1-b)log1p(-x)); damped while u stays
          # < 0 (0 < x <= 1). Engaged when the leading-order seed
          # x0 = (tgt*a*B(a,b))^(1/a) is in (0, 1e-13) - R's own cut
          # is seed < e^-5; ours is conservative so every verified
          # mid cell keeps the bracket path. Seeds that underflow to
          # 0 (the a<~1e-4 medians, whose true root is below the
          # double grid) stay on the bracket - the documented floor.
          (let [tgt pp
                lb (log-beta a b)
                x0 (math/exp (/ (+ (math/log tgt) (math/log a) lb) a))]
            (if (and (> x0 0) (< x0 1e-13))
              (do
                (def la (math/log tgt))
                (def r (- 1 a))
                (def t (- 1 b))
                (var u (math/log x0))
                (var done false)
                (for _ 0 100
                  (when (not done)
                    (def x (math/exp u))
                    (def l1mx (math/log1p (- x)))
                    (def y (+ (- lb) (* a u) (* t l1mx)
                              (math/log (beta-cf a b x))
                              (- (math/log a))))
                    (def p1 (- y la))
                    (when (< (math/abs p1) (* 1e-15 (math/abs la)))
                      (set done true) (break))
                    (def w (* p1 (math/exp (+ (- y u) lb (* r u) (* t l1mx)))))
                    (when (or (not= w w) (= w math/inf) (= w (- math/inf)))
                      (set done true) (break))
                    # R's damped step: shrink g until u stays < 0
                    (var gn 1)
                    (var un (- u w))
                    (while (>= un 0)
                      (set gn (/ gn 2))
                      (set un (- u (* gn w)))
                      (when (< gn 1e-8) (break)))
                    (when (>= un 0) (set done true) (break))
                    # improvement check on the log-CDF scale
                    (def x2 (math/exp un))
                    (def y2 (+ (- lb) (* a un) (* t (math/log1p (- x2)))
                               (math/log (beta-cf a b x2))
                               (- (math/log a))))
                    (when (>= (math/abs (- y2 la)) (math/abs p1))
                      (set done true) (break))
                    (set u un)))
                (math/exp u))
              (do
                (var x (solve-monotone (fn [q] (beta-inc a b q)) 0 1 tgt))
                (for _ 0 4
                  (def fx (- (beta-inc a b x) tgt))
                  (def dfx (dbeta x a b {}))
                  (if (or (= dfx 0) (not= dfx dfx)) (break))
                  (def nx (- x (/ fx dfx)))
                  (if (or (< nx 0) (> nx 1) (not= nx nx)) (break))
                  (set x nx))
                x))))))))
    p))

#==------------------------------------------------------------------------==#
# gamma - RATE parameterisation: dgamma(x, shape, rate), pgamma lower
# = P(shape, rate*x) via the regularized incomplete gamma (utils).
#==------------------------------------------------------------------------==#

(defn dgamma
  "Gamma density in the RATE parameterisation (R dgamma(x, shape,
  rate); log with {:log true}). Kernel-backed on vectors.
  Opts: {:log false :impl}."
  [x shape rate &opt opts]
  (default opts {})
  # kernel-backed: log-density = shape·log(rate) + (shape−1)·
  # log x − rate·x − lgamma(shape). NaN domain contract: x < 0 → NaN
  (def log? (get opts :log))
  # the FULL-marker gate (the dlnorm marker
  # shape): a partial-marker struct is actionable misuse,
  # never a num/stage staging crash
  (when (and (or (struct? x) (table? x)) (not (packed-handle? x)))
    (error "stats: dgamma expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (if (and (= :num (num/resolve-impl opts)) (not (number? x)))
    (let [k (k*)
          xf (num/stage x k)
          lv ((k 'vmap1) :log xf)
          t ((k 'vmap1-scalar) :mul lv (- shape 1))
          rv ((k 'vmap1-scalar) :mul xf (- rate))
          ld ((k 'vmap1-scalar) :add ((k 'vmap2) :add t rv)
              (- (+ (lgamma shape) (* (- shape) (math/log rate)))))]
      (if log?
        ((k 'unpack-f64) ld)
        ((k 'unpack-f64) ((k 'vmap1) :exp ld))))
    (do
      (num/check-impl-landed opts x)
      (vapply
        (fn [v]
          # poles to R - negative or inf
          # -> 0 (log: -inf)
          (if (or (< v 0) (= v math/inf) (= v (- math/inf)))
            (if log? (- math/inf) 0.0)
            (let [ld (- (+ (* shape (math/log rate))
                          (* (- shape 1) (math/log v)))
                      (+ (* rate v) (lgamma shape)))]
              (if log? ld (math/exp ld)))))
        x))))

(defn pgamma
  "Gamma CDF = the regularized lower incomplete gamma P(shape,
  rate*x) (log with {:log true}); poles follow R. Kernel-backed on
  vectors. Opts: {:log false :impl}."
  [x shape rate &opt opts]
  (default opts {})
  (def log? (get opts :log))
  # the FULL-marker gate (the dlnorm marker
  # shape): a partial-marker struct is actionable misuse,
  # never a num/stage staging crash
  (when (and (or (struct? x) (table? x)) (not (packed-handle? x)))
    (error "stats: pgamma expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (if (and (= :num (num/resolve-impl opts)) (not (number? x)))
    (let [k (k*)
          n (num/count x)
          rx ((k 'vmap1-scalar) :mul (num/stage x k) rate)
          out ((k 'unpack-f64)
               ((k 'gamma-p-v) ((k 'pack-f64) (array/new-filled n shape)) rx))]
      (if log? (map math/log out) out))
    (do
      (num/check-impl-landed opts x)
      (vapply
        (fn [v]
          # R poles: x=+inf -> 1, x<0 (incl -inf) -> 0
          (def p (cond (= v math/inf) 1.0
                       (< v 0) 0.0
                       (gamma-p shape (* rate v))))
          (if log? (math/log p) p))
        x))))

(defn qchisq-appr
  "R nmath qgamma.c qchisq_appr, decoded verbatim: the AS 91
  Phase I starting approximation on the raw p + flag. Internal to
  qgamma."
  [p nu g lower-tail]
  # the AS 91 Phase I starting approximation. The
  # RAW p + the lower_tail flag per the source - R
  # threads the flags through R_DT_log/R_DT_Clog/qnorm
  # (qchisq_appr(p, nu, g, lower_tail, log_p)); the materialized
  # pp swamped to 1.0 on the upper lane below v ~ 5.5e-17 and the
  # WH seed degenerated (qnorm(1) = inf). Lower-lane arithmetic is
  # unchanged (log p = log pp, log1p(-p) = log(1-pp) there).
  (def alpha (/ nu 2))
  (def c (- alpha 1))
  # R_DT_log(p): lower_tail ? log(p) : log1p(-p)
  (def plog (if lower-tail (math/log p) (math/log1p (- p))))
  # R_DT_Clog(p): lower_tail ? log1p(-p) : log(p)
  (def pclog (if lower-tail (math/log1p (- p)) (math/log p)))
  (var ch 0)
  (cond
    # small chi-squared
    (< nu (* -1.24 plog))
    (do
      (def lgam1pa
        (if (< alpha 0.5) (lgamma (+ alpha 1)) (+ (math/log alpha) g)))
      (set ch (math/exp (+ (/ (+ lgam1pa plog) alpha)
                           (math/log 2)))))
    # Wilson-Hilferty
    (> nu 0.32)
    (do
      (def x (qnorm p 0 1 (if lower-tail {} {:lower-tail false})))
      (def p1 (/ 2 (* 9 nu)))
      (set ch (* nu (math/pow (+ 1 (- p1) (* x (math/sqrt p1))) 3))))
      # approximation for p tending to 1 (R: ch = -2*(Clog(p) -
      # c*log(0.5*ch) + g)) fires inside the same branch below
    # small nu: the C7-C10 iteration
    (do
      (set ch 0.4)
      (def a (+ pclog g (* c (math/log 2))))
      (var q 0)
      (while true
        (set q ch)
        (def p1 (/ 1 (* 1 (+ ch (* (+ 4.67 ch) ch)))))
        (def p2 (* ch (+ 6.73 (* ch (+ 6.66 ch)))))
        (def t (+ -0.5 (* (+ 4.67 (* 2 ch)) p1)
                    (- (/ (+ 6.73 (* ch (+ 13.32 (* 3 ch)))) p2))))
        (set ch (- ch (/ (- 1 (* (math/exp (+ a (* 0.5 ch))) p2 p1)) t)))
        (when (<= (math/abs (- q ch)) (* 1e-2 (math/abs ch))) (break)))))
  # the near-1 correction (R applies it when ch > 2.2*nu+6, inside
  # the WH branch; hoisted here so the upper lane keeps it too -
  # note: the g term is ADDED (the first port negated it)
  (when (and (> ch 0) (> ch (+ (* 2.2 nu) 6)) (> nu 0.32))
    (set ch (* -2 (+ (- pclog (* c (math/log (* 0.5 ch)))) g))))
  ch)

(defn qgamma-halley
  "AS 91 Phase II (R nmath/qgamma.c): the Halley iteration on the
  chisq scale. Returns [ch max-newton]. Internal to qgamma."
  [pp shape g ch]
  (var ch2 ch)
  (var max-newton 1)
  (def c (- shape 1))
  (def s6 (/ (+ 120 (* c (+ 346 (* 127 c)))) 5040))
  (var q 0)
  (var i 0)
  (while (< i 1000)
    (++ i)
    (set q ch2)
    (def p1 (* 0.5 ch2))
    (def p2 (- pp (gamma-p shape p1)))
    (if (or (not= p2 p2) (<= ch2 0))
      (do (set ch2 q) (set max-newton 27) (break))
      (do
        (def t (* p2 (math/exp (+ (* shape (math/log 2)) g p1
                                  (- (* c (math/log ch2)))))))
        (def b (/ t ch2))
        (def a (- (* 0.5 t) (* b c)))
        (def s1 (/ (+ 210 (* a (+ 140 (* a (+ 105 (* a (+ 84 (* a (+ 70 (* 60 a)))))))))) 420))
        (def s2 (/ (+ 420 (* a (+ 735 (* a (+ 966 (* a (+ 1141 (* 1278 a)))))))) 2520))
        (def s3 (/ (+ 210 (* a (+ 462 (* a (+ 707 (* 932 a)))))) 2520))
        (def s4 (/ (+ 252 (* a (+ 672 (* 1182 a))) (* c (+ 294 (* a (+ 889 (* 1740 a)))))) 5040))
        (def s5 (/ (+ 84 (* 2264 a) (* c (+ 1175 (* 606 a)))) 2520))
        (set ch2 (+ ch2 (* t (+ 1 (* 0.5 t s1)
                              (- (* b c (+ s1 (- (* b (+ s2 (- (* b (+ s3 (- (* b (+ s4 (- (* b (- s5 (* b s6))))))))))))))))))))
        (when (< (math/abs (- q ch2)) (* 5e-7 ch2)) (break))
        (when (> (math/abs (- q ch2)) (* 0.1 ch2))
          (set ch2 (if (< ch2 q) (* 0.9 q) (* 1.1 q)))))))
  [ch2 max-newton])

(defn qgamma-newton
  "AS 91 final Newton (R nmath/qgamma.c END block) - the log-scale
  port with the complement discipline (pp > 0.5 via log1p(-Q);
  the upper lane targets the raw v). Internal to qgamma."
  [pp shape rate x max-it v lower-tail]
  # the linear complement was a
  # documented bridge, losing 6-25% in the 1e-14..1e-16 band; R
  # always converts to log scale - p1 = logP - logp differences
  # stay full-precision in the far tails). pp > 0.5: logP via
  # log1p(-Q) (the complement identity, full precision); pp <= 0.5:
  # logP via log(P) linear (P > 1e-300 representable in the
  # acceptance band). delta-x = p1*exp(p_ - g) (R's log form);
  # t = x - delta (lower-tail, R's branch). The
  # negative-step guard ORDER - t <= 0 branches BEFORE
  # any gamma-q/gamma-p call (their x>=0 VALIDATION threw before
  # the no-improvement guard could fire); damped half-step.
  # The upper lane targets the
  # RAW v - R's END block takes log(p) on the raw upper prob and
  # calls pgamma(x, ..., lower_tail, log) flag-threaded, so the
  # solver never sees the double-rounded 1-(1-v); the materialized
  # pp had silently biased the lane (rel 2.3e-3 at p=1e-16) and
  # swallowed p < ~1.1e-16 to inf (1-p == 1.0). Step x + delta per
  # R's `t = lower_tail ? x - t : x + t`; p_ = log(gamma-q) direct
  # (the gcf upper tail), mirroring pgamma's flag threading.
  (var x2 x)
  (var i 0)
  (while (< i max-it)
    (++ i)
    (if (not lower-tail)
      (do
        (def target-log (math/log v))
        (var p-log2 (math/log (gamma-q shape (* rate x2))))
        (when (= p-log2 (- math/inf)) (set x2 0) (break))
        (def p1 (- p-log2 target-log))
        (when (< (math/abs p1) (* 1e-15 (math/abs target-log))) (break))
        (def g (dgamma x2 shape rate {:log true}))
        (when (or (= g (- math/inf)) (= g 0)) (break))
        (var t (+ x2 (* p1 (math/exp (- p-log2 g)))))
        (when (<= t 0) (set t (/ x2 2)))
        (def p2-log (math/log (gamma-q shape (* rate t))))
        (when (or (> (math/abs (- p2-log target-log)) (math/abs p1))
                  (and (> i 1) (= (math/abs (- p2-log target-log)) (math/abs p1))))
          (break))
        (set x2 t))
      (if (> pp 0.5)
       (do
        # the log1p(-Q) route is GATED to the near-1 far
        # tail it was ported for (the 1-1e-12..1-1e-16 band);
        # pp in (0.5, 1-1e-6] rides the linear complement -
        # at pp just above 0.5 the log form iterates toward gamma-q's
        # series-truncation noise at huge shape (qgamma(0.5+1e-10,
        # 2e5, 0.5) was rel 7.01e-5; the old guard-rescue held 2e-12)
        (if (> pp (- 1 1e-6))
         (do
           (def target (- 1 pp))
           (def target-log (math/log1p (- target)))
           (var p-log2 (math/log1p (- (gamma-q shape (* rate x2)))))
           (def p1 (- p-log2 target-log))
           (when (< (math/abs p1) (* 1e-15 (math/abs target-log))) (break))
           (def g (dgamma x2 shape rate {:log true}))
           (when (or (= g (- math/inf)) (= g 0)) (break))
           (var t (- x2 (* p1 (math/exp (- p-log2 g)))))
           (when (<= t 0)
             (set t (/ x2 2))
             (when (< x2 1e-300) (break)))
           (def p2-log (math/log1p (- (gamma-q shape (* rate t)))))
           (when (or (> (math/abs (- p2-log target-log)) (math/abs p1))
                     (and (> i 1) (= (math/abs (- p2-log target-log)) (math/abs p1))))
             (break))
           (set x2 t))
         (do
           # the linear complement (the middle
           # domain, R's non-log_p Newton form)
           (def target (- 1 pp))
           (def qc (- (gamma-q shape (* rate x2)) target))
           (when (< (math/abs qc) (* 1e-15 target)) (break))
           (def dq (* -1 rate (dgamma x2 shape rate {})))
           (when (= dq 0) (break))
           (def t (- x2 (/ qc dq)))
           (def qc2 (- (gamma-q shape (* rate t)) target))
           (when (or (> (math/abs qc2) (math/abs qc))
                     (and (> i 1) (= (math/abs qc2) (math/abs qc))))
             (break))
           (set x2 t))))
      (do
        # the log-scale route is
        # GATED to the far-tail domain it was ported for (the
        # 1e-12..1e-16 band, where the linear complement lost
        # 6-25%); at pp in [1e-6, 0.5] - including the p=0.5 median -
        # the LINEAR branch runs instead: at large shape the
        # log form's early break strands the iterate on gamma-p's
        # lgamma-difference noise (qgamma(0.5,2e5,0.5) was rel 7.01e-5
        # where the linear branch held 7.79e-13; the border sits exactly at 0.5)
        (if (< pp 1e-6)
          (do
            (def pp-log (math/log pp))
            (def pv (gamma-p shape (* rate x2)))
            (if (<= pv 0) (break))
            (var p-log2 (math/log pv))
            (def p1 (- p-log2 pp-log))
            (when (< (math/abs p1) (* 1e-15 (math/abs pp-log))) (break))
            (def g (dgamma x2 shape rate {:log true}))
            (when (or (= g (- math/inf)) (= g 0)) (break))
            (var t (- x2 (* p1 (math/exp (- p-log2 g)))))
            (when (<= t 0)
              (set t (/ x2 2))
              (when (< x2 1e-300) (break)))
            (def pv2 (gamma-p shape (* rate t)))
            (if (<= pv2 0) (break))
            (def p2-log (math/log pv2))
            (when (or (> (math/abs (- p2-log pp-log)) (math/abs p1))
                      (and (> i 1) (= (math/abs (- p2-log pp-log)) (math/abs p1))))
              (break))
            (set x2 t))
          (do
            # the linear branch - the
            # median/interior domain, R's non-log_p Newton form
            (def pc (- (gamma-p shape (* rate x2)) pp))
            (when (< (math/abs pc) (* 1e-15 pp)) (break))
            (def dp (* rate (dgamma x2 shape rate {})))
            (when (= dp 0) (break))
            (def t (- x2 (/ pc dp)))
            (def pc2 (- (gamma-p shape (* rate t)) pp))
            (when (or (> (math/abs pc2) (math/abs pc))
                      (and (> i 1) (= (math/abs pc2) (math/abs pc))))
              (break))
            (set x2 t)))))))
  x2)

(defn qgamma
  "Gamma quantile: the AS 91 port (Best & Roberts, decoded from
  R 4.5.0 nmath/qgamma.c) - Phase I approximation, Phase II
  Halley, final Newton with the complement discipline. The upper
  lane targets the RAW p. At p=0.5 from shape ~4e6 the solver
  diverges: a scoped huge-shape error (the documented boundary -
  see api-contract tail precision). Returns NaN elementwise for
  p outside [0,1]. Opts: {:lower-tail true :impl}."
  [p shape rate &opt opts]
  # the AS 91 port (Best & Roberts,
  # decoded from R 4.5.0 nmath/qgamma.c) - Phase I approximation
  # (qchisq-appr), Phase II Halley (qgamma-halley), final Newton
  # with the complement discipline (qgamma-newton). Replaces the
  # Wilson-Hilferty + Newton + bisection whose ok-flag bypassed
  # the safety net (qgamma(1-1e-12,2,1) was 4.99e35 vs R 31.1).
  (default opts {})
  (comp-only opts p)
(def lower-tail (q-lower-tail "qgamma" opts))
  (vapply
    (fn [v]
      (def pp (if lower-tail v (- 1 v)))
      (if (or (not= v v) (< v 0) (> v 1)) math/nan (do
        # the endpoints key on the RAW v per R's
        # R_Q_P01_boundaries (p==0 -> lower_tail ? LEFT : RIGHT) -
        # the materialized pp hit 1.0 on the upper lane below
        # v ~ 5.5e-17 and answered inf for finite R values
        (cond
          (= v 0) (if lower-tail 0.0 math/inf)
          (= v 1) (if lower-tail math/inf 0.0)
          (= shape 0) 0
          (do
            (def g (lgamma shape))
            # Phase I - raw p + flag (the S1 threading)
            (def ch0 (qchisq-appr v (* 2 shape) g lower-tail))
            # Phase II or straight to Newton (R's skip conditions)
            (var ch ch0)
            (var max-newton 1)
            (cond
              (or (not= ch ch) (= ch math/inf)) (set max-newton 0)
              (< ch 5e-7) (set max-newton 20)
              # R's pMAX skip keys on the materialized p_
              # (unreliable near 1). On the upper lane the materialized
              # pp = 1-v carries ulp(1) absolute error - rel 1e-4 at
              # v=1e-12 - so the skip (Newton 20 from the seed, no
              # Phase II on the rounded complement) fires on the RAW
              # v below 1e-12; R's own boundary (v < 1e-14, pMAX
              # materialized) is ulp-arbitrary at the knife edge
              (or (> pp (- 1 1e-14)) (< pp 1e-100)
                  (and (not lower-tail) (<= v 1e-12))) (set max-newton 20)
              (let [r (qgamma-halley pp shape g ch)]
                (set ch (get r 0))
                (set max-newton (get r 1))))
            # final Newton on the gamma scale: x = 0.5*ch/rate
            (var x (/ (* 0.5 ch) rate))
            (when (> max-newton 0)
              (set x (qgamma-newton pp shape rate x max-newton v lower-tail)))
            # a diverged solver (NaN/negative x at huge shape) is a
            # LOUD scoped error, never a laundered confident 0 - the
            # old clamp returned 0 for finite R values (qchisq(0.5,
            # 8e6+), qf(0.025,1e9,1e9)). Mechanism: gamma-p's series
            # needs > 1000 iterations from shape ~2e4 at x=shape; R
            # routes the x~alph region to ppois_asymp/bd0 (pgamma.c
            # decoded) - our series-only port is the recorded
            # boundary, api-contract tail-precision (b)
            (if (or (not= x x) (< x 0))
              (errorf "stats: qgamma: huge-shape boundary - the solver diverged at shape %q (pgamma's series truncates near x=shape; R uses its bd0/stirlerr asymptotics there); documented in api-contract tail-precision (b)" shape)
              x))))))
    p))


(defn qchisq
  "Chi-squared quantile: composes qgamma(p, df/2, scale 2) with
  the same flags (R qchisq.c). Returns NaN elementwise for p
  outside [0,1]; 0/1 endpoints exact. Opts: {:lower-tail true
  :impl}."
  [p df &opt opts]
  (default opts {})
  (comp-only opts p)
  # the flip - R's qchisq composes
  # qgamma(p, df/2, scale 2) with the SAME flags (qchisq.c)
(def lower-tail (q-lower-tail "qchisq" opts))
  (vapply
    (fn [v]
      # R endpoints: p=0 -> 0, p=1 -> +inf.
      # NaN-p + out-of-domain p ->
      # NaN elementwise (was a whole-vector THROW class)
      (if (or (not= v v) (< v 0) (> v 1))
        math/nan
        (do
          (cond (= v 0) (if lower-tail 0.0 math/inf)
                (= v 1) (if lower-tail math/inf 0.0)
                # qchisq = qgamma(p, df/2,
                # scale 2) composed over the ported AS 91 (verified
                # exact vs live R: 58.9198 at 1-1e-12, df 3) - the
                # plateaued bracket solver is gone
                (qgamma v (/ df 2) 0.5
                        (if lower-tail {} {:lower-tail false}))))))
    p))

(defn qf
  "F quantile: R qf.c decoded - the extreme-df fudge past
  4e5, else the (1/qbeta(p, d2/2, d1/2, upper) - 1)*d2/d1
  composition (the upper-tail flip solves the SMALL beta side).
  Returns NaN elementwise for p outside [0,1]. Opts:
  {:lower-tail true :impl}."
  [p df1 df2 &opt opts]
  (default opts {})
  (comp-only opts p)
  # the flip - composed per-lane below
(def lower-tail (q-lower-tail "qf" opts))
  (vapply
    (fn [v]
      # endpoints exact (R qf(0)=0,
      # qf(1)=Inf) - no bracket saturation.
      # NaN-p + out-of-domain p ->
      # NaN elementwise (was bracket-endpoint garbage)
      (if (or (not= v v) (< v 0) (> v 1))
        math/nan
        (do
          (cond (= v 0) (if lower-tail 0.0 math/inf)
                (= v 1) (if lower-tail math/inf 0.0)
                # R qf.c's extreme-df fudge,
                # decoded from the source (lines 43-51) - "qbeta
                # doesn't do this well" past df 4e5. df1==Inf (so
                # df2==Inf too) -> 1. The
                # flags now compose per the source - lane 1 passes
                # qf's OWN lower_tail to qchisq, lane 2 INVERTS it
                # (qchisq(p, df2, !lower_tail)), the qbeta lane
                # INVERTS it (qbeta(p, d2/2, d1/2, !lower_tail))
                (and (<= df1 df2) (> df2 4e5))
                (if (= df1 math/inf) 1.0
                  (/ (qchisq v df1 (if lower-tail {} {:lower-tail false})) df1))
                (> df1 4e5)
                (/ df2 (qchisq v df2 (if lower-tail {:lower-tail false} {})))
                # R's qf composition -
                # (1/qbeta(p, d2/2, d1/2, upper) - 1)*d2/d1 - the
                # upper-tail flip solves the SMALL beta side (no
                # CDF-rounds-to-1 plateau; verified vs live R:
                # 97300.45 at 1-1e-16, df 4,7)
                (* (- (/ (qbeta v (/ df2 2) (/ df1 2)
                             (if lower-tail {:lower-tail false} {}))) 1)
                   (/ df2 df1))))))
    p))

(defn qt
  "Student-t quantile: the AS 396 port (Hill 1970, R qt.c) - the
  df~2/df~1 closed forms, the general P>P0 split + 2-term Taylor
  polish (skipped at df >= 1e7, a recorded boundary). Returns NaN
  elementwise for p outside [0,1]. Opts: {:lower-tail true
  :impl}."
  [p df &opt opts]
  (default opts {})
  (comp-only opts p)
  # the family-wide :lower-tail flip
(def lower-tail (q-lower-tail "qt" opts))
  (vapply
    (fn [v]
      # NaN-p + out-of-domain p ->
      # NaN elementwise (R's ISNAN check; janet comparisons are
      # false on NaN, so the old guard let it through)
      (if (or (not= v v) (< v 0) (> v 1))
        math/nan
        (cond
          # R_Q_P01_boundaries (qt.c:55)
          # with the tail flag
          (= v 0) (if lower-tail (- math/inf) math/inf)
          (= v 1) (if lower-tail math/inf (- math/inf))
          # R qt.c:100: ndf > 1e20 -> qnorm(p, flags) verbatim
          (> df 1e20) (qnorm v 0 1 (if lower-tail {} {:lower-tail false}))
          (do
            # R qt.c's neg/P decode
            # (lines 102-110, dpq.h's 0.5-p+0.5 complement): P =
            # 2*min(p', 1-p') without materializing 1-p on the small
            # side - subnormal p survives (qt(5e-324,5) =
            # -7.189486e+64, R-exact; was -Inf: the old v<0.5
            # recursion's (- 1 5e-324) swallowed to 1.0)
            (def cval (+ (- 0.5 v) 0.5))
            (def neg (and (or (not lower-tail) (< v 0.5))
                          (or lower-tail (> v 0.5))))
            (def Pp (if neg
                     (* 2 (if lower-tail v cval))
                     (* 2 (if lower-tail cval v))))
            (def q
              (cond
                # df ~= 2: R's closed forms (qt.c:112-127) as direct q
                (< (math/abs (- df 2)) 1e-12)
                (cond
                  (<= Pp 2.2250738585072014e-308) math/inf
                  (< (* 3 Pp) 2.2204460492503131e-16)
                  (/ 1 (math/sqrt Pp))
                  (> Pp 0.9)
                  (* (- 1 Pp) (math/sqrt (/ 2 (* Pp (- 2 Pp)))))
                  (math/sqrt (- (/ 2 (* Pp (- 2 Pp))) 2)))
                # df ~= 1 (df < 1 rides here too - R bisects, a
                # recorded deviation): Cauchy q = 1/tan(pi*P/2)
                (< df (+ 1 1e-12))
                (cond
                  (= Pp 1) 0
                  (> Pp 0) (/ 1 (math/tan (* math/pi (/ Pp 2))))
                  math/inf)
                # general df: the AS 396
                # port (Hill 1970, qt.c's usual case) - the P>P0
                # asymptotic-inverse split + the middle branch + the
                # 2-term Taylor polish (dt/pt iterations); replaced
                # the t^2=F(1,df) composition (215x garbage at 1e20)
                (let [ndf df
                      a (/ 1 (- ndf 0.5))
                      b (/ 48 (* a a))
                      # c = ((20700*a/b - 98)*a - 16)*a + 96.36
                      c (+ (* (- (* (- (/ (* 20700 a) b) 98) a) 16) a) 96.36)
                      # d = ((94.5/(b+c) - 3)/b + 1)*sqrt(a*pi/2)*ndf
                      d (* (+ (/ (- (/ 94.5 (+ b c)) 3) b) 1)
                           (math/sqrt (* a (/ math/pi 2))) ndf)]
                  # R's P_ok route
                  # (qt.c:146-157) - (d*P)^(2/df) QUANTIZES on the
                  # subnormal grid when d*P is subnormal (p <~ 1e-308);
                  # the log route (x = (log d + ln2 + log(P/2))/df,
                  # y0 = exp(2x)) keeps full tail precision.
                  # half-log = log(P/2) per dpq.h: is_neg_lower ?
                  # log(p) : log1p(-p) - the half two-sided mass
                  (def half-log (if (= lower-tail neg)
                                (math/log v)
                                (math/log1p (- v))))
                  (var y0 (math/pow (* d Pp) (/ 2.0 ndf)))
                  (var pok true)
                  (var xlog 0)
                  (when (< y0 2.2204460492503131e-16)
                    (set pok false)
                    (set xlog (/ (+ (+ (math/log d) 0.69314718055994531)
                                     half-log) ndf))
                    (set y0 (math/exp (* 2 xlog))))
                  (var q 0)
                  (cond
                    # P > P0(df): asymptotic inverse about normal
                    (or (and (< ndf 2.1) (> Pp 0.5)) (> y0 (+ 0.05 a)))
                    (do
                      (def x (qnorm (* 0.5 Pp) 0 1 {}))
                      (var cc c)
                      (when (< ndf 5)
                        (set cc (+ cc (* 0.3 (- ndf 4.5) (+ x 0.6)))))
                      # c = (((0.05*d*x - 5)*x - 7)*x - 2)*x + b + c
                      (set cc (+ (* (- (* (- (* (- (* 0.05 d x) 5) x) 7) x) 2) x) b cc))
                      (def y2 (* x x))
                      # y = (((((0.4y+6.3)y+36)y+94.5)/c - y - 3)/b + 1)*x
                      (def poly (+ (* (+ (* (+ (* 0.4 y2) 6.3) y2) 36) y2) 94.5))
                      (def y4 (* (+ (/ (- (- (/ poly cc) y2) 3) b) 1) x))
                      (def y5 (math/expm1 (* a y4 y4)))
                      (set q (math/sqrt (* ndf y5))))
                    # R qt.c:174-177: the underflow branch - !P_ok
                    # and x < -ln2*53: q = sqrt(ndf)*exp(-x)
                    (and (not pok) (< xlog -36.736800569677102))
                    (set q (* (math/sqrt ndf) (math/exp (- xlog))))
                    # middle: re-use y0
                    (do
                      (def denom (* (- (/ (+ ndf 6) (* ndf y0)) (+ (* 0.089 d) 0.822))
                                    (* (+ ndf 2) 3)))
                      (def y5 (+ (* (- (* (+ (/ 1 denom) (/ 0.5 (+ ndf 4))) y0) 1)
                                   (/ (+ ndf 1) (+ ndf 2)))
                                (/ 1 y0)))
                      (set q (math/sqrt (* ndf y5)))))
                  # the 2-term Taylor polish (Hill 1981; R's it<10
                  # loop). A recorded deviation: R polishes at every
                  # df, but our dt (lgamma
                  # difference) and the polish's upper tail (beta-inc
                  # near-1 complement, a measured deviation class) hold R-precision
                  # only to df~1e7; beyond that the polish corrupts
                  # the seed, which alone tracks R to <1e-10 at the
                  # bar cells (verified: rel 2.5e-9 at 1e8, 2.5e-11
                  # at 1e10, ulp at 1e16+).
                  (var it 0)
                  (var go true)
                  (while (and go (< it 10) (< ndf 1e7))
                    (++ it)
                    (def dy (dt q ndf {}))
                    (when (<= dy 0) (set go false) (break))
                    # R uses pt(q, ndf, upper=TRUE): P(T>q) =
                    # 0.5*I_{df/(df+q^2)}(df/2, 1/2), computed direct
                    # (1 - lower) loses the tail to cancellation
                    (def up (* 0.5 (beta-inc (/ ndf 2) 0.5 (/ ndf (+ ndf (* q q))))))
                    (def x (/ (- up (/ Pp 2)) dy))
                    (when (or (not= x x)
                              (<= (math/abs x) (* 1e-14 (math/abs q))))
                      (set go false) (break))
                    (def F (if (< (math/abs q) (- (math/sqrt (/ 1.7976931348623157e308 2)) ndf))
                             (* q (/ (+ ndf 1) (* 2 (+ (* q q) ndf))))
                             (/ (+ ndf 1) (* 2 (+ q (/ ndf q))))))
                    (def del (* x (+ 1 (* x F))))
                    (cond
                      (and (not= del math/inf) (not= (+ q del) math/inf) (not= (+ q del) (- math/inf)))
                      (set q (+ q del))
                      (and (not= x math/inf) (not= (+ q x) math/inf) (not= (+ q x) (- math/inf)))
                      (set q (+ q x))
                      (do (set go false))))
                  q)))
            (if neg (- q) q)))))
    p))


(defn dexp
  "Exponential density (RATE parameterisation; log with
  {:log true}); x < 0 -> 0. Kernel lane stages via the pure
  elementwise path when any element is out of domain. Opts:
  {:log false :impl}."
  [x rate &opt opts]
  (default opts {})
  # R dexp (rate interface; nmath computes with scale = 1/rate):
  # exp(-x/scale)/scale, log form -x/scale - log(scale)
  (def log? (get opts :log))
  # the FULL-marker gate (the dlnorm marker
  # shape) - a partial-marker struct is actionable misuse
  # (the reproduced cell: dexp({:a 1}) died staging with "could not
  # find method :* for 8 or :r* for nil")
  (when (and (or (struct? x) (table? x)) (not (packed-handle? x)))
    (error "stats: dexp expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (when (<= rate 0) (errorf "stats: dexp rate > 0 required, got %q" rate))
  (def scale (/ 1 rate))
  (if (and (= :num (num/resolve-impl opts)) (not (number? x)))
    (do
      # pole to R - x < 0 -> 0 (log: -inf);
      # the kernel lane stages via the pure elementwise path when any
      # element is out of domain
      (def pure-fallback (some |(< $ 0) x))
      (if pure-fallback
        (vapply
          (fn [v]
            (if (< v 0)
              (if log? (- math/inf) 0.0)
              (let [lq (- (/ v scale))]
                (if log? (- lq (math/log scale)) (/ (math/exp lq) scale)))))
          x)
        (let [k (k*)
              xf (num/stage x k)
              neg-q ((k 'vmap1) :neg ((k 'vmap1-scalar) :div xf scale))]
          (if log?
            ((k 'unpack-f64) ((k 'vmap1-scalar) :add neg-q (- (math/log scale))))
            ((k 'unpack-f64) ((k 'vmap1-scalar) :div ((k 'vmap1) :exp neg-q) scale))))))
    (do
      (comp-only opts x)
      (vapply
        (fn [v]
          # pole to R - x < 0 -> 0
          (if (< v 0)
            (if log? (- math/inf) 0.0)
            (let [lq (- (/ v scale))]
              (if log? (- lq (math/log scale)) (/ (math/exp lq) scale)))))
        x))))

(defn pexp
  "Exponential CDF: -expm1(-q/rate) (no cancellation); q <= 0
  -> 0, inf -> 1. Log with {:log true}. Opts: {:log false
  :impl}."
  [q rate &opt opts]
  (default opts {})
  # R pexp: -expm1(-q/scale) lower tail; q <= 0 -> 0, inf -> 1.
  # COMPOSITION: the vmap enum has no expm1 (core math/expm1 is libm)
  (def log? (get opts :log))
  (when (<= rate 0) (errorf "stats: pexp rate > 0 required, got %q" rate))
  (comp-only opts q)
  (def scale (/ 1 rate))
  (vapply
    (fn [v]
      (cond
        (<= v 0) (if log? (- math/inf) 0)
        (= v math/inf) (if log? 0 1)
        (do
          (def p (- (math/expm1 (- (/ v scale)))))
          (if log? (math/log p) p))))
    q))

(defn qexp
  "Exponential quantile: -scale*log(p) on the raw p (never a
  materialized 1-p); boundaries p=0 -> 0, p=1 -> +inf.
  Opts: {:lower-tail true :impl}."
  [p rate &opt opts]
  (default opts {})
  # R qexp: -scale * log(p); boundaries p=0 -> 0, p=1 -> +inf
  (comp-only opts p)
  (when (<= rate 0) (errorf "stats: qexp rate > 0 required, got %q" rate))
  (def scale (/ 1 rate))
  # R's -R_DT_Clog(p) - the upper lane
  # is -log(p) on the RAW p, never 1-p
(def lower-tail (q-lower-tail "qexp" opts))
  (vapply
    (fn [v]
      (if (or (not= v v) (< v 0) (> v 1)) math/nan (do
      (cond
        (= v 0) (if lower-tail 0 math/inf)
        (= v 1) (if lower-tail math/inf 0)
        (if lower-tail
          (- (/ (math/log1p (- v)) rate))
          (- (/ (math/log v) rate)))))))
    p))

(defn/typed dweibull {:args [:any :any :any :any] :ret :any}
  [x shape scale opts]
  "Weibull density (shape/scale; log with {:log true}).
  Opts: {:log false :impl}."
  # R dweibull (shape + SCALE interface): shape * tmp1 * exp(-tmp2) /
  # scale with tmp1 = pow(x/scale, shape-1), tmp2 = tmp1*(x/scale);
  # the x==0 && shape<1 edge -> +inf (R verbatim). COMPOSITION: the
  # vmap enum has no pow - math/pow per element (libm) both impls
  (def log? (get opts :log))
  (when (or (<= shape 0) (<= scale 0))
    (errorf "stats: dweibull needs shape > 0 and scale > 0, got %q/%q" shape scale))
  (comp-only opts x)
  (vapply
    (fn [v]
      (cond
        (< v 0) 0
        (= v math/inf) (if log? (- math/inf) 0)
        (and (= v 0) (< shape 1)) math/inf
        (do
          (def t1 (math/pow (/ v scale) (- shape 1)))
          (def t2 (* t1 (/ v scale)))
          (if log?
            (- (math/log (/ (* shape t1) scale)) t2)
            (/ (* shape t1 (math/exp (- t2))) scale)))))
    x))

(defn/typed pweibull {:args [:any :any :any :any] :ret :any}
  [x shape scale opts]
  "Weibull CDF: 1 - exp(-(x/scale)^shape); log with {:log
  true}. Opts: {:log false :impl}."
  # R pweibull: -expm1(-pow(x/scale, shape)); x <= 0 -> 0, inf -> 1.
  # COMPOSITION (no expm1 in the vmap enum)
  (def log? (get opts :log))
  (when (or (<= shape 0) (<= scale 0))
    (errorf "stats: pweibull needs shape > 0 and scale > 0, got %q/%q" shape scale))
  (comp-only opts x)
  (vapply
    (fn [v]
      (cond
        (<= v 0) (if log? (- math/inf) 0)
        (= v math/inf) (if log? 0 1)
        (do
          (def p (- (math/expm1 (- (math/pow (/ v scale) shape)))))
          (if log? (math/log p) p))))
    x))

  (defn/typed qweibull {:args [:any :any :any :any] :ret :any}
  [p shape scale opts]
  "Weibull quantile: scale * (-log p)^(1/shape) on the raw p
  (never a materialized 1-p); p=0 -> 0, p=1 -> +inf. Opts:
  {:lower-tail true :impl}."
  # R qweibull: scale * pow(-log(p), 1/shape)
  (comp-only opts p)
  (when (or (<= shape 0) (<= scale 0))
    (errorf "stats: qweibull needs shape > 0 and scale > 0, got %q/%q" shape scale))
  # R's -R_DT_Clog(p) - the upper lane
  # is -log(p) on the RAW p, never 1-p
(def lower-tail (q-lower-tail "qweibull" opts))
  (vapply
    (fn [v]
      (if (or (not= v v) (< v 0) (> v 1)) math/nan (do
      (cond
        (= v 0) (if lower-tail 0 math/inf)
        (= v 1) (if lower-tail math/inf 0)
        (* scale (math/pow (if lower-tail
                          (- (math/log1p (- v)))
                          (- (math/log v)))
                          (/ 1 shape)))))))
    p))

(defn/typed dlnorm {:args [:any :any :any :any] :ret :any}
  [x meanlog sdlog opts]
  "Log-normal density: delegates to dnorm on log x; x <= 0 -> 0
  (log: -inf). Opts: {:log false :impl}."
  [x meanlog sdlog opts]
  # R dlnorm: M_1_SQRT_2PI * exp(-0.5*y^2)/(x*sdlog), y = (log x -
  # meanlog)/sdlog - delegates to the LANDED dnorm on log x (the
  # composition the plan names); x <= 0 -> 0 (R returns 0;
  # an earlier throw here contradicted the comment).
  # sdlog < 0 -> error. Non-positive elements fall to elementwise.
  (def log? (get opts :log))
  (when (< sdlog 0) (errorf "stats: dlnorm sdlog >= 0 required, got %q" sdlog))
  (when (< sdlog 0) (errorf "stats: dlnorm sdlog >= 0 required, got %q" sdlog))
  # the handle detector - packed handles
  # carry the :n/:bytes marker; any OTHER struct is caller misuse,
  # actionable (never a staging crash)
  # The FULL marker (packed-handle?) - the
  # old :n/:bytes disjunction let partial markers through to a
  # staging crash
  (when (and (or (struct? x) (table? x)) (not (packed-handle? x)))
    (error "stats: dlnorm expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (if (and (= :num (num/resolve-impl opts)) (or (indexed? x) (struct? x)))
    # handles are INTENDED kernel inputs
    # (the ppois/dpois gate shape); arrays with non-positive or
    # non-number elements fall elementwise
    (let [pure-fallback (if (indexed? x)
                         (some |(or (<= $ 0) (not (number? $))) x)
                         false)]
      (if pure-fallback
        (vapply
          (fn [v]
            (if (<= v 0)
              (if log? (- math/inf) 0.0)
              (if log?
                (- (dnorm (math/log v) meanlog sdlog {:log true}) (math/log v))
                (/ (dnorm (math/log v) meanlog sdlog opts) v))))
          x)
        (let [k (k*)
              xf (num/stage x k)
              lf ((k 'vmap1) :log xf)
              dn (num/stage (dnorm lf meanlog sdlog opts) k)]
          (if log?
            ((k 'unpack-f64) ((k 'vmap2) :sub dn ((k 'vmap1) :log xf)))
            ((k 'unpack-f64) ((k 'vmap2) :div dn xf))))))
    (do
      (comp-only opts x)
      (vapply
        (fn [v]
          # x <= 0 -> 0 (log: -inf), R
          (if (<= v 0)
            (if log? (- math/inf) 0.0)
            (if log?
              (- (dnorm (math/log v) meanlog sdlog {:log true}) (math/log v))
              (/ (dnorm (math/log v) meanlog sdlog opts) v))))
        x))))

(defn/typed plnorm {:args [:any :any :any :any] :ret :any}
  [x meanlog sdlog opts]
  "Log-normal CDF: pnorm(log(x), meanlog, sdlog); x <= 0 -> 0.
  Opts: {:log false :impl}."
  [x meanlog sdlog opts]
  # R plnorm: pnorm(log(x), meanlog, sdlog); x <= 0 -> 0. Delegates
  # to the landed pnorm (kernel path inherited)
  (when (< sdlog 0) (errorf "stats: plnorm sdlog >= 0 required, got %q" sdlog))
  # the handle detector (dlnorm's shape)
  # The FULL marker (packed-handle?)
  (when (and (or (struct? x) (table? x)) (not (packed-handle? x)))
    (error "stats: plnorm expected numbers, an array, or a packed handle ({:n :bytes :ptr} from pack-f64) - a plain struct or table, or a partial marker, is caller misuse"))
  (if (and (= :num (num/resolve-impl opts)) (not (number? x)))
    (let [pure-fallback (some |(<= $ 0) x)]
      (if pure-fallback
        (vapply
          (fn [v] (if (<= v 0) 0.0 (pnorm (math/log v) meanlog sdlog opts)))
          x)
        (let [k (k*)
              lf ((k 'vmap1) :log (num/stage x k))]
          (pnorm lf meanlog sdlog opts))))
    (do
      (comp-only opts x)
      (vapply
        (fn [v]
          (if (<= v 0) 0 (pnorm (math/log v) meanlog sdlog opts)))
        x))))

(defn/typed qlnorm {:args [:any :any :any :any] :ret :any}
  [p meanlog sdlog opts]
  "Log-normal quantile: exp(qnorm(p, ...)) - the endpoints ride
  the flag-aware qnorm's -inf/inf. Opts: {:lower-tail true
  :impl}."
  [p meanlog sdlog opts]
  # R qlnorm: exp(qnorm(p, meanlog, sdlog)) - the endpoints ride the
  # flag-aware qnorm's -inf/inf (the branches were
  # once lower-lane-only, the upper lane returned 0/inf inverted)
  (when (< sdlog 0) (errorf "stats: qlnorm sdlog >= 0 required, got %q" sdlog))
  (comp-only opts p)
(def lower-tail (q-lower-tail "qlnorm" opts))
  (vapply
    (fn [v]
      (if (or (not= v v) (< v 0) (> v 1)) math/nan (do
      (cond
        (= v 0) (if lower-tail 0 math/inf)
        (= v 1) (if lower-tail math/inf 0)
        (math/exp (qnorm v meanlog sdlog opts))))))
    p))

(defn dnbinom
  "Negative-binomial PMF (size + prob parameterisation; log with
  {:log true}): composed over dbinom's conventions, the lgamma
  direct form for non-integer size (R allows it); x=0 edge
  R-verbatim. Opts: {:log false :impl}."
  [x size prob &opt opts]
  (default opts {})
  # R dnbinom (size + prob): size/(size+x) * C(size+x, size) *
  # prob^size * (1-prob)^x - composed over the LANDED dbinom (R: p =
  # size/(size+x) times dbinom_raw(size, size+x, prob, 1-prob));
  # the x==0 edge: size==0 -> 1, else prob^size (R verbatim);
  # integer prepass like dbinom (the documented VALIDATION exception)
  (def log? (get opts :log))
  (when (or (<= prob 0) (> prob 1) (< size 0))
    (errorf "stats: dnbinom needs size >= 0 and 0 < prob <= 1, got %q/%q" size prob))
  (comp-only opts x)
  (vapply
    (fn [v]
      # the relative tolerance decode
      # (r-int-tol); beyond tolerance = NaN (OUR convention)
      (def kk (r-int-tol v))
      (cond
        # poles before the tolerance decode; beyond-tolerance -> 0
        (or (< v 0) (= v math/inf)) (if log? (- math/inf) 0.0)
        (nil? kk) (if log? (- math/inf) 0.0)
        (= kk 0) (if (= size 0)
                  (if log? 0 1)
                  (do
                    (def d (math/pow prob size))
                    (if log? (math/log d) d)))
        (= size 0) (if log? (- math/inf) 0)
        (do
          # direct lgamma form (R's interface allows REAL size; the
          # dbinom composition would trip its integer prepass): ld =
          # lgamma(s+x+1) - lgamma(s+1) - lgamma(x+1) + s*log(p) +
          # x*log1p(-p) + log(s/(s+x))  (R: dbinom_raw * s/(s+x))
          (def ld (+ (lgamma (+ size kk 1))
                     (- (lgamma (+ size 1)))
                     (- (lgamma (+ kk 1)))
                     (* size (math/log prob))
                     (* kk (math/log1p (- prob)))
                     (math/log (/ size (+ size kk)))))
          (if log? ld (math/exp ld)))))
    x))

(defn pnbinom
  "Negative-binomial CDF: pbeta(prob, size, floor(x + 1e-7) + 1)
  (the pbinom pattern); size=0 -> the point mass at 0; poles
  follow R. Opts: {:log false :impl}."
  [x size prob &opt opts]
  (default opts {})
  # R pnbinom: pbeta(prob, size, floor(x + 1e-7) + 1) - the beta-inc
  # composition (the pbinom pattern); size==0 -> point mass at 0;
  # x < 0 -> 0, +inf -> 1
  (when (or (<= prob 0) (> prob 1) (< size 0))
    (errorf "stats: pnbinom needs size >= 0 and 0 < prob <= 1, got %q/%q" size prob))
  (comp-only opts x)
  (vapply
    (fn [v]
      (cond
        (= size 0) (if (>= v 0) 1 0)
        (< v 0) 0
        (= v math/inf) 1
        (pbeta prob size (+ (math/floor (+ v 1e-7)) 1) opts)))
    x))

(defn/typed qnbinom {:args [:any :any :any :any] :ret :any}
  [p size prob opts]
  "Negative-binomial quantile: R's q_DISCRETE_BODY search
  (qDiscrete_search.h) - the 8eps fuzz, the Cornish-Fisher seed,
  the upper lane walking I_{1-prob}(k+1, size) directly. Returns
  NaN elementwise for p outside [0,1]. Opts: {:lower-tail true
  :impl}."
  [p size prob opts]
  # R qnbinom via the qbinom PRECEDENT: normal-approx start (mu =
  # size*(1-prob)/prob, sigma = sqrt(size*(1-prob))/prob), walk out
  # on pnbinom - the smallest k with P(K<=k) >= p. prob=1 or size=0
  # -> 0; p boundaries 0 / +inf
  (comp-only opts p)
  (when (or (<= prob 0) (> prob 1) (< size 0))
    (errorf "stats: qnbinom needs size >= 0 and 0 < prob <= 1, got %q/%q" size prob))
  # R's q_DISCRETE_BODY
  # search (qDiscrete_search.h) - the 8eps fuzz on the RAW p (the
  # 32eps = _xf_*8eps guard), the Cornish-Fisher seed on the
  # flag-aware qnorm, and do_search on the FLAG-THREADED CDF:
  # lower walks P(K<=k) = pnbinom; upper walks P(K>k) =
  # I_{1-prob}(k+1, size) DIRECTLY (the old materialized target
  # rounded to 1.0 at subnormal p and the walk ran to inf,
  # qnbinom(1e-300,10,0.3,up) vs R 2063; and the fuzzed-window
  # off-by-ones, 166 vs 167)
(def lower-tail (q-lower-tail "qnbinom" opts))
  (vapply
    (fn [v]
      (if (or (not= v v) (< v 0) (> v 1)) math/nan (do
      (cond
        (or (= prob 1) (= size 0) (= v 0)) (if lower-tail 0 math/inf)
        (= v 1) (if lower-tail math/inf 0)
        (do
          (def e 1.7763568394002505e-15)
          (def pt (if lower-tail (* v (- 1 e))
                    (if (> (- 1 v) (* 4 e)) (* v (+ 1 e)) v)))
          (defn pdist [ky]
            (if lower-tail
              (pnbinom ky size prob {})
              (beta-inc (+ ky 1) size (- 1 prob))))
          (def mu (/ (* size (- 1 prob)) prob))
          (def sigma (/ (math/sqrt (* size (- 1 prob))) prob))
          (var k0 (math/round (+ mu (* (qnorm v 0 1
                                         (if lower-tail {} {:lower-tail false}))
                                       sigma))))
          (set k0 (max 0 k0))
          (def zc (pdist k0))
          (if (if lower-tail (>= zc pt) (< zc pt))
            (while (and (> k0 0)
                        (if lower-tail
                          (>= (pdist (- k0 1)) pt)
                          (< (pdist (- k0 1)) pt)))
              (-- k0))
            (while (if lower-tail
                     (< (pdist k0) pt)
                     (>= (pdist k0) pt))
              (++ k0)))
          k0)))))
    p))

#==------------------------------------------------------------------------==#
# The studentized range (2026-10-05): R's
# ptukey.c + qtukey.c VERBATIM (Copenhaver-Holland 1988; the
# AS 70 qinv), decoded from the r-source mirror. PURE only - the
# C kernel + the vector :num lane land after this
# verifies (the recorded sequencing).
#==------------------------------------------------------------------------==#

(def tukey-xleg @[0.981560634246719250690549090149
                  0.904117256370474856678465866119
                  0.769902674194304687036893833213
                  0.587317954286617447296702418941
                  0.367831498998180193752691536644
                  0.125233408511468915472441369464])
(def tukey-aleg @[0.047175336386511827194615961485
                  0.106939325995318430960254718194
                  0.160078328543346226334652529543
                  0.203167426723065921749064455810
                  0.233492536538354808760849898925
                  0.249147045813402785000562436043])
(def tukey-xlegq @[0.989400934991649932596154173450
                  0.944575023073232576077988415535
                  0.865631202387831743880467897712
                  0.755404408355003033895101194847
                  0.617876244402643748446671764049
                  0.458016777657227386342419442984
                  0.281603550779258913230460501460
                  0.0950125098376374401853193354250])
(def tukey-alegq @[0.0271524594117540948517805724560
                   0.0622535239386478928628438369944
                   0.0951585116824927848099251076022
                   0.124628971255533872052476282192
                   0.149595988816576732081501730547
                   0.169156519395002538189312079030
                   0.182603415044923588866763667969
                   0.189450610455068496285396723208])

(defn/typed tukey-wprob {:args [:number :number :number] :ret :number}
  [w rr cc]
  "R's wprob verbatim: Hartley's range integral - 12-point
  Legendre over [w/2, 8] in 2-3 intervals. Internal."
  # R's wprob verbatim: Hartley's range integral - 12-point
  # Legendre over [w/2, 8] in 2-3 intervals; the C1/C2/C3 gates
  (def qsqz (* w 0.5))
  (if (>= qsqz 8)
    1
    (let [phi (fn [x] (* 0.5 (erfc (/ (- x) (math/sqrt 2)))))
          pr-w0 (- (* 2 (phi qsqz)) 1)
          pr-w (if (>= pr-w0 (math/exp (/ -50 cc)))
                 (math/pow pr-w0 cc)
                 0)
          wincr (if (> w 3) 2 3)
          cc1 (- cc 1)
          binc (/ (- 8 qsqz) wincr)]
      (var blb qsqz)
      (var bub (+ blb binc))
      (var einsum 0)
      (for wi 1 (+ wincr 1)
        (var elsum 0)
        (def a (* 0.5 (+ bub blb)))
        (def b (* 0.5 (- bub blb)))
        (for jj 1 13
          (def j (if (>= 6 jj) jj (- 13 jj)))
          (def xx (if (>= 6 jj) (- (tukey-xleg (- j 1))) (tukey-xleg (- j 1))))
          (def ac (+ a (* b xx)))
          (def qexpo (* ac ac))
          (when (> qexpo 60) (break))
          (def pplus (* 2 (phi ac)))
          (def pminus (* 2 (phi (- ac w))))
          (def rinsum (- (* pplus 0.5) (* pminus 0.5)))
          (when (>= rinsum (math/exp (/ -30 cc1)))
            (+= elsum (* (* (tukey-aleg (- j 1))
                            (math/exp (* -0.5 qexpo)))
                         (math/pow rinsum cc1)))))
        (*= elsum (* (* (* 2 b) cc) 0.398942280401432677939946059934))
        (+= einsum elsum)
        (set blb bub)
        (+= bub binc))
      (def pr-w2 (+ pr-w einsum))
      (if (<= pr-w2 (math/exp (/ -30 rr)))
        0
        (let [pw (math/pow pr-w2 rr)]
          (if (>= pw 1) 1 pw))))))

(defn/typed tukey-cdf {:args [:number :number :number :number] :ret :number}
  [q rr cc df]
  "R's ptukey core verbatim: df > 25000 -> wprob directly, else
  the chi-density integral (16-point Legendre, up to 50 outer
  intervals); NaN in ANY arg -> NaN out (the 4-arg ISNAN gate).
  Internal."
  # R's ptukey core verbatim: df > 25000 -> wprob directly; else
  # the chi-density integral (16-point Legendre, up to 50 outer
  # intervals). The break rule DISCARDS the final sub-eps
  # interval (R's ans += AFTER the break - mirrored). NaN in ANY arg
  # (q/rr/cc/df) -> a NaN out (R's ISNAN gate, the full 4-arg form).
  # R's ptukey.c tests
  # ISNAN(q) || ISNAN(rr) || ISNAN(cc) || ISNAN(df) FIRST. TWO branches
  # follow: q-NaN -> q (preserve the caller's NaN payload/sign);
  # param-NaN (rr/cc/df) -> math/nan (a canonical NaN, R's
  # ML_WARN_return_NAN = return R_NaN - NOT q, which would be a FINITE
  # value for finite q, violating "NaN in, NaN out"). The kernel's
  # jn_tukey_cdf mirrors.
  (cond
    (not= q q) q
    (or (not= rr rr) (not= cc cc) (not= df df)) math/nan
    (<= q 0) 0
    (not (and (< q math/inf) (> q (- math/inf)))) 1
    (> df 25000) (tukey-wprob q rr cc)
    (let [f2 (* df 0.5)
          f2lf (+ (- (- (* f2 (math/log df)) (* df 0.693147180559945309417232121458))
                    (lgamma f2))
                  (math/log (cond
                              (<= df 100) 1
                              (<= df 800) 0.5
                              (<= df 5000) 0.25
                              0.125)))
          f21 (- f2 1)
          ff4 (* df 0.25)
          ulen (cond
                (<= df 100) 1
                (<= df 800) 0.5
                (<= df 5000) 0.25
                0.125)]
      (var ans 0)
      (var otsum 0)
      (for i 1 51
        (set otsum 0)
        (def twa1 (* (- (* 2 i) 1) ulen))
        (for jj 1 17
          (def j (if (>= 8 jj) (- jj 1) (- jj 9)))
          (def xu (* (tukey-xlegq j) ulen))
          (def t1 (if (>= 8 jj)
                  (+ f2lf (* f21 (math/log (- twa1 xu)))
                     (* (- xu twa1) ff4))
                  (+ f2lf (* f21 (math/log (+ twa1 xu)))
                     (- (* (+ xu twa1) ff4)))))
          (when (>= t1 -30)
            (def qsqz (* q (math/sqrt (* (if (>= 8 jj) (- twa1 xu) (+ twa1 xu)) 0.5))))
            (+= otsum (* (tukey-wprob qsqz rr cc)
                         (tukey-alegq j)
                         (math/exp t1)))))
        (when (and (>= (* i ulen) 1) (<= otsum 1e-14)) (break))
        (+= ans otsum))
      (min ans 1))))

(defn/typed ptukey {:args [:any :any :any :any] :ret :any}
  [q nmeans df opts]
  "Studentized-range CDF, R ptukey (Copenhaver-Holland 1988,
  ptukey.c verbatim). VALIDATION ours (R NaNs with a warning):
  df >= 2, nmeans >= 2, nranges >= 1; the NaN pole q -> the
  caller's NaN. Kernel-backed on vectors (ptukey-v, parity bar
  1e-9 rel); Opts: {:nranges 1 :log false :impl}."
  [q nmeans df opts]
  # R ptukey: the studentized-range CDF.
  # VALIDATION ours (R NaNs with a warning): df >= 2, nmeans >= 2,
  # nranges >= 1 - CALLER-side for both lanes (the ABI row)
  (def log? (get opts :log))
  (def rr (get opts :nranges 1))
  (when (or (< df 2) (< nmeans 2) (< rr 1))
    (errorf "stats: ptukey needs nmeans >= 2, df >= 2, nranges >= 1, got %q/%q/%q"
            nmeans df rr))
  # the vector :num lane through the ptukey-v kernel
  # - arrays stage, packed handles chain (RESIDENT),
  # outputs unpacked. The NaN pole AGREES across lanes
  # (janet-num 1.2.0+, the 4-arg ISNAN gate returning the caller's
  # NaN): NaN q -> NaN in the kernel and the pure - R's own value
  (if (and (= :num (num/resolve-impl opts)) (not (number? q)))
    (let [k (num/kernels)
          out ((k 'unpack-f64)
               ((k 'ptukey-v) (num/stage q k) nmeans df rr))]
      (if log? (map math/log out) out))
    (do
      (num/check-impl-landed opts q)
      (vapply
        (fn [v]
          (def val (tukey-cdf v rr nmeans df))
          (if log? (math/log val) val))
        q))))

(defn/typed qtukey-qinv {:args [:number :number :number] :ret :number}
  [p c v]
  "R's qinv (AS 70) verbatim - the secant's initial estimate.
  Internal."
  # R's qinv (AS 70) verbatim - the secant's initial estimate.
  # The rational is Horner form, highest degree first (R's nested
  # yi*p4 + p3 ... form)
  (def ps (- 0.5 (* 0.5 p)))
  (def yi (math/sqrt (math/log (/ 1 (* ps ps)))))
  (defn horner [coefs x]
    (var h (first coefs))
    (each cc (slice coefs 1) (set h (+ (* h x) cc)))
    h)
  (var t (+ yi
           (/ (horner [0.322232421088 -1 -0.342242088547
                       -0.204231210125 -4.53642210148e-05] yi)
              (horner [0.0993484626060 0.588581570495 0.531103462366
                       0.103537752850 0.0038560700634] yi))))
  (when (< v 120) (+= t (/ (+ (* t t t) t) v 4)))
  (var qv (- 0.8832 (* 0.2368 t)))
  (when (< v 120) (set qv (+ qv (- (/ 1.214 v)) (/ (* 1.208 t) v))))
  (* t (+ (* qv (math/log (- c 1))) 1.4142)))

(defn qtukey
  "Studentized-range quantile, R qtukey: the secant method (eps
  1e-4, maxiter 50) over ptukey with the AS 70 start; the p
  argument is a probability in [0,1] (R convention), out-of-domain
  -> NaN (R's NaN + warning). Near p=1 the documented R
  solver-floor divergences live in api-contract. Opts: {:nranges
  1 :lower-tail true :impl}."
  [p &opt nmeans df opts]
  # R qtukey: the secant method (eps 1e-4,
  # maxiter 50) over ptukey with the AS 70 start. Ours: p outside
  # (0,1) -> VALIDATION (R NaNs); the degenerate-secant guard
  # breaks (R's dead valx1 = -p assignment not mirrored)
  (default opts {})
  (def rr (get opts :nranges 1))
  # the flip - R qtukey.c:140 p = R_DT_qIv(p)
  # (the 0.5-p+0.5 complement), then the same secant on the lower lane
(def lower-tail (q-lower-tail "qtukey" opts))
  (when (not (number? p))
    (errorf "stats: qtukey p must be a number, got %q" p))
  (when (or (< df 2) (< nmeans 2) (< rr 1))
    (errorf "stats: qtukey needs nmeans >= 2, df >= 2, nranges >= 1, got %q/%q/%q"
            nmeans df rr))
  # out-of-domain or NaN p -> NaN
  # (R's NaN + warning; the old THROW is withdrawn - the q-pole
  # convention, family-wide now)
  (cond
    (or (not= p p) (< p 0) (> p 1)) math/nan
    (= p 0) (if lower-tail 0 math/inf)
    (= p 1) (if lower-tail math/inf 0)
    (do
      (def pp (if lower-tail p (+ (- 0.5 p) 0.5)))
      (def x0 (qtukey-qinv pp nmeans df))
      (var valx0 (- (tukey-cdf x0 rr nmeans df) pp))
      (var x1 (if (> valx0 0) (max 0 (- x0 1)) (+ x0 1)))
      (var valx1 (- (tukey-cdf x1 rr nmeans df) pp))
      (var xx0 x0)
      (var ans 0)
      (for iter 0 50
        (if (= valx1 valx0) (break))
        (set ans (- x1 (/ (* valx1 (- x1 xx0)) (- valx1 valx0))))
        (set valx0 valx1)
        (set xx0 x1)
        (when (< ans 0) (set ans 0))
        (set valx1 (- (tukey-cdf ans rr nmeans df) pp))
        (set x1 ans)
        (when (< (math/abs (- x1 xx0)) 1e-4) (break)))
      ans)))

#==------------------------------------------------------------------------==#
# draw-from - Box-Muller for normals (O(1) per draw; qnorm inversion is
# ~100 bisection iterations and too slow for sampler inner loops);
# inversion via the q family for gamma/beta (no closed form).
# math/rng is deterministic per seed string, so draws
# replay identically. RNG paths are never compared across languages
# (charter) - conformance is moments + determinism.
#==------------------------------------------------------------------------==#

(defn/typed rnorm-pair {:args [:any] :ret :tuple} [rng]
  "Two independent standard normals from two uniforms
  (Box-Muller; O(1) per draw - qnorm inversion is too slow for
  sampler inner loops). Internal."
  # two independent standard normals from two uniforms (Box-Muller)
  (def u1 (math/rng-uniform rng))
  (def u2 (math/rng-uniform rng))
  (def r (math/sqrt (* -2 (math/log (max u1 1e-300)))))
  (def a (* 2 math/pi u2))
  [(* r (math/cos a)) (* r (math/sin a))])

(defn/typed rnorm {:args [:any :number :number :number] :ret :array} [rng n mean sd]
  "n normal draws (mean, sd) from a math/rng. Deterministic per
  rng state; conformance is moments + determinism (RNG paths are
  never cross-compared)."
  (def out (array/new n))
  (var i 0)
  (while (< i n)
    (def [z1 z2] (rnorm-pair rng))
    (array/push out (+ mean (* sd z1)))
    (++ i)
    (when (< i n)
      (array/push out (+ mean (* sd z2)))
      (++ i)))
  out)

(defn/typed rgamma-one {:args [:any :number :number] :ret :number} [rng shape rate]
  "One gamma draw (Marsaglia-Tsang direct; the shape < 1 boost).
  Deterministic per rng state. Internal."
  # Marsaglia-Tsang direct gamma draw (shape >= 1 branch + the shape < 1
  # boost) - O(1), unlike qgamma inversion (~100 bisections). Deterministic
  # per rng state; conformance is moments (charter: never compare RNG paths).
  (when (or (<= shape 0) (<= rate 0)) (error "stats: rgamma needs shape, rate > 0"))
  (defn mt [sh]
    (def d (- sh 0.3333333333333333))
    (def c (/ 1 (math/sqrt (* 9 d))))
    (var v nil)
    (while (nil? v)
      (def [z1 z2] (rnorm-pair rng))
      (def x z1)
      (def vv (* (+ 1 (* c x)) (+ 1 (* c x)) (+ 1 (* c x))))
      (when (> vv 0)
        (def u (math/rng-uniform rng))
        (when (or (< u (- 1 (* 0.0331 x x x x)))
                  (< (math/log (max u 1e-300))
                     (+ (* 0.5 x x) (* d (+ 1 (- vv) (math/log vv))))))
          (set v (* d vv)))))
    v)
  (def g
    (if (>= shape 1)
      (mt shape)
      # shape < 1: draw Gamma(shape+1) and multiply by U^{1/shape}
      (* (mt (+ shape 1)) (math/pow (max (math/rng-uniform rng) 1e-300) (/ 1 shape)))))
  (/ g rate))

(defn/typed rgamma {:args [:any :number :number :number] :ret :array} [rng n shape rate]
  "n gamma draws (RATE parameterisation) from a math/rng.
  Deterministic per rng state; conformance is moments."
  (seq [_ :range [0 n]] (rgamma-one rng shape rate)))

(defn/typed rbeta {:args [:any :number :number :number] :ret :array} [rng n a b]
  "n beta draws from a math/rng (gamma pair, normalised).
  Deterministic per rng state; conformance is moments."
  (seq [_ :range [0 n]] (qbeta (math/rng-uniform rng) a b {})))
