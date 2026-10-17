# stats/regsurv.janet - regression + survival.
# OLS full inference (R lm), GLM IRLS (binomial-logit +
# poisson-log, Wald CIs), Cox PH (Breslow + Efron ties), CLOGIT
# (stratified partial likelihood; survival::clogit - the standard
# matched-sets form, exact for 1-case strata, Breslow-style for
# multi-case).
# lm/glm are kernel-backed: X'X/X'WX through jn_matmul (the GEMM -
# the measured top evidence), X'y/X'Wz + fitted through jn_matvec,
# residual/rss/tss through the vector kernels; the p x p
# solve/inverse ride the linalg module's own num path. coxph/clogit
# stay pure.
(import jennystats/typed :prefix "")
(import ./utils :prefix "")
(import ./dist :prefix "")
(import ./descriptive :prefix "")
(import ./linalg :prefix "")
(import ./num)

#==------------------------------------------------------------------------==#
# small dense-matrix core (p is small in the stats regime)
#==------------------------------------------------------------------------==#

(defn/typed m-rows {:args [:any] :ret :number} [m] "Row-major matrix row count." (length m))
(defn/typed m-cols {:args [:any] :ret :number} [m] "Row-major matrix column count." (length (m 0)))

(defn/typed solve-sp {:args [:any :any] :ret :array} [A b]
  "Solve A x = b; thin delegation to the linalg module (single
  owner). Internal."
  (solve A b {}))

(defn/typed inverse-sp {:args [:any] :ret :array} [A]
  "Matrix inverse; thin delegation to the linalg module (single
  owner). Internal."
  (inverse A {}))

(defn/typed zeros-m {:args [:number :number] :ret :array} [rows cols]
  "Zero matrix with distinct row arrays - array/new-filled shares
  the fill value between rows."
  (def out @[])
  (for i 0 rows (array/push out (array/new-filled cols 0)))
  out)

(defn/typed xtwx {:args [:any :any] :ret :array} [X w]
  "X'WX for per-row weights w. Internal."
  (def n (m-rows X))
  (def p (m-cols X))
  (def out (zeros-m p p))
  (for i 0 p
    (for j 0 p
      (var s 0)
      (for r 0 n (+= s (* (get w r) (get (get X r) i) (get (get X r) j))))
      (put (get out i) j s)))
  out)

(defn/typed xtwv {:args [:any :any :any] :ret :array} [X w v]
  "X'W v as a p-vector. Internal."
  (def n (m-rows X))
  (def p (m-cols X))
  (def out (array/new-filled p 0))
  (for i 0 p
    (var s 0)
    (for r 0 n (+= s (* (get w r) (get (get X r) i) (get v r))))
    (put out i s))
  out)

(defn/typed mmul {:args [:any :any] :ret :array} [A B]
  "Row-major matrix product. Internal."
  (def n (m-rows A))
  (def m (m-cols B))
  (def k (m-cols A))
  # zeros-m, NOT (array/new-filled n (array/new-filled m 0)) - the latter
  # shares one row array across all n rows and silently corrupts the result
  (def out (zeros-m n m))
  (for i 0 n
    (for j 0 m
      (var s 0)
      (for t 0 k (+= s (* (get (get A i) t) (get (get B t) j))))
      (put (get out i) j s)))
  out)

(defn/typed matvec {:args [:any :any] :ret :array} [A v]
  "A (n x k) times v (k-vector) -> n-vector. Internal."
  (def n (m-rows A))
  (def k (length v))
  (def out (array/new-filled n 0))
  (for i 0 n
    (var s 0)
    (for j 0 k (+= s (* (get (get A i) j) (get v j))))
    (put out i s))
  out)

(defn/typed eta-of {:args [:any :any] :ret :number} [X-row beta]
  "The linear predictor X-row . beta. Internal."
  (var eta 0)
  (for j 0 (length beta) (+= eta (* (get beta j) (get X-row j))))
  eta)

#==------------------------------------------------------------------------==#
# num-path staging - the glue collapse: the cross-product
# accumulation is a GEMM; the p x p solve/inverse stay on the linalg
# module's own num path (p is small in the stats regime)
#==------------------------------------------------------------------------==#

(defn/typed xt-flat {:args [:any :number :number] :ret :array} [X n p]
  "X' flat row-major (p x n) - the one O(np) transpose; the
  O(np^2) accumulation behind it is the GEMM. Internal."
  (def out (array/new (* p n)))
  (for j 0 p
    (for i 0 n (put out (+ (* j n) i) (get (get X i) j))))
  out)

(defn/typed tile-each {:args [:any :number] :ret :array} [v times]
  "v[i] repeated `times` consecutively - the row-broadcast of w
  for XW = w .* X (flat n*p: index i*p+j wants w[i]). Internal."
  (def out (array/new (* (length v) times)))
  (each x v (for _ 0 times (array/push out x)))
  out)

(defn/typed lm-num-core {:args [:any :any :number :number] :ret :table} [X y n p]
  "lm core via the kernels: X'X via GEMM, X'y + fitted via matvec,
  resid/rss/tss via the vector ops - pack once, chain kernel calls
  (resident-packed). Internal."
  (def k (num/kernels))
  (def nopts {:impl :num})
  (def Xf (num/stage (flat-row-major X) k))
  (def Xt (num/stage (xt-flat X n p) k))
  (def yf (num/stage y k))
  (def xtx (rebuild-rows ((k 'unpack-f64) ((k 'matmul) Xt Xf p n p)) p p))
  (def xty ((k 'unpack-f64) ((k 'matvec) Xt yf p n)))
  (def beta (solve xtx xty nopts))
  (def fitted-h ((k 'matvec) Xf (num/stage beta k) n p))
  (def resid-h ((k 'vmap2) :sub yf fitted-h))
  (def rss ((k 'vdot) resid-h resid-h))
  (def inv (inverse xtx nopts))
  (def ybar (/ ((k 'vreduce) :sum yf) n))
  (def yc ((k 'vmap1-scalar) :sub yf ybar))
  (def tss ((k 'vdot) yc yc))
  @{:beta beta :fitted ((k 'unpack-f64) fitted-h)
    :resid ((k 'unpack-f64) resid-h) :rss rss :inv inv
    :ybar ybar :tss tss})

(defn/typed lm-pure-core {:args [:any :any :number :number] :ret :table} [X y n p]
  "lm core via the pure loops (the reference lane). Internal."
  (def ones (array/new-filled n 1))
  (def xtx (xtwx X ones))
  (def xty (xtwv X ones y))
  (def beta (solve-sp xtx xty))
  (def fitted (array/new n))
  (def resid (array/new n))
  (for i 0 n
    (def f (eta-of (get X i) beta))
    (put fitted i f)
    (put resid i (- (get y i) f)))
  (var rss 0)
  (each r resid (+= rss (* r r)))
  (def inv (inverse-sp xtx))
  (var ybar 0)
  (each v y (+= ybar v))
  (set ybar (/ ybar n))
  (var tss 0)
  (each v y (+= tss (* (- v ybar) (- v ybar))))
  @{:beta beta :fitted fitted :resid resid :rss rss :inv inv
    :ybar ybar :tss tss})

#==------------------------------------------------------------------------==#
# OLS - R lm()
#==----------------------------------------------------------------==========#

(defn lm
  "Ordinary least squares, R lm(). X = n x p WITH the intercept
  column included (1s first is R's order); returns {:coefficients
  :se :t :p :ci :sigma :df :r2 :adj-r2 :f :f-p :residuals
  :fitted}. Opts: {:conf-level 0.95 :impl}; throws when n <= p."
  [X y &opt opts]
  # X: n x p WITH intercept column included (1s first is R's order).
  (default opts {})
  (def n (m-rows X))
  (def p (m-cols X))
  (when (<= n p) (error "stats: lm needs n > p"))
  (def core
    (if (= :num (num/resolve-impl opts))
      (lm-num-core X y n p)
      (do (num/check-impl-landed opts) (lm-pure-core X y n p))))
  (def beta (core :beta))
  (def fitted (core :fitted))
  (def resid (core :resid))
  (def rss (core :rss))
  (def inv (core :inv))
  (def ybar (core :ybar))
  (def tss (core :tss))
  (def dof (- n p))
  (def sigma2 (/ rss dof))
  (def ses (array/new p))
  (def tvals (array/new p))
  (def pvals (array/new p))
  (def ci-lo (array/new p))
  (def ci-hi (array/new p))
  (def conf (get opts :conf-level 0.95))
  (def tcrit (qt (- 1 (* 0.5 (- 1 conf))) dof {}))
  (for j 0 p
    (def se (math/sqrt (* sigma2 (get (get inv j) j))))
    (put ses j se)
    (def t (/ (get beta j) se))
    (put tvals j t)
    (put pvals j (* 2 (- 1 (pt (math/abs t) dof {}))))
    (put ci-lo j (- (get beta j) (* tcrit se)))
    (put ci-hi j (+ (get beta j) (* tcrit se))))
  (def r2 (- 1 (/ rss tss)))
  (def adj-r2 (- 1 (* (- 1 r2) (/ (- n 1) dof))))
  (def f (* (/ (- tss rss) (- p 1)) (/ dof rss)))
  {:coefficients beta :se ses :t tvals :p pvals
   :ci [ci-lo ci-hi] :sigma (math/sqrt sigma2) :df dof
   :r2 r2 :adj-r2 adj-r2 :f f
   :f-p (- 1 (pf f (- p 1) dof {}))
   :residuals resid :fitted fitted})

#==------------------------------------------------------------------------==#
# GLM IRLS - R glm(binomial(logit)) / glm(poisson(log))
#==----------------------------------------------------------------==========#

(defn glm
  "GLM by IRLS, R glm(binomial(logit)) / glm(poisson(log)):
  full Wald inference off the LAST IRLS weighted least squares.
  Returns {:coefficients :se :z :p :ci :family :deviance :iters
  :converged}. Convergence follows R glm.fit (deviance relative
  change < epsilon); throws on a non-finite deviance step or
  non-convergence. Opts: {:epsilon 1e-8 :max-it 25 :conf-level 0.95
  :impl}; throws when n <= p."
  [X y family &opt opts]
  (default opts {})
  (def impl (num/resolve-impl opts))
  (def n (m-rows X))
  (def p (m-cols X))
  (when (<= n p) (error "stats: glm needs n > p"))
  (def num? (= :num impl))
  (when (= :pure impl) (num/check-impl-landed opts))
  # staging: X/X' packed ONCE (the resident-packed core); w/z
  # pack per iteration (O(n) glue next to the O(np^2) GEMM)
  (def k (if num? (num/kernels)))
  (def Xf (if num? (num/stage (flat-row-major X) k)))
  (def Xt (if num? (num/stage (xt-flat X n p) k)))
  (def eps (get opts :epsilon 1e-8))
  (def max-it (get opts :max-it 25))
  (def beta (array/new-filled p 0))
  (def w (array/new n))
  (def z (array/new n))
  (def mu (array/new n))
  (def eta (array/new n))
  (var iters 0)
  (var converged false)
  # R glm.fit conventions (decoded from R 4.5 glm.fit source):
  # - mustart init: binomial (y+0.5)/2, poisson y+0.1 - NOT intercept-only;
  #   the iteration trajectory matters because summary's covariance comes
  #   from the LAST IRLS QR (weights at the second-to-last beta)
  # - deviance = sum of dev.resids (DECREASE = improvement)
  # - z-form solve returns the ABSOLUTE new beta (assign, never add)
  # - converge when |dev - devold| / (0.1 + |dev|) < epsilon (default 1e-8)
  (for i 0 n
    (case family
      :binomial (let [m (/ (+ (get y i) 0.5) 2)]
                  (put mu i m)
                  (put eta i (math/log (/ m (- 1 m)))))
      :poisson (let [m (+ (get y i) 0.1)]
                 (put mu i m)
                 (put eta i (math/log m)))
      (errorf "stats: glm family %q not supported" family)))
  (defn deviance-mu []
    (var dv 0)
    (for i 0 n
      (def m (get mu i))
      (def yv (get y i))
      (case family
        :binomial (+= dv (* 2 (+ (if (> yv 0) (* yv (math/log (/ yv m))) 0)
                                 (if (< yv 1) (* (- 1 yv) (math/log (/ (- 1 yv) (- 1 m)))) 0))))
        :poisson (+= dv (* 2 (- (if (> yv 0) (* yv (math/log (/ yv m))) 0)
                               (- yv m))))))
    dv)
  (defn fin? [v] (and (= v v) (< (math/abs v) math/inf)))
  (var dev-old nil)
  (var last-xtwx nil)
  (while (and (not converged) (< iters max-it))
    (++ iters)
    (for i 0 n
      (def e (get eta i))
      (def m (get mu i))
      (def yv (get y i))
      (case family
        :binomial (let [wi (max (* m (- 1 m)) 1e-10)]
                    (put w i wi)
                    (put z i (+ e (/ (- yv m) wi))))
        :poisson (let [wi (max m 1e-10)]
                   (put w i wi)
                   (put z i (+ e (/ (- yv m) wi))))))
    (def step
      (if num?
        (let [wt (num/stage (tile-each w p) k)
              xw ((k 'vmap2) :mul Xf wt)
              xtwx (rebuild-rows ((k 'unpack-f64) ((k 'matmul) Xt xw p n p)) p p)
              wz ((k 'vmap2) :mul (num/stage w k) (num/stage z k))
              xtwz ((k 'unpack-f64) ((k 'matvec) Xt wz p n))]
          (set last-xtwx xtwx)
          (solve xtwx xtwz {:impl :num}))
        (let [xtwx (xtwx X w)]
          (set last-xtwx xtwx)
          (solve-sp xtwx (xtwv X w z)))))
  # z-form IRLS: the solve returns the ABSOLUTE new beta, not a delta -
  # z = eta + (y-mu)/w makes X'Wz = X'WX*beta + X'(y-mu). Assign, don't add.
    (for j 0 p (put beta j (get step j)))
    (if num?
      (let [ea ((k 'unpack-f64) ((k 'matvec) Xf (num/stage beta k) n p))]
        (for i 0 n (put eta i (get ea i))))
      (for i 0 n (put eta i (eta-of (get X i) beta))))
    (for i 0 n
      (def e (get eta i))
      (case family
        :binomial (put mu i (/ 1 (+ 1 (math/exp (- e)))))
        :poisson (put mu i (math/exp e))))
    (def dev-new (deviance-mu))
    (unless (fin? dev-new)
      (error "stats: glm step produced non-finite deviance (boundary)"))
    (when (and dev-old
               (< (math/abs (- dev-old dev-new)) (* eps (+ 0.1 (math/abs dev-new)))))
      (set converged true))
    (set dev-old dev-new))
  (unless converged (errorf "stats: glm did not converge in %d iters" max-it))
  (def inv (if num? (inverse last-xtwx {:impl :num}) (inverse-sp last-xtwx)))
  (def ses (array/new p))
  (def z-vals (array/new p))
  (def pvals (array/new p))
  (def conf (get opts :conf-level 0.95))
  (def zcrit (qnorm (- 1 (* 0.5 (- 1 conf))) 0 1 {}))
  (def ci-lo (array/new p))
  (def ci-hi (array/new p))
  (for j 0 p
    (def se (math/sqrt (get (get inv j) j)))
    (put ses j se)
    (def zv (/ (get beta j) se))
    (put z-vals j zv)
    (put pvals j (* 2 (- 1 (pnorm (math/abs zv) 0 1 {}))))
    (put ci-lo j (- (get beta j) (* zcrit se)))
    (put ci-hi j (+ (get beta j) (* zcrit se))))
  {:coefficients beta :se ses :z z-vals :p pvals :ci [ci-lo ci-hi]
   :family family :deviance dev-old :iters iters :converged true})

#==------------------------------------------------------------------------==#
# Cox PH - Breslow + Efron ties, analytic score + information.
# Flat style: the four accumulation blocks share two helpers.

(defn/typed risk-sums! {:args [:any :any :any :any :any :any :any] :ret :number}
  [X beta order k-start k-end s1-acc s2-acc]
  "exp-weighted sums over rows order[k-start..k-end) into the
  accumulators; returns s0, mutates s1-acc / s2-acc (s2: p x p
  rows). Internal to coxph."
  (var s0 0)
  (def p (length s1-acc))
  (for k k-start k-end
    (def row (get X (get order k)))
    (def e (math/exp (eta-of row beta)))
    (+= s0 e)
    (for c 0 p
      (put s1-acc c (+ (get s1-acc c) (* e (get row c))))
      (for c2 0 p
        (put (get s2-acc c) c2
             (+ (get (get s2-acc c) c2)
                (* e (get row c) (get row c2)))))))
  s0)

(defn/typed add-event! {:args [:any :any :any :any :any :any :any :any :any] :ret :number}
  [X beta order k s0 s1 s2 g I]
  "One event row's score/info/loglik contribution against given
  sums; returns the ll increment, mutates g / I. Internal to coxph."
  (def row (get X (get order k)))
  (def p (length s1))
  (var ll 0)
  (+= ll (- (eta-of row beta) (math/log s0)))
  (for c 0 p
    (put g c (+ (get g c) (- (get row c) (/ (get s1 c) s0))))
    (for c2 0 p
      (put (get I c) c2
           (+ (get (get I c) c2)
              (- (/ (get (get s2 c) c2) s0)
                 (/ (* (get s1 c) (get s1 c2)) (* s0 s0)))))))
  ll)

(defn coxph
  "Cox proportional hazards, survival::coxph: Newton-Raphson on
  the partial likelihood with analytic score + information. X =
  n x p; time = survival times; event = 0/1 (1 = event) or 1/2
  (2 = event, R Surv convention). Returns {:coefficients :se :z :p
  :ties :iters :loglik :converged}. Opts: {:ties :efron|:breslow
  :epsilon 1e-10 :max-it 200 :impl}; throws on NR non-convergence.
  SPEC: the risk sets are NESTED PREFIXES of the descending-time
  order - the num path's running sums associate identically to
  fresh re-sums, so both impls agree bitwise."
  [X time event &opt opts]
  # X: n x p; time: survival time; event: 0/1 (1=event) or 1/2 (2=event,
  # R Surv convention). Newton-Raphson on the partial likelihood;
  # ties :efron (R coxph default) or :breslow. Conformance:
  # survival::coxph.
  # the risk sets are NESTED PREFIXES of the descending-time
  # order - the num path takes ONE vscan per accumulator (O(n*p^2)/pass)
  # instead of fresh re-sums per event time (O(n^2*p^2)); running sums
  # associate identically to fresh left-to-right re-sums, so the pass
  # is BITWISE equal to pure and the NR trajectory matches exactly
  (default opts {})
  (num/check-impl-landed opts)
  (def tie-method (get opts :ties :efron))
  (def n (m-rows X))
  (def p (m-cols X))
  (def dead-value (if (some |(= $ 2) event) 2 1))
  (def dead (map (fn [e] (if (= e dead-value) 1 0)) event))
  (def order (sorted (seq [i :range [0 n]] i) (fn [a b] (> (get time a) (get time b)))))
  (def eps (get opts :epsilon 1e-10))
  (def max-it (get opts :max-it 200))

  (defn pass [beta]
    # one full pass over distinct times -> [g I ll]
    (def num? (= :num (num/resolve-impl opts)))
    (def pre
      # ordered gathers + exp + p^2 prefix scans - staged once per pass
      (when num?
        (let [kk (num/kernels)
              xo (seq [i :range [0 n]] (get X (get order i)))
              xot (seq [c :range [0 p]]
                     (seq [i :range [0 n]] (get (get X (get order i)) c)))
              xof (num/stage (flat-row-major xo) kk)
              e-o ((kk 'vmap1) :exp
                   ((kk 'matvec) xof (num/stage beta kk) n p))
              colf (seq [c :range [0 p]] (num/stage (get xot c) kk))
              excf (seq [c :range [0 p]] ((kk 'vmap2) :mul e-o (get colf c)))]
          {:s0 ((kk 'unpack-f64) ((kk 'vscan) :sum e-o))
           :s1 (seq [c :range [0 p]]
                 ((kk 'unpack-f64) ((kk 'vscan) :sum (get excf c))))
           :s2 (seq [c :range [0 p]]
                 (seq [c2 :range [0 p]]
                   ((kk 'unpack-f64)
                    ((kk 'vscan) :sum ((kk 'vmap2) :mul (get excf c) (get colf c2))))))})))
    (defn prefix-vals [j]
      # risk set order[0..j): the running sums at index j-1
      (def ji (dec j))
      [((pre :s0) ji)
       (map |(get $ ji) (pre :s1))
       (map |(map |(get $ ji) $) (pre :s2))])
    (def g (array/new-filled p 0))
    (def I (zeros-m p p))
    (var ll 0)
    (var i 0)
    (while (< i n)
      (def t (get time (get order i)))
      (var j i)
      (while (and (< j n) (>= (get time (get order j)) t)) (++ j))
      (def deaths @[])
      (for k i j
        (when (= (get dead (get order k)) 1)
          (array/push deaths k)))
      (when (pos? (length deaths))
        (def d (length deaths))
        # risk set = ALL rows with time >= t = order[0..j) (descending time)
        (def [s0 s1 s2]
          (if num?
            (prefix-vals j)
            (do
              (def s1p (array/new-filled p 0))
              (def s2p (zeros-m p p))
              [(risk-sums! X beta order 0 j s1p s2p) s1p s2p])))
        (if (= tie-method :breslow)
          (each k deaths
            (+= ll (add-event! X beta order k s0 s1 s2 g I)))
          # Efron: l-th death uses s0 - (l/d)*tie-sums
          (do
            (def s1-tie (array/new-filled p 0))
            (def s2-tie (zeros-m p p))
            (def s0-tie
              (do
                (var s0t 0)
                (each k deaths
                  (def row (get X (get order k)))
                  (def e (math/exp (eta-of row beta)))
                  (+= s0t e)
                  (for c 0 p
                    (put s1-tie c (+ (get s1-tie c) (* e (get row c))))
                    (for c2 0 p
                      (put (get s2-tie c) c2
                           (+ (get (get s2-tie c) c2)
                              (* e (get row c) (get row c2)))))))
                s0t))
            (for l 0 d
              (def frac (/ l d))
              (def s0-l (- s0 (* frac s0-tie)))
              (def s1-l (array/new-filled p 0))
              (def s2-l (zeros-m p p))
              (for c 0 p
                (put s1-l c (- (get s1 c) (* frac (get s1-tie c))))
                (for c2 0 p
                  (put (get s2-l c) c2
                       (- (get (get s2 c) c2)
                          (* frac (get (get s2-tie c) c2))))))
              (def k (get deaths l))
              (+= ll (add-event! X beta order k s0-l s1-l s2-l g I))))))
      (set i j))
    [g I ll])
  (def beta (array/new-filled p 0))
  (var iters 0)
  (var converged false)
  (var ll-final 0)
  (while (and (not converged) (< iters max-it))
    (++ iters)
    (def [g I ll] (pass beta))
    (set ll-final ll)
    (def step (matvec (inverse-sp I) g))
    (var max-delta 0)
    (for j 0 p
      (put beta j (+ (get beta j) (get step j)))
      (when (> (math/abs (get step j)) max-delta)
        (set max-delta (math/abs (get step j)))))
    (when (< max-delta eps) (set converged true)))
  (unless converged
    (errorf "stats: coxph NR did not converge in %d iters" max-it))
  (def [_ I-final _] (pass beta))
  (def inv-I (inverse-sp I-final))
  (def ses (array/new p))
  (def z-vals (array/new p))
  (def pvals (array/new p))
  (for j 0 p
    (def se (math/sqrt (get (get inv-I j) j)))
    (put ses j se)
    (put z-vals j (/ (get beta j) se))
    (put pvals j (* 2 (- 1 (pnorm (math/abs (/ (get beta j) se)) 0 1 {})))))
  {:coefficients beta :se ses :z z-vals :p pvals :ties tie-method
   :iters iters :loglik ll-final :converged true})

#==------------------------------------------------------------------------==#
# CLOGIT - conditional logistic, stratified partial likelihood
# (survival::clogit encoding; exact for 1-case strata)
#==------------------------------------------------------------------------==#

(defn clogit
  "Conditional logistic regression, survival::clogit: stratified
  partial likelihood, exact for 1-case strata. X = n x p; y = 0/1
  outcomes; strata = the matched-set id per row. Returns
  {:coefficients :se :z :p :iters :loglik :converged}. Opts:
  {:epsilon 1e-10 :max-it 200 :impl}; throws on NR non-convergence.
  Composition (matched sets are small - per-stratum staging would
  cost more than the per-set loops)."
  [X y strata &opt opts]
  (default opts {})
  # matched-set strata are small
  # (the design) - per-stratum staging would cost more than the
  # per-set loops; :num runs the identical janet code
  (num/check-impl-landed opts)
  (def n (m-rows X))
  (def p (m-cols X))
  (def eps (get opts :epsilon 1e-10))
  (def max-it (get opts :max-it 200))
  # group rows by stratum
  (def strata-map @{})
  (for i 0 n
    (def s (get strata i))
    (unless (get strata-map s) (put strata-map s @[]))
    (array/push (get strata-map s) i))
  (def groups (map (fn [kv] (kv 1)) (pairs strata-map)))

  (defn pass [beta]
    (def g (array/new-filled p 0))
    (def I (zeros-m p p))
    (var ll 0)
    (each rows groups
      (var ncase 0)
      (each i rows (when (= (get y i) 1) (++ ncase)))
      (when (pos? ncase)
        # stratum risk set = all its rows (cases are the events)
        (def s1 (array/new-filled p 0))
        (def s2 (zeros-m p p))
        (var s0 0)
        (each i rows
          (def row (get X i))
          (def e (math/exp (eta-of row beta)))
          (+= s0 e)
          (for c 0 p
            (put s1 c (+ (get s1 c) (* e (get row c))))
            (for c2 0 p
              (put (get s2 c) c2
                   (+ (get (get s2 c) c2)
                      (* e (get row c) (get row c2)))))))
        (each i rows
          (when (= (get y i) 1)
            (+= ll (- (eta-of (get X i) beta) (math/log s0)))
            (for c 0 p
              (put g c (+ (get g c) (- (get (get X i) c) (/ (get s1 c) s0))))
              (for c2 0 p
                (put (get I c) c2
                     (+ (get (get I c) c2)
                        (- (/ (get (get s2 c) c2) s0)
                           (/ (* (get s1 c) (get s1 c2)) (* s0 s0)))))))))))
    [g I ll])

  (def beta (array/new-filled p 0))
  (var iters 0)
  (var converged false)
  (var ll-final 0)
  (while (and (not converged) (< iters max-it))
    (++ iters)
    (def [g I ll] (pass beta))
    (set ll-final ll)
    (def step (matvec (inverse-sp I) g))
    (var max-delta 0)
    (for j 0 p
      (put beta j (+ (get beta j) (get step j)))
      (when (> (math/abs (get step j)) max-delta)
        (set max-delta (math/abs (get step j)))))
    (when (< max-delta eps) (set converged true)))
  (unless converged
    (errorf "stats: clogit NR did not converge in %d iters" max-it))
  (def [_ I-final] (pass beta))
  (def inv-I (inverse-sp I-final))
  (def ses (array/new p))
  (def z-vals (array/new p))
  (def pvals (array/new p))
  (for j 0 p
    (def se (math/sqrt (get (get inv-I j) j)))
    (put ses j se)
    (put z-vals j (/ (get beta j) se))
    (put pvals j (* 2 (- 1 (pnorm (math/abs (/ (get beta j) se)) 0 1 {})))))
  {:coefficients beta :se ses :z z-vals :p pvals
   :iters iters :loglik ll-final :converged true})

#==------------------------------------------------------------------------==#
# Kaplan-Meier + two-group log-rank. SPEC km-logrank.qnt. THE tie
# convention: deaths BEFORE censorings - every risk set is
# #{time >= t}, INCLUDING those censored at the same time
# (confirmed: survdiff.fit's order(strat, time, -status) sort +
# survfit's n.risk, R 4.5.0 / survival 3.8.6). Dual-impl: the
# block walks stay janet (sort-bound composition); the S curve =
# vscan :prod of the factors (the coxph pattern), the log-rank
# accumulations = vreduce :sum - identical sequential order both
# paths. NO new kernels.
#==------------------------------------------------------------------------==#

(defn/typed norm-events {:args [:any] :ret :array} [statuses]
  "Statuses 0/1 (1 = event) or 1/2 (2 = event, R Surv convention)
  normalised to 0/1; throws on anything else."
  (def out (array/new-filled (length statuses) 0))
  (def mx (reduce max 0 statuses))
  (cond
    (or (= mx 1) (= mx 0))
    (for i 0 (length statuses)
      (unless (or (= (statuses i) 0) (= (statuses i) 1))
        (errorf "stats: status must be 0/1 or 1/2, got %q at index %d"
                (statuses i) i))
      (put out i (statuses i)))
    (= mx 2)
    (for i 0 (length statuses)
      (unless (or (= (statuses i) 1) (= (statuses i) 2))
        (errorf "stats: status must be 0/1 or 1/2, got %q at index %d"
                (statuses i) i))
      (put out i (- (statuses i) 1)))
    (errorf "stats: status must be 0/1 or 1/2, got %q" mx))
  out)

(defn surv-check
  "The survival guard pair: composition (sort-bound, janet arrays
  only) + VALIDATION (equal lengths, non-empty, finite times).
  Internal."
  [opts times statuses]
  # the composition guard (sort-bound) + VALIDATION
  (num/check-impl-landed opts)
  (each x [times statuses]
    (unless (indexed? x)
      (error "stats: composition path consumes janet arrays; unpack packed handles first")))
  (unless (= (length times) (length statuses))
    (error "stats: times and statuses must have equal length"))
  (when (zero? (length times))
    (error "stats: survival needs at least one observation"))
  (each t times
    (unless (and (< t math/inf) (> t (- math/inf)))
      (error "stats: survival times must be finite")))
  nil)

(defn/typed km-blocks {:args [:any :any] :ret :any} [times events]
  "The distinct-time blocks (sort-bound): [time size deaths
  censors]. Internal to the survival fns."
  (def pts (sorted (do (def acc @[])
                        (for i 0 (length times)
                          (array/push acc [(times i) (events i)]))
                        acc)
                   (fn [p q] (< (p 0) (q 0)))))
  (def blocks @[])
  (var i 0)
  (def n (length pts))
  (while (< i n)
    (def t ((pts i) 0))
    (var j i)
    (var d 0)
    (while (and (< j n) (= ((pts j) 0) t))
      (when (= 1 ((pts j) 1)) (++ d))
      (++ j))
    (array/push blocks [t (- j i) d (- (- j i) d)])
    (set i j))
  blocks)

(defn survfit-km
  "Kaplan-Meier estimator, R survfit (stype=1): S(t) =
  prod(1 - d_t/Y_t) over event times, deaths processed before
  same-time censorings (Y includes them). Returns {:n :time :surv
  :n-risk :n-event :n-censor} - aligned arrays over ALL distinct
  times; censor-only times carry the curve forward. Opts:
  {:impl}. SPEC km-logrank.qnt."
  [times statuses &opt opts]
  # R survfit (KM, stype=1); the observed survfit output contract
  (default opts {})
  (surv-check opts times statuses)
  (def events (norm-events statuses))
  (def blocks (km-blocks times events))
  (def n (length times))
  (def factors (array/new-filled (length blocks) 1.0))
  (def n-risk (array/new-filled (length blocks) 0))
  (var remaining n)
  (eachp [bi b] blocks
    (put n-risk bi remaining)
    (when (pos? (b 2))
      (put factors bi (/ (- remaining (b 2)) remaining)))
    (-= remaining (b 1)))
  (def surv
    (if (= :num (num/resolve-impl opts))
      (let [k (num/kernels)]
        ((k 'unpack-f64) ((k 'vscan) :prod (num/stage factors k))))
      (do (var s 1.0)
          (seq [f :in factors] (set s (* s f)) s))))
  {:n n
   :time (map (fn [b] (b 0)) blocks)
   :surv surv
   :n-risk n-risk
   :n-event (map (fn [b] (b 2)) blocks)
   :n-censor (map (fn [b] (b 3)) blocks)})

(defn/typed lr-check-groups {:args [:any] :ret :any} [groups]
  "Assert exactly 2 groups; returns the sorted levels. Internal."
  (def uniq (sorted (distinct groups)))
  (unless (= 2 (length uniq))
    (errorf "stats: survdiff-logrank needs exactly 2 groups, got %d" (length uniq)))
  uniq)

(defn survdiff-logrank
  "Two-group log-rank test, R survdiff (rho = 0): O-E with the
  hypergeometric variance; chi = (O1-E1)^2/V over the pooled
  event-time blocks, p = pchisq(chi, 1, upper). Returns
  {:statistic :p-value :observed :expected :variance :n :df}.
  Throws unless exactly 2 groups. Opts: {:impl}. SPEC
  km-logrank.qnt."
  [times statuses groups &opt opts]
  # R survdiff (rho = 0, two groups); Y = Y1 + Y2 at every event time;
  # O1 + O2 == d
  (default opts {})
  (surv-check opts times statuses)
  (unless (indexed? groups)
    (error "stats: composition path consumes janet arrays; unpack packed handles first"))
  (unless (= (length groups) (length times))
    (error "stats: groups must match the times length"))
  (def events (norm-events statuses))
  (def levels (lr-check-groups groups))
  (def g-id (fn [g] (if (= g (levels 0)) 0 1)))
  (def n1 (count (fn [g] (= g (levels 0))) groups))
  # pooled blocks at EVENT times with the per-group decomposition; the
  # at-risk Y_g = #{time >= t} (remaining BEFORE this block's
  # decrement - same-time censorings INCLUDED, the spec's convention)
  (def pts (sorted (do (def acc @[])
                        (for i 0 (length times)
                          (array/push acc [(times i) (events i) (g-id (groups i))]))
                        acc)
                   (fn [p q] (< (p 0) (q 0)))))
  (def n (length pts))
  (def rows @[])
  (var i 0)
  (var rem1 n1)
  (var rem2 (- (length groups) n1))
  (while (< i n)
    (def t ((pts i) 0))
    (var j i)
    (def d @[0 0])
    (def sz @[0 0])
    (while (and (< j n) (= ((pts j) 0) t))
      (def g ((pts j) 2))
      (+= (sz g) 1)
      (when (= 1 ((pts j) 1)) (+= (d g) 1))
      (++ j))
    (when (pos? (+ (d 0) (d 1)))
      (array/push rows [(d 0) (d 1) rem1 rem2]))
    (-= rem1 (sz 0))
    (-= rem2 (sz 1))
    (set i j))
  # per-event-time terms (the janet block walk), then the O/E/V sums
  (def e1-terms (seq [r :in rows] (/ (* (+ (r 0) (r 1)) (r 2)) (+ (r 2) (r 3)))))
  (def e2-terms (seq [r :in rows] (/ (* (+ (r 0) (r 1)) (r 3)) (+ (r 2) (r 3)))))
  (def v-terms
    (seq [r :in rows]
      (def ytot (+ (r 2) (r 3)))
      (def dtot (+ (r 0) (r 1)))
      (if (< 1 ytot)
        (* dtot (/ (r 2) ytot) (/ (r 3) ytot) (/ (- ytot dtot) (- ytot 1)))
        0)))
  (defn sum-of [terms]
    (if (= :num (num/resolve-impl opts))
      (let [k (num/kernels)]
        ((k 'vreduce) :sum (num/stage terms k)))
      (do (var acc 0) (each v terms (+= acc v)) acc)))
  (def o1 (sum-of (seq [r :in rows] (r 0))))
  (def o2 (sum-of (seq [r :in rows] (r 1))))
  (def e1 (sum-of e1-terms))
  (def e2 (sum-of e2-terms))
  (def v (sum-of v-terms))
  (def chi (if (pos? v) (/ (* (- o1 e1) (- o1 e1)) v) 0))
  {:statistic chi :p-value (- 1 (pchisq chi 1 {}))
   :observed @[o1 o2] :expected @[e1 e2] :variance v
   :n @[n1 (- (length groups) n1)] :df 1})

#==------------------------------------------------------------------------==#
# One-way ANOVA + Tukey HSD. R aov/oneway.test/TukeyHSD
# conventions decoded from the r-source mirror: oneway.test.R,
# TukeyHSD.R verbatim; the aov one-way table is the orthogonal
# decomposition (SS_g = TSS - RSS, R computes it via proj/qr -
# same values within the fixture bar). Dual-impl: aov rides lm
# (the fit routes :impl); tukey-hsd's pairwise p's are ONE vector
# ptukey call (scalar k/df broadcast over q). Poles: janet
# division THROWS where R yields Inf/NaN - every reachable
# 0-denominator carries the explicit guard.
#==------------------------------------------------------------------------==#

(defn ao-sum-of
  "Impl-threaded scalar accumulation (vreduce :sum under :num).
  Internal."
  [terms &opt opts]
  (default opts {})
  # impl-threaded scalar accumulation (the survdiff sum-of pattern)
  (if (= :num (num/resolve-impl opts))
    (let [k (num/kernels)]
      ((k 'vreduce) :sum (num/stage terms k)))
    (do (var acc 0) (each v terms (+= acc v)) acc)))

(defn aov
  "One-way ANOVA, R aov(): the lm fit PLUS the ANOVA table -
  groups = array of k sample arrays; returns everything lm
  returns plus {:k :means :group-n :df-residual :mse :table},
  the table {:group {:df :ss :ms :f :p} :residual
  {:df :ss :ms}} (the orthogonal decomposition; the table F ==
  lm's overall F). R's Inf/NaN F poles are mirrored explicitly
  (janet division throws). Opts: {:impl} + lm's. Throws on
  < 2 groups, empty groups, non-finite observations or n <= k."
  [groups &opt opts]
  # groups = array of k sample arrays (the kruskal-wallis convention;
  # janet-side only - the design-matrix build is the composition
  # class). Returns the lm fit plus the ANOVA table + the tukey-hsd
  # inputs (:k :means :group-n :df-residual :mse :table).
  (default opts {})
  (def k (length groups))
  (when (< k 2) (errorf "stats: aov needs >= 2 groups, got %d" k))
  (def ns @[])
  (def ys @[])
  (def X @[])
  (eachp [gi g] groups
    (unless (indexed? g)
      (error "stats: aov consumes janet arrays per group; unpack packed handles first"))
    (when (empty? g) (error "stats: aov needs non-empty groups"))
    (array/push ns (length g))
    (each v g
      (unless (and (number? v) (< (math/abs v) math/inf))
        (error "stats: aov needs finite observations"))
      (array/push ys v)
      (def row (array/new-filled k 0))
      (put row 0 1)
      (when (pos? gi) (put row gi 1))
      (array/push X row)))
  (def n (length ys))
  (when (<= n k) (error "stats: aov needs n > k (lm's own check)"))
  (def fit (lm X ys opts))
  (def ms (seq [g :in groups] (mean g opts)))
  (def ybar (mean ys opts))
  (def ssb (ao-sum-of (seq [i :range [0 k]] (* (ns i) (* (- (ms i) ybar) (- (ms i) ybar)))) opts))
  (def ssw (ao-sum-of (seq [r :in (fit :residuals)] (* r r)) opts))
  (def df-g (- k 1))
  (def df-r (- n k))
  (def ms-g (/ ssb df-g))
  (def mse (/ ssw df-r))
  # R's F at mse = 0: Inf (unequal means) / NaN (equal means) - janet
  # division throws, so the explicit guards
  (def f (cond
           (zero? mse) (if (zero? ms-g) math/nan math/inf)
           (zero? ms-g) 0
           (/ ms-g mse)))
  (def p (- 1 (pf f df-g df-r {})))
  (merge fit
         {:k k :means ms :group-n ns :df-residual df-r :mse mse
          :table {:group {:df df-g :ss ssb :ms ms-g :f f :p p}
                  :residual {:df df-r :ss ssw :ms mse}}}))

(defn tukey-hsd
  "Tukey honestly-significant-differences, R TukeyHSD.aov
  verbatim; fit = the aov result. Returns {:pairs :conf-level},
  each pair {:i :j :label :diff :lwr :upr :p} - pair order = the
  lower-triangle COLUMN-MAJOR (k=3: 2-1, 3-1, 3-2; later level
  minus earlier, level order not sorted unless :ordered); p adj =
  1 - ptukey(|est|, k, df.res). Opts: {:conf-level 0.95 :labels
  :ordered false :impl}; throws on a malformed fit or a label
  count mismatch."
  [fit &opt opts]
  # R TukeyHSD.aov verbatim; fit = the aov result. Pair order =
  # lower-triangle COLUMN-MAJOR (k=3: 2-1, 3-1, 3-2 - later level
  # minus earlier, level order not sorted unless :ordered). p adj =
  # ptukey upper tail = 1 - ptukey(|est|, k, df.res). ONE vector
  # ptukey call for the whole table; qtukey stays
  # scalar composition (all pairs share the critical value).
  (default opts {})
  # validate the required fit keys up front
  # - a malformed struct HANGS downstream (nil k loops) otherwise
  (each key [:k :df-residual :mse :means :group-n]
    (when (nil? (get fit key))
      (errorf "stats: tukey-hsd needs an aov fit (missing :%s)" key)))
  (def k (fit :k))
  (def df-r (fit :df-residual))
  (def mse (fit :mse))
  (def conf (get opts :conf-level 0.95))
  (def raw-labels (get opts :labels (seq [i :range [0 k]] (string (+ i 1)))))
  (unless (= (length raw-labels) k)
    (errorf "stats: tukey-hsd needs %d labels, got %d" k (length raw-labels)))
  (def order
    (if (get opts :ordered false)
      # sort by means, ties by original index (R's order is stable)
      (map |(get $ 1)
           (sorted (seq [i :range [0 k]] [((fit :means) i) i])
                   (fn [a b] (or (< (a 0) (b 0)) (and (= (a 0) (b 0)) (< (a 1) (b 1)))))))
      (seq [i :range [0 k]] i)))
  (def means (map |((fit :means) $) order))
  (def ns (map |((fit :group-n) $) order))
  (def labels (map |(raw-labels $) order))
  (def qcrit (qtukey conf k df-r {}))
  (def diffs @[])
  (def ses @[])
  (def pairs @[])
  (def ests @[])
  (for j 0 (- k 1)
    (for i (+ j 1) k
      (def d (- (means i) (means j)))
      (def se (math/sqrt (* 0.5 mse (+ (/ (ns i)) (/ (ns j))))))
      # mse = 0 (groups internally constant): R's est = d/se is
      # ±Inf / NaN - the explicit guards (janet division throws)
      (def est (if (zero? se)
                 (if (zero? d) math/nan (if (pos? d) math/inf (- math/inf)))
                 (/ d se)))
      (array/push pairs {:i (order i) :j (order j)
                         :label (string (labels i) "-" (labels j))})
      (array/push diffs d)
      (array/push ses se)
      (array/push ests (math/abs est))))
  (def cdfs (ptukey ests k df-r {:impl (get opts :impl :impl)}))
  (def rows
    (map (fn [p d se c]
           (merge p {:diff d
                     :lwr (- d (* qcrit se))
                     :upr (+ d (* qcrit se))
                     :p (- 1 c)}))
         pairs diffs ses cdfs))
  {:pairs rows :conf-level conf})
