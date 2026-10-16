# bayes.janet - bayesian v1 (conjugate posteriors, MH, Gibbs,
# NIG bayes-lm, posterior toolkit, split R-hat + ESS) + EB
# beta-binomial + Nelder-Mead.
# Conformance policy: exact for analytic (the R q family; the
# flat-prior NIG limit = the R-verified lm table); moments +
# determinism for simulation; conventions OURS where no base or
# recommended R exists (rhat/ess, HPD, equal-tailed draw CI).
# draws-summary/rhat/ess moment sums ride the kernels (the
# half-chain moments; autocovariance = fused vdot on shifted slices,
# centered ONCE - the pure loop re-centers every lag); gibbs' s2-
# conditional sum kernels per step (staged once; the RNG sequence is
# untouched - draws identical, sums bitwise); MH/optimizers/bayes-lm
# are composition (scalar chain state, user fns, small-p glue).
(import typed :prefix "")
(import ./utils :prefix "")
(import ./dist :prefix "")
(import ./regsurv :prefix "")
(import ./num)

(defn arrays-only
  "Composition guard: resolves/validates the impl, then asserts
  janet-array inputs - sorts (hpd/order statistics) and the chain
  structures need janet-side values; handles unpack first."
  [opts & args]
  (num/check-impl-landed opts)
  (each x args
    (unless (indexed? x)
      (error "stats: composition path consumes janet arrays; unpack packed handles first")))
  nil)

(defn mean-of
  "Mean (the 1.0-force keeps int data honest); kernel path under
  :num. Internal."
  [xs &opt opts]
  (default opts {})
  # kernel moment under :num; the 1.0-force keeps int data honest
  (if (= :num (num/resolve-impl opts))
    (let [k (num/kernels)]
      (/ ((k 'vreduce) :sum (num/stage xs k)) (* 1.0 (num/count xs))))
    (do
      (var s 0)
      (each v xs (+= s v))
      (/ s (* 1.0 (length xs))))))

(defn var1-of
  "Sample variance, ddof=1 (R var); kernel path under :num.
  Internal."
  [xs &opt opts]
  (default opts {})
  # sample variance, ddof=1 (R var); kernel: centered pass + fused vdot
  (if (= :num (num/resolve-impl opts))
    (let [k (num/kernels)
          n (num/count xs)
          xf (num/stage xs k)
          m (/ ((k 'vreduce) :sum xf) (* 1.0 n))
          c ((k 'vmap1-scalar) :sub xf m)]
      (/ ((k 'vdot) c c) (* 1.0 (- n 1))))
    (do
      (def n (length xs))
      (def m (mean-of xs))
      (var s 0)
      (each v xs (+= s (* (- v m) (- v m))))
      (/ s (* 1.0 (- n 1))))))

(defn/typed log-choose {:args [:number :number] :ret :number} [n k]
  "log C(n,k) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)."
  (- (lgamma (inc n)) (lgamma (inc k)) (lgamma (inc (- n k)))))

#==------------------------------------------------------------------------==#
# conjugate posteriors - analytic; equal-tailed CIs at opts :conf.
#==------------------------------------------------------------------------==#

(defn post-beta-binomial
  "Beta(a,b) prior on the success probability, k successes of n
  trials: posterior Beta(a+k, b+n-k). Returns {:a :b :mean :sd :ci
  :log-marginal} - the equal-tailed CI at :conf (default 0.95) and
  the FULL beta-binomial log-marginal (log C(n,k) included, the
  VGAM dbetabinom.ab convention, so log-marginals compare across
  languages). Opts: {:conf :impl}."
  [k n a b &opt opts]
  (default opts {})
  # closed-form conjugate scalars
  (num/check-impl-landed opts)
  (def conf (get opts :conf 0.95))
  (def pa (+ a k))
  (def pb (+ b (- n k)))
  (def lo (qbeta (/ (- 1 conf) 2) pa pb {}))
  (def hi (qbeta (/ (+ 1 conf) 2) pa pb {}))
  {:a pa :b pb
   :mean (/ pa (* 1.0 (+ pa pb)))
   :sd (math/sqrt (/ (* pa pb) (* 1.0 (+ pa pb) (+ pa pb 1))))
   :ci [lo hi]
   # log P(k) = log C(n,k) + lbeta(a+k, b+n-k) - lbeta(a,b) - the FULL
   # beta-binomial PMF (VGAM dbetabinom.ab convention; binomial
   # coefficient INCLUDED so log-marginals compare across languages)
   :log-marginal (+ (log-choose n k)
                    (- (log-beta pa pb) (log-beta a b)))})

(defn post-normal-normal
  "Normal N(mu0, tau2) prior on the mean, known data variance
  sigma2: posterior over the mean given observations xs. Returns
  {:mean :var :sd :ci :log-marginal} - the equal-tailed CI at :conf
  and the marginal likelihood of the data via xbar ~ N(mu0,
  tau2 + sigma2/n). Opts: {:conf :impl}."
  [xs mu0 tau2 sigma2 &opt opts]
  (default opts {})
  (num/check-impl-landed opts xs)
  (def conf (get opts :conf 0.95))
  (def n (length xs))
  (def xbar (mean-of xs opts))
  (def prec (+ (/ 1 tau2) (/ n sigma2)))
  (def v2 (/ 1 prec))
  (def m2 (* v2 (+ (/ mu0 tau2) (/ (* n xbar) sigma2))))
  {:mean m2 :var v2 :sd (math/sqrt v2)
   :ci [(- m2 (* (qnorm (/ (- 1 conf) 2) 0 1 {}) (math/sqrt v2)))
        (+ m2 (* (qnorm (/ (+ 1 conf) 2) 0 1 {}) (math/sqrt v2)))]
   # xbar ~ N(mu0, tau2 + sigma2/n): log marginal likelihood of the data
   :log-marginal (dnorm xbar mu0 (math/sqrt (+ tau2 (/ sigma2 n))) {:log true})})

(defn post-gamma-poisson
  "Gamma(a,b) prior (rate parameterisation) on the Poisson rate:
  posterior Gamma(a + sum ks, b + sum exposure) with ks counts over
  exposures (a scalar exposure applies to every count). Returns
  {:shape :rate :mean :sd :ci} - the equal-tailed CI at :conf.
  Opts: {:conf :impl}."
  [ks exposure a b &opt opts]
  (default opts {})
  # the data sum kernels under :num (handles chain); exposure
  # stays scalar/array glue
  (num/check-impl-landed opts ks)
  (def conf (get opts :conf 0.95))
  (def sk
    (if (= :num (num/resolve-impl opts))
      (let [k (num/kernels)]
        ((k 'vreduce) :sum (num/stage ks k)))
      (do
        (var s 0)
        (each k ks (+= s k))
        s)))
  (def se (if (indexed? exposure)
           (do (var s 0) (each e exposure (+= s e)) s)
           (* exposure (* 1.0 (length ks)))))
  (def shape (+ a sk))
  (def rate (+ b se))
  {:shape shape :rate rate
   :mean (/ shape (* 1.0 rate))
   :sd (math/sqrt (/ shape (* 1.0 (* rate rate))))
   :ci [(qgamma (/ (- 1 conf) 2) shape rate {})
        (qgamma (/ (+ 1 conf) 2) shape rate {})]})

#==------------------------------------------------------------------------==#
# MH - random-walk Metropolis-Hastings, 1-D, NORMAL proposal (inversion
# via qnorm). Adaptive: Robbins-Monro on log proposal sd during burn-in
# toward 0.44 acceptance (1-D optimal band). Deterministic per :seed.
# Conventions OURS (documented); conformance = moments + determinism.
#==------------------------------------------------------------------------==#

(defn mh
  "1-D random-walk Metropolis-Hastings with a NORMAL proposal and
  optional Robbins-Monro adaptation of the proposal sd toward
  0.44 acceptance during burn-in. log-post is the user log-posterior
  fn; returns {:draws :accept-rate :post-burn-accept} - draws of
  length steps - burn. Deterministic per :seed. Opts: {:steps 5000
  :burn 1000 :adaptive true :prop-sd 1 :seed 'stats :impl}. SPEC
  mc-chain.qnt."
  [log-post init &opt opts]
  (default opts {})
  # 1-D random-walk chain - scalar state, user
  # log-posterior; the loop IS the algorithm
  (num/check-impl-landed opts)
  (def steps (get opts :steps 5000))
  (def burn (get opts :burn 1000))
  (def adaptive (get opts :adaptive true))
  (def rng (math/rng (string (get opts :seed "stats"))))
  (var x init)
  (var lx (log-post x))
  (var log-sd (math/log (get opts :prop-sd 1)))
  (var accepts 0)
  (var acc-burn 0)
  (def draws (array/new (max 0 (- steps burn))))
  (var acc-post 0)
  (var npost 0)
  (var z-spare nil)
  (for i 0 steps
    (when (nil? z-spare) (set z-spare (rnorm-pair rng)))
    (def z (if (number? z-spare) z-spare (z-spare 0)))
    (if (number? z-spare)
      (set z-spare nil)
      (set z-spare (z-spare 1)))
    (def prop (+ x (* (math/exp log-sd) z)))
    (def lprop (log-post prop))
    (when (< (math/log (math/rng-uniform rng)) (- lprop lx))
      (set x prop)
      (set lx lprop)
      (++ accepts)
      (when (< i burn) (++ acc-burn))
      (when (>= i burn) (++ acc-post)))
    (when (and adaptive (< i burn))
      # Robbins-Monro on log sd toward 0.44 acceptance during burn-in
      (set log-sd (- log-sd (* 0.05 (- 0.44 (/ acc-burn (* 1.0 (+ i 1))))))))
    (when (>= i burn)
      (++ npost)
      (array/push draws x)))
  {:draws draws :accept-rate (/ accepts (* 1.0 steps))
   :post-burn-accept (if (pos? npost) (/ acc-post (* 1.0 npost)) 0)})

#==------------------------------------------------------------------------==#
# Gibbs - normal model, unknown mean + variance (canonical 2-block).
# Prior: mu|s2 ~ N(m0, s2/k0); s2 ~ InvGamma(a0/2, b0/2).
# Conditionals: mu|s2,x ~ N(m_n, s2/k_n), k_n = k0+n,
# m_n = (k0*m0 + n*xbar)/k_n; s2|mu,x ~ InvGamma((a0+n+1)/2,
# (b0 + k0*(mu-m0)^2 + sum (x-mu)^2)/2) - the mu-prior term k0(mu-m0)^2
# IS part of the s2-conditional because the prior mu-variance scales
# with s2 (without it E[s2] lands at 0.36 instead of 0.69, caught vs
# numerical integration). Drawn via Gamma + reciprocal.
#==------------------------------------------------------------------------==#

(defn gibbs-normal
  "Gibbs sampler for the normal model with unknown mean and
  variance (canonical 2-block): prior mu|s2 ~ N(m0, s2/k0), s2 ~
  InvGamma(a0/2, b0/2). Returns {:mu-draws :s2-draws :mu-mean
  :mu-sd :s2-mean}. Deterministic per :seed. Opts: {:steps 4000
  :burn 1000 :m0 0 :k0 1 :a0 2 :b0 1 :seed 'stats :impl}. SPEC
  mc-chain.qnt."
  [xs &opt opts]
  (default opts {})
  (num/check-impl-landed opts xs)
  (def steps (get opts :steps 4000))
  (def burn (get opts :burn 1000))
  (def rng (math/rng (string (get opts :seed "stats"))))
  (def n (length xs))
  (def xbar (mean-of xs opts))
  (def m0 (get opts :m0 0))
  (def k0 (get opts :k0 1))
  (def a0 (get opts :a0 2))
  (def b0 (get opts :b0 1))
  (var mu m0)
  (var s2 (/ b0 (* 1.0 a0)))
  (def mu-draws (array/new (max 0 (- steps burn))))
  (def s2-draws (array/new (max 0 (- steps burn))))
  (def kn (+ k0 n))
  (def mn (/ (+ (* k0 m0) (* n xbar)) kn))
  (var z-spare nil)
  # the s2-conditional sum - kernel branch stages xs ONCE and
  # centers per step (2 kernel calls vs n interpreted mul-adds); the
  # RNG sequence is untouched either way - draws identical, the sum
  # bitwise (same products, same order)
  (def sumsq
    (if (= :num (num/resolve-impl opts))
      (let [k (num/kernels)
            xf (num/stage xs k)]
        (fn [mu]
          (def c ((k 'vmap1-scalar) :sub xf mu))
          (+ ((k 'vdot) c c) (* k0 (- mu m0) (- mu m0)))))
      (fn [mu]
        (var ss 0)
        (for j 0 n (+= ss (* (- (get xs j) mu) (- (get xs j) mu))))
        (+ ss (* k0 (- mu m0) (- mu m0))))))
  (for i 0 steps
    (when (nil? z-spare) (set z-spare (rnorm-pair rng)))
    (def z (if (number? z-spare) z-spare (z-spare 0)))
    (if (number? z-spare)
      (set z-spare nil)
      (set z-spare (z-spare 1)))
    (set mu (+ mn (* z (math/sqrt (/ s2 kn)))))
    (set s2 (/ 1 (rgamma-one rng (/ (+ a0 n 1) 2) (/ (+ b0 (sumsq mu)) 2))))
    (when (>= i burn)
      (array/push mu-draws mu)
      (array/push s2-draws s2)))
  {:mu-draws mu-draws :s2-draws s2-draws
   :mu-mean (mean-of mu-draws opts) :mu-sd (math/sqrt (var1-of mu-draws opts))
   :s2-mean (mean-of s2-draws opts)})

#==------------------------------------------------------------------------==#
# Bayes linear regression - conjugate NIG(m0, L0^-1, a0/2, b0/2):
# beta|s2 ~ N(m0, s2*L0^-1) (L0 = prior PRECISION matrix),
# s2 ~ InvGamma(a0/2, b0/2). Posterior: Ln = L0 + X'X;
# mn = Ln^-1 (L0 m0 + X'y); an = a0 + n;
# bn = b0 + y'y + m0'L0 m0 - mn'Ln mn. Marginal of beta_j:
# t_{an}, scale sqrt((bn/an) * (Ln^-1)_jj) - :se is that scale (the
# lm-table quantity), :marginal-sd the t sd sqrt(an/(an-2)) * se.
# DEFAULT = flat/Jeffreys limit (l0 = zero matrix, m0 = 0, a0 = -p,
# b0 = 0): reproduces the classical OLS table EXACTLY (R-verified:
# coef, se, t, p, CI, df = n - p, RSS scale).
#==------------------------------------------------------------------------==#

(defn bayes-lm
  "Conjugate NIG Bayes linear regression; the flat/Jeffreys default
  reproduces the classical OLS table exactly (coef, se, t, p, CI, df
  = n - p). X = n x p row-major design; returns {:beta-mean :se :t
  :p :ci :df :sigma2-mean :marginal-sd}. Opts: {:m0 :l0 :a0 :b0
  :conf 0.95 :impl} - :l0 is the p x p prior PRECISION matrix.
  Throws when n <= p."
  [X y &opt opts]
  (default opts {})
  # the p x p NIG machinery is small-p glue (the
  # GEMMs live in lm/glm); the janet matrix code IS the algorithm
  (num/check-impl-landed opts)
  (def conf (get opts :conf 0.95))
  (def n (m-rows X))
  (def p (m-cols X))
  (when (<= n p) (error "stats: bayes-lm needs n > p"))
  (def ones (array/new-filled n 1))
  (def xtx (xtwx X ones))
  (def xty (xtwv X ones y))
  (var yty 0)
  (each v y (+= yty (* v v)))
  (def m0 (get opts :m0 (array/new-filled p 0)))
  (def l0 (get opts :l0))
  (def a0 (get opts :a0 (- p)))
  (def b0 (get opts :b0 0))
  (when (and l0 (not= (length l0) p))
    (error "stats: bayes-lm :l0 must be p x p"))
  # Ln = L0 + X'X; mn = Ln^-1 (L0 m0 + X'y)
  (def ln-full
    (if l0
      (let [acc (map |(array ;$) l0)]
        (for i 0 p
          (for j 0 p
            (put (get acc i) j (+ (get (get acc i) j) (get (get xtx i) j)))))
        acc)
      xtx))
  (def rhs (if l0 (mmul l0 m0) (array/new-filled p 0)))
  (for j 0 p (put rhs j (+ (get rhs j) (get xty j))))
  (def mn (solve-sp ln-full rhs))
  (def an (+ a0 n))
  (def lmn-mn (matvec ln-full mn))
  (var qf 0)
  (for j 0 p (+= qf (* (get mn j) (get lmn-mn j))))
  (var bn (+ b0 (- yty qf)))
  (when l0
    (def l0-m0 (matvec l0 m0))
    (var q2 0)
    (for j 0 p (+= q2 (* (get m0 j) (get l0-m0 j))))
    (+= bn q2))
  (def inv (inverse-sp ln-full))
  (def s2n (/ bn (* 1.0 an)))
  (def ses (array/new p))
  (def tvals (array/new p))
  (def pvals (array/new p))
  (def ci-lo (array/new p))
  (def ci-hi (array/new p))
  (def marg-sd (array/new p))
  (def tcrit (qt (- 1 (* 0.5 (- 1 conf))) an {}))
  (for j 0 p
    (def se (math/sqrt (* s2n (get (get inv j) j))))
    (put ses j se)
    (put tvals j (/ (get mn j) se))
    (put pvals j (* 2 (- 1 (pt (math/abs (/ (get mn j) se)) an {}))))
    (put ci-lo j (- (get mn j) (* tcrit se)))
    (put ci-hi j (+ (get mn j) (* tcrit se)))
    (put marg-sd j (* se (math/sqrt (/ an (* 1.0 (- an 2)))))))
  {:beta-mean mn :se ses :t tvals :p pvals :ci [ci-lo ci-hi]
   :df an :sigma2-mean (/ bn (* 1.0 (- an 2)))
   :marginal-sd marg-sd})

#==------------------------------------------------------------------------==#
# posterior toolkit - moments from draws, equal-tailed CI (empirical
# order statistics), HPD (shortest window of ceil(conf*n) order stats;
# ties -> lowest). Conventions OURS (documented).
#==------------------------------------------------------------------------==#

(defn hpd
  "Highest-posterior-density interval of raw draws: the shortest
  window of ceil(conf*n) order statistics (ties resolve to the
  lowest window). Returns a [lo hi] tuple; throws on fewer than 2
  draws."
  [draws conf]
  (def xs (sorted draws))
  (def n (length xs))
  (when (< n 2) (error "stats: hpd needs >= 2 draws"))
  (def m (max 2 (math/ceil (* conf n))))
  (var best-lo (get xs 0))
  (var best-hi (get xs (dec m)))
  (var best-w (- best-hi best-lo))
  (for i 0 (- (inc n) m)
    (def lo (get xs i))
    (def hi (get xs (+ i (dec m))))
    (def w (- hi lo))
    (when (< w best-w)
      (set best-w w)
      (set best-lo lo)
      (set best-hi hi)))
  [best-lo best-hi])

(defn draws-summary
  "Posterior summary of raw draws: {:mean :sd :n :ci :hpd} - the
  equal-tailed CI from empirical order statistics and the HPD
  window. Opts: {:conf 0.95 :impl}; throws on fewer than 2 draws."
  [draws &opt opts]
  (default opts {})
  # moment sums via the kernels; the order-statistic CI + hpd
  # sorts stay janet - arrays only
  (arrays-only opts draws)
  (def conf (get opts :conf 0.95))
  (def xs (sorted draws))
  (def n (length xs))
  (when (< n 2) (error "stats: draws-summary needs >= 2 draws"))
  {:mean (mean-of draws opts)
   :sd (math/sqrt (var1-of draws opts))
   :n n
   :ci [(get xs (max 0 (math/floor (* n (/ (- 1 conf) 2)))))
        (get xs (min (dec n) (math/ceil (* n (/ (+ 1 conf) 2)))))]
   :hpd (hpd draws conf)})

#==------------------------------------------------------------------------==#
# convergence - SPLIT R-hat (BDA3) and ESS (Stan's Geyer initial
# monotone sequence). Conventions OURS (no base/recommended R exists);
# verified on synthetic chains with known truth: iid -> rhat ~ 1,
# ess ~ m*n; AR(1) rho -> ess ~ m*n*(1-rho)/(1+rho).
#==------------------------------------------------------------------------==#

(defn/typed halves {:args [:any] :ret :array} [chains]
  "Each chain split into two halves (odd length: first half
  shorter). Internal to rhat/ess."
  (def out @[])
  (each c chains
    (def n (length c))
    (def h (div n 2))
    (array/push out (slice c 0 h))
    (array/push out (slice c h (+ h h))))
  out)

(defn b-and-w
  "[B W N m] from half-chains (equal length, >= 2 halves) - the
  split-R-hat between/within variances. Internal."
  [hs &opt opts]
  # per-chain moments thread the impl
  (default opts {})
  (def m (length hs))
  (when (< m 2) (error "stats: need >= 2 half-chains"))
  (def n0 (length (get hs 0)))
  (each c hs (unless (= (length c) n0)
               (error "stats: chains must have equal even-ish lengths")))
  (def means (map |(mean-of $ opts) hs))
  (def vars (map |(var1-of $ opts) hs))
  (def grand (mean-of means))
  (var b 0)
  (each mu means (+= b (* (- mu grand) (- mu grand))))
  [(set b (* n0 (/ b (* 1.0 (- m 1))))) (mean-of vars) n0 m])

(defn rhat
  "Split R-hat (BDA3): each chain halved, the between/within
  variance ratio sqrt'd. chains = array of draw arrays. Returns inf
  for constant chains that DISAGREE across; throws on zero
  within-chain variance. Opts: {:impl}."
  [chains &opt opts]
  (default opts {})
  (arrays-only opts chains)
  (when (empty? chains) (error "stats: rhat needs chains"))
  (def [b w n0 m] (b-and-w (halves chains) opts))
  (cond
    # constant chains that DISAGREE across: maximally diverged
    (and (zero? w) (pos? b)) math/inf
    (zero? w) (error "stats: rhat: zero within-chain variance")
    (do
      (def var-plus (+ (* (/ (- n0 1) n0) w) (/ b n0)))
      (math/sqrt (/ var-plus w)))))

(defn autocov-biased
  "Biased (1/N) autocovariance at lags 0..tmax. num: the series is
  centered ONCE (the pure loop re-subtracts m at every lag -
  O(n*tmax) subtractions), then each lag is a fused vdot on
  shifted slices (bitwise parity). Opts: {:impl}."
  [xs tmax &opt opts]
  (default opts {})
  # the series is centered ONCE (the pure loop re-subtracts m at
  # every lag - O(n*tmax) subtractions), then each lag is a fused
  # vdot on shifted slices (slice+pack per lag is memcpy-class; the
  # win is C multiply-accumulate over the interpreted loop; same
  # products in the same order - bitwise parity)
  (def n (length xs))
  (def m (mean-of xs opts))
  (def out (array/new-filled (inc tmax) 0))
  (if (= :num (num/resolve-impl opts))
    (let [k (num/kernels)
          xc ((k 'unpack-f64) ((k 'vmap1-scalar) :sub (num/stage xs k) m))]
      (for t 0 (inc tmax)
        (def np (- n t))
        (put out t (/ ((k 'vdot) (num/stage (slice xc 0 np) k)
                                   (num/stage (slice xc t n) k))
                      (* 1.0 n)))))
    (for t 0 (inc tmax)
      (var s 0)
      (for i 0 (- n t) (+= s (* (- (get xs i) m) (- (get xs (+ i t)) m))))
      (put out t (/ s (* 1.0 n)))))
  out)

(defn ess
  "Effective sample size, Stan's Geyer initial monotone sequence
  estimator on split half-chains: tau = -1 + 2*sum(Gamma_k) over
  positive running-min pair sums, ESS = m*n/tau (floored at m*n).
  chains = array of draw arrays; throws on zero variance.
  Opts: {:impl}."
  [chains &opt opts]
  (default opts {})
  (arrays-only opts chains)
  (when (empty? chains) (error "stats: ess needs chains"))
  (def hs (halves chains))
  (def [b w n0 m] (b-and-w hs opts))
  (def var-plus (+ (* (/ (- n0 1) n0) w) (/ b n0)))
  (when (or (zero? var-plus) (zero? w)) (error "stats: ess: zero variance"))
  (def tmax (dec n0))
  (def acovs (map |(autocov-biased $ tmax opts) hs))
  # rho_t = 1 - (W - mean_j acov_{j,t}) / var+
  (def rho (seq [t :range [0 (inc tmax)]]
              (- 1 (/ (- w (mean-of (map |(get $ t) acovs))) var-plus))))
  # Geyer initial-positive + running-min over pair sums
  # Gamma_k = rho_{2k} + rho_{2k+1}; stop at first non-positive Gamma.
  # tau = -1 + 2 * sum(Gamma_k), floored at 1.
  (var gsum 0)
  (var prev math/inf)
  (var k 0)
  (var done false)
  (while (not done)
    (def i2 (* 2 k))
    (def i2p1 (inc i2))
    (if (>= i2p1 (length rho))
      (set done true)
      (let [g (+ (get rho i2) (get rho i2p1))]
        (if (<= g 0)
          (set done true)
          (do
            (def gm (min g prev))
            (set prev gm)
            (+= gsum gm)
            (++ k))))))
  (def tau (max 1 (- (* 2 gsum) 1)))
  (/ (* (* 1.0 m) n0) tau))

#==------------------------------------------------------------------------==#
# Nelder-Mead - standard (alpha=1, gamma=2, rho=0.5, sigma=0.5);
# initial simplex x0 and x0 + step*e_i; convergence when
# |f_worst - f_best| <= tol*(1+|f_best|). Conformance with base R
# optim(method='Nelder-Mead'): optimum agreement, not path agreement.
#==------------------------------------------------------------------------==#

(defn/typed nm-score {:args [:any :any] :ret :array} [f simplex]
  "Simplex scored and sorted by f, carrying the ORIGINAL index (the
  array itself is never reordered - replacements must target the
  worst INDEX). Internal to nelder-mead."
  (sorted (seq [i :range [0 (length simplex)]]
           {:f (f (get simplex i)) :v (get simplex i) :i i})
       (fn [a b] (< (a :f) (b :f)))))

(defn nelder-mead
  "Standard Nelder-Mead minimization (alpha=1, gamma=2, rho=0.5,
  sigma=0.5); the initial simplex is x0 and x0 + step*e_i.
  Returns {:x :f :iters :converged} - never errors on
  non-convergence (the flag family). Opts: {:step 1 :tol 1e-12
  :max-it 500 :impl}."
  [f x0 &opt opts]
  (default opts {})
  # the simplex walk is scalar state + a user fn
  (num/check-impl-landed opts)
  (def step (get opts :step 1))
  (def tol (get opts :tol 1e-12))
  (def max-it (get opts :max-it 500))
  (def p (length x0))
  (when (< p 1) (error "stats: nelder-mead needs a non-empty start"))
  (def simplex @[(array ;x0)])
  (for j 0 p
    (def v (array ;x0))
    (put v j (+ (get v j) step))
    (array/push simplex v))
  (var iters 0)
  (var converged false)
  (while (and (not converged) (< iters max-it))
    (++ iters)
    (def scored (nm-score f simplex))
    (def fb ((get scored 0) :f))
    (def worst (last scored))
    (def fw (worst :f))
    (when (<= (math/abs (- fw fb)) (* tol (+ 1 (math/abs fb))))
      (set converged true)
      (break))
    (def n (length scored))
    (def best-v (get (get scored 0) :v))
    (def worst-v (worst :v))
    (def fw2 ((get scored (- n 2)) :f))
    (def cen (array/new-filled p 0))
    (for i 0 (- n 1)
      (for j 0 p (put cen j (+ (get cen j) (get (get (get scored i) :v) j)))))
    (for j 0 p (put cen j (/ (get cen j) (* 1.0 (- n 1)))))
    (defn edge-point [a b frac]
      # cen + frac*(a - b)
      (def v (array/new p))
      (for j 0 p (put v j (+ (get cen j) (* frac (- (get a j) (get b j))))))
      v)
    (def refl (edge-point cen worst-v 1)) # 2*cen - worst
    (def fr (f refl))
    (defn replace-worst [v] (put simplex (worst :i) v))
    (cond
      (and (>= fr fb) (< fr fw2))
      (replace-worst refl)
      (< fr fb)
      (let [expv (edge-point refl cen 2) # cen + 2*(refl - cen)
            fe (f expv)]
        (replace-worst (if (> fe fr) expv refl)))
      (let [cc (if (< fr fw)
               (edge-point refl cen 0.5)     # outside
               (edge-point worst-v cen 0.5)) # inside
            fc (f cc)]
        (if (< fc (min fr fw))
          (replace-worst cc)
          (do # shrink toward best
            (for i 1 (length simplex)
              (def v (array/new p))
              (for j 0 p
                (put v j (+ (get best-v j) (* 0.5 (- (get (get simplex i) j) (get best-v j))))))
              (put simplex i v)))))))
  (def best (get (nm-score f simplex) 0))
  {:x (best :v) :f (best :f) :iters iters :converged converged})

#==------------------------------------------------------------------------==#
# EB beta-binomial - groups k_i ~ Bin(n_i, pi_i), pi_i ~ Beta(a,b);
# (a,b) by marginal-likelihood MLE (Nelder-Mead on [log a, log b]),
# then per-group posterior Beta(a+k_i, b+n_i-k_i).
#==------------------------------------------------------------------------==#

(defn eb-beta-binomial
  "Empirical-Bayes beta-binomial: (a, b) fitted by marginal-likelihood
  MLE over the groups (ks/ns = per-group successes/trials), then a
  per-group Beta posterior. Returns {:a :b :log-marginal :converged
  :groups} - each group {:k :n :a :b :mean :ci}. Opts: {:a0 1 :b0 1
  :conf 0.95 :impl}; throws on length mismatch."
  [ks ns &opt opts]
  (default opts {})
  # the EM loop is scalar state per group
  (num/check-impl-landed opts)
  (when (not= (length ks) (length ns))
    (error "stats: eb-beta-binomial needs equal-length ks and ns"))
  (def a0 (get opts :a0 1))
  (def b0 (get opts :b0 1))
  (defn neg-log-marginal [ab]
    # ab = [log a, log b] - log scale keeps the search positive. The
    # objective is the FULL beta-binomial log-likelihood (log C(n,k)
    # included - same convention as the PMF and as ebbr/VGAM; the MLE
    # is unaffected by the constant)
    (def a (math/exp (get ab 0)))
    (def b (math/exp (get ab 1)))
    (var ll 0)
    (for i 0 (length ks)
      (def k (get ks i))
      (def n (get ns i))
      (+= ll (+ (log-choose n k)
                (- (log-beta (+ a k) (+ b (- n k))) (log-beta a b)))))
    (- ll))
  (def res (nelder-mead neg-log-marginal [(math/log a0) (math/log b0)]
                        {:step 0.5 :tol 1e-12 :max-it 1000}))
  (def a (math/exp (get (res :x) 0)))
  (def b (math/exp (get (res :x) 1)))
  (def conf (get opts :conf 0.95))
  (def groups
    (seq [i :range [0 (length ks)]]
      (def post (post-beta-binomial (get ks i) (get ns i) a b {:conf conf}))
      {:k (get ks i) :n (get ns i)
       :a (post :a) :b (post :b) :mean (post :mean) :ci (post :ci)}))
  {:a a :b b :log-marginal (- (res :f)) :converged (res :converged)
   :groups groups})
