# suite-stats.janet - the janet-stats conformance suite (testament;
# stats is pure janet, no libjn required). The spec-mirror home:
# BOTH forward contracts ACTIVATED here (ecdf-ks, km-logrank).
# Default :impl resolution throughout (the machine-dynamic path);
# explicit :pure/:num rows where the contract distinguishes them.
(unless (get (os/environ) "JANET_PATH")
  (setdyn *syspath* (string (or (get (os/environ) "HOME") (error "HOME"))
                            "/.local/lib/janet")))
(import typed :prefix "")
(import kit :prefix "")
(import testament :prefix "")
(import stats/num)
(import stats/tests)
(import stats/utils)
(import stats/ranksum)
(import stats/descriptive)
(import stats/dist)
(import stats/regsurv)
(import stats/resample)
(import stats/bayes)

(defn close [a b]
  #  (the close() inf-degeneracy audit): a finite-pinned
  # row passed vacuously on an inf actual ((<= inf inf) via the rel
  # clause). Finiteness guards first - inf/NaN actuals fail; rows
  # that pin infinity use (= v math/inf), never close.
  (and (< (math/abs a) math/inf) (< (math/abs b) math/inf)
       (or (<= (math/abs (- a b)) 1e-12)
           (<= (math/abs (- a b)) (* 1e-12 (math/abs b))))))

#==------------------------------------------------------------------------==#
# SPEC MIRRORS - the Quint contracts executed against the real code
# (the impl-resolution pattern). Retroactive (moved from the
# libjn-gated janet-num suite, 2026-10-04): iterative-solver,
# mc-chain, perm-enumeration. FORWARD (activated at landing): ecdf-ks
# and km-logrank.
#==------------------------------------------------------------------------==#

(deftest solver-error-paths
  # iterative-solver.qnt mirror: the error-terminating family (glm/
  # coxph/clogit) fails with EXACTLY the two defined causes; results
  # only when converged. NM is the flag family - never errors on budget
  (put num/state :presence :present)
  (def X (seq [i :range [0 40]] @[1.0 (math/sin (* i 0.1)) (math/cos (* i 0.05))]))
  (def y (seq [i :range [0 40]] (+ (* 2 (math/sin (* i 0.1))) (* 0.1 i) 1)))
  (var caught nil)
  (try (regsurv/glm X y :poisson {:max-it 1}) ([e] (set caught e)))
  (is (and caught (string/find "did not converge in 1" caught))
      "glm max-it 1 -> the iteration-bound error, actionable")
  (var caught2 nil)
  (try (regsurv/coxph X y (seq [i :range [0 40]] (if (even? i) 1 0)) {:max-it 1})
    ([e] (set caught2 e)))
  (is (and caught2 (string/find "did not converge in 1" caught2))
      "coxph max-it 1 -> the iteration-bound error")
  (var caught3 nil)
  (try (regsurv/clogit X (seq [i :range [0 40]] (if (> (y i) 2) 1 0))
                       (seq [i :range [0 40]] (% i 4)) {:max-it 1})
    ([e] (set caught3 e)))
  (is (and caught3 (string/find "did not converge" caught3))
      "clogit max-it 1 -> the iteration-bound error")
  # the flag family: NM returns :converged, never errors on budget
  (var nm-caught nil)
  (def nm (try (bayes/nelder-mead (fn [x] (+ (* 1.0 (get x 0) (get x 0)) 1.0)) @[2.0] {:max-it 1})
             ([e] (set nm-caught e) nil)))
  (is (and (nil? nm-caught) (not (nm :converged)))
      "NM non-convergence is a FLAG, not an error")
  # converged runs: iters within the bound, results present
  (def r (regsurv/glm X y :poisson {:max-it 25}))
  (is (and (>= (r :iters) 1) (<= (r :iters) 25) (r :converged))
      "converged glm: iters within the bound"))

(deftest chain-burn-boundary
  # mc-chain.qnt mirror: draws appended ONLY post-burn, exactly one
  # per iteration; accept accounting in [0,1]
  (def lp (fn [x] (* -0.5 x x)))
  (def m (bayes/mh lp 0 {:steps 120 :burn 30 :seed "mb"}))
  (is (= 90 (length (m :draws))) "draws == steps - burn")
  (is (and (>= (m :post-burn-accept) 0) (<= (m :post-burn-accept) 1))
      "post-burn accept in [0,1]")
  (def m0 (bayes/mh lp 0 {:steps 50 :burn 0 :seed "mb"}))
  (is (= 50 (length (m0 :draws))) "burn 0 -> every step draws")
  (def g (bayes/gibbs-normal @[5.1 4.9 5.2 5.0 4.8 5.05]
                            {:steps 60 :burn 10 :seed "gb"}))
  (is (and (= 50 (length (g :mu-draws))) (= 50 (length (g :s2-draws))))
      "gibbs draws == steps - burn, both chains"))

(deftest perm-two-lane-contract
  # perm-enumeration.qnt mirror: exact lane B == C(n1+n2,n1) (every
  # split once), p in [0,1]; MC add-one p in (0,1] (OUR convention)
  (def a @[1 2 3 4])
  (def b @[5 6 7])
  (def r (resample/perm-test-two a b {}))
  (is (r :exact) "auto-exact for tiny samples")
  (is (= 35 (r :B)) "B == C(7,4) = 35 - all splits, once")
  (is (and (>= (r :p-value) 0) (<= (r :p-value) 1)) "exact p in [0,1]")
  (def m (resample/perm-test-two a b {:exact false :B 19 :seed "pb"}))
  (is (and (> (m :p-value) 0) (<= (m :p-value) 1))
      "MC add-one p in (0,1]"))

#==------------------------------------------------------------------------==#
# ecdf-ks.qnt - FORWARD CONTRACT ACTIVATED
# The dormant mirror's fn names finalized: stats/ecdf (proportions,
# the R convention - the draft's 2.0 values were the spec's scaled
# numerators), stats/ks-test (struct result; D at :statistic).
#==------------------------------------------------------------------------==#

(deftest ecdf-fold-contract
  # ecdf-ks.qnt invariants: F monotone, in [0,1]; TIES SHARE ONE VALUE
  # (full-block jumps - the fold consumed each tie block whole);
  # F(max) == 1 at exhaustion
  (def F (descriptive/ecdf @[1 1 2 3 3 3] {}))
  (is (close (/ 2 6) (F 1)) "ties share one value: F(1) = 2/6")
  (is (close 0.5 (F 2)) "F(2) = 3/6")
  (is (= 1 (F 3)) "F(3) = 6/6 at the max block")
  (is (= 0 (F 0.5)) "below the data F = 0")
  (is (= 1 (F 10)) "above the data F = 1")
  (is (deep= @[0 (/ 2 6) 1] (F @[0.5 1 3])) "grid evaluation (R's vectorized F)"))

(deftest ks-statistic-contract
  # ecdf-ks.qnt: D over the UNION grid, bounded [0,1], SYMMETRIC in its
  # arguments; the spec's concrete pair A=[1,1,2] B=[1,3,3] -> D = 2/3
  (def d1 ((tests/ks-test @[1 1 2] @[1 3 3] {}) :statistic))
  (is (close (/ 2 3) d1) "the concrete pair D = 2/3")
  (is (= d1 ((tests/ks-test @[1 3 3] @[1 1 2] {}) :statistic))
      "D symmetric in its arguments")
  (is (and (>= d1 0) (<= d1 1)) "D in [0,1]"))

(deftest ks-exact-switch-decode
  # the decoded switches (R ks.test): two-sample
  # nx*ny < 10000; one-sample n < 100 && no ties (ties force the
  # asymptotic lane one-sample; the two-sample exact is tie-aware)
  (def a (seq [i :range [0 99]] (math/sin (* i 0.37))))
  (def b (seq [i :range [0 100]] (math/cos (* i 0.29))))
  (is ((tests/ks-test a b {}) :exact) "99*100 = 9900 < 10000 -> exact")
  (def a2 (seq [i :range [0 100]] (math/sin (* i 0.37))))
  (is (not ((tests/ks-test a2 b {}) :exact)) "100*100 = 10000 -> asymptotic")
  (def z (seq [i :range [0 99]] (math/sin (* i 0.41))))
  (is ((tests/ks-test z :pnorm {}) :exact) "one-sample n=99 < 100 -> exact")
  (def zt (map (fn [v] (/ (math/round (* 10 v)) 10)) z))
  (is (not ((tests/ks-test zt :pnorm {}) :exact))
      "one-sample ties force the asymptotic lane"))

(deftest ks-guards-and-tie-saturation
  # the psmirnov/pkolmogorov guards + tie saturation (all expected
  # values verified vs R 2026-10-04)
  (def r1 (tests/ks-test @[1.5 2.5] @[1.5 2.5] {}))
  (is (and (zero? (r1 :statistic)) (= 1 (r1 :p-value)))
      "identical samples -> D=0, p=1 (the q<=0 guard)")
  (def r2 (tests/ks-test @[5 5 5] @[7 7 7] {}))
  (is (and (= 1 (r2 :statistic)) (r2 :ties)) "all-equal pools: D=1, ties flagged")
  (is (close 0.1 (r2 :p-value)) "D=1 exact p = 0.1 (R parity)")
  (def r3 (tests/ks-test @[1 2 2 4] @[1 2 3 3] {}))
  (is (and (close 0.25 (r3 :statistic)) (= 1 (r3 :p-value)))
      "mixed ties: D=0.25, p=1 (R parity)"))

(deftest ks-error-taxonomy
  # api-contract S3: the VALIDATION message shapes; the COMPOSITION
  # unpack-first error for the sort-bound paths
  (var caught nil)
  (try (tests/ks-test @[] @[] {}) ([e] (set caught e)))
  (is (and caught (string/find "needs at least one observation" caught))
      "empty arrays -> the VALIDATION shape")
  (set caught nil)
  (try (tests/ks-test @[1 2 math/nan] @[1 2] {}) ([e] (set caught e)))
  (is (and caught (string/find "needs finite observations" caught))
      "non-finite input -> the VALIDATION shape")
  (set caught nil)
  (try (tests/ks-test @[1 2] @[3 4] {:alternative :both}) ([e] (set caught e)))
  (is (and caught (string/find ":alternative must be" caught))
      "bad :alternative -> the VALIDATION shape")
  (set caught nil)
  (try (tests/ks-test @[1 2] "pnorm" {}) ([e] (set caught e)))
  (is (and caught (string/find "y must be an array, a keyword" caught))
      "unknown y kind -> the VALIDATION shape")
  (set caught nil)
  (try (tests/ks-test {:nrow 2 :ncol 1} @[1 2] {}) ([e] (set caught e)))
  (is (and caught (string/find "unpack packed handles first" caught))
      "packed handle into the sort-bound path -> the COMPOSITION shape"))

(deftest ks-two-sample-walk-parity
  # S1.1 rule 5: the chartered vscan walk == the pure walk - D and p
  # IDENTICAL (both strict left folds over the same ±1/n sequence)
  (def x (seq [i :range [0 30]] (math/sin (* i 0.3))))
  (def y (seq [i :range [0 25]] (+ 0.5 (math/cos (* i 0.21)))))
  (def rp (tests/ks-test x y {:impl :pure}))
  (def rn (tests/ks-test x y {:impl :num}))
  (is (= (rp :statistic) (rn :statistic))
      "vscan walk D == pure walk D (bitwise)")
  (is (= (rp :p-value) (rn :p-value)) "exact p identical both impls"))

#==------------------------------------------------------------------------==#
# km-logrank.qnt - FORWARD CONTRACT ACTIVATED (KM +
# two-group log-rank). THE pinned convention: at tied times DEATHS are processed
# FIRST - the factor's risk set INCLUDES same-time censorings
# (confirmed: survdiff.fit's order(strat, time, -status) sort +
# survfit's n.risk, R 4.5.0 / survival 3.8.6).
#==------------------------------------------------------------------------==#

(deftest km-risk-set-conservation
  # km-logrank.qnt invariants: S non-increasing in [0,1]; cohort
  # conservation; the deaths-before-censorings factor (R survfit
  # parity: n-risk(2) = 5 with the same-time censoring INCLUDED,
  # S(2) = 5/6 * 4/5 = 2/3)
  (def r (regsurv/survfit-km @[1.0 2.0 2.0 3.0 4.0 5.0] @[1 0 1 0 1 0] {}))
  (is (= 6 (r :n)) "cohort size")
  (def s (r :surv))
  (var s-prev 1.0)
  (var mono true)
  (each v s (when (> v s-prev) (set mono false)) (set s-prev v))
  (is mono "S non-increasing")
  (is (all |(and (<= 0 $) (<= $ 1)) s) "S in [0,1]")
  (is (= (length s) (length (r :time))) "curve points aligned")
  (is (= 5 (get (r :n-risk) 1)) "n-risk(2) = 5 - same-time censoring at risk")
  (is (close (/ 2 3) (get (r :surv) 1)) "S(2) = 2/3 (the factor with censor included)")
  (def r-n (regsurv/survfit-km @[1.0 2.0 2.0 3.0 4.0 5.0]
                               @[1 0 1 0 1 0] {:impl :num}))
  (is (deep= (r :surv) (r-n :surv)) "vscan :prod curve == pure loop (bitwise)"))

(deftest logrank-at-risk-decomposition
  # km-logrank.qnt invariants: the two-group decomposition conserves
  # (O total == E total == total deaths); chi/p present. Expected
  # values verified vs survival::survdiff 2026-10-04 (chisq
  # 0.79178082191780785, p 0.37356218672091324)
  (def r (regsurv/survdiff-logrank
           @[1 2 2 3 4 5 2 3 3 6] @[1 0 1 0 1 0 1 1 0 1]
           @[1 1 1 1 1 2 2 2 2 2] {}))
  (is (and (number? (r :statistic)) (number? (r :p-value)))
      "chi-square statistic + p-value present")
  (is (close 0.79178082191780785 (r :statistic)) "chisq vs R")
  (is (close 0.37356218672091324 (r :p-value)) "p vs R")
  (is (= (+ ((r :observed) 0) ((r :observed) 1))
         (+ ((r :expected) 0) ((r :expected) 1)))
      "O total == E total (the conservation: O1+O2 == d at every time)"))

(deftest km-error-taxonomy
  # api-contract S3: the VALIDATION shapes for the survival pair
  (var caught nil)
  (try (regsurv/survfit-km @[1 2] @[1 3] {}) ([e] (set caught e)))
  (is (and caught (string/find "status must be" caught))
      "bad status -> the VALIDATION shape")
  (set caught nil)
  (try (regsurv/survfit-km @[1 2] @[1] {}) ([e] (set caught e)))
  (is (and caught (string/find "equal length" caught))
      "length mismatch -> the VALIDATION shape")
  (set caught nil)
  (try (regsurv/survdiff-logrank @[1 2] @[1 1] @[1 1] {}) ([e] (set caught e)))
  (is (and caught (string/find "exactly 2 groups" caught))
      "1-group log-rank -> the VALIDATION shape")
  (set caught nil)
  (try (regsurv/survfit-km @[1 math/inf] @[1 0] {}) ([e] (set caught e)))
  (is (and caught (string/find "finite" caught))
      "non-finite time -> the VALIDATION shape")
  (set caught nil)
  (try (regsurv/survfit-km {:nrow 2 :ncol 1} @[1 0] {}) ([e] (set caught e)))
  (is (and caught (string/find "janet arrays" caught))
      "packed handle into the sort-bound path -> the COMPOSITION shape"))

(deftest dist-5b-family-contracts
  # p monotone + support bounds; q(p(x)) ~ id on a
  # grid; densities integrate ~1 (coarse sum); vector == scalar under
  # both impls. R-convention contracts.
  (def grid (seq [i :range [0 21]] (* 0.6 i)))
  (defn mono? [xs]
    (var ok true)
    (var prev (- math/inf))
    (each v xs (when (< v prev) (set ok false)) (set prev v))
    ok)
  # exp
  (is (mono? (map |(dist/pexp $ 0.7 {}) grid)) "pexp monotone")
  (is (and (zero? (dist/pexp (- 1) 0.7 {})) (= 1 (dist/pexp math/inf 0.7 {})))
      "pexp support bounds")
  (is (all |(< (math/abs (- (dist/pexp (dist/qexp $ 0.7 {}) 0.7 {}) $)) 1e-12)
           [0.01 0.2 0.5 0.8 0.99])
      "qexp/pexp round trip")
  (var mass-e 0)
  (each x (seq [i :range [1 4000]] (* 0.01 i))
    (+= mass-e (* 0.01 (dist/dexp x 0.7 {}))))
  (is (< (math/abs (- 1 mass-e)) 0.01) "dexp integrates to ~1 (coarse left-sum)")
  # weibull
  (is (mono? (map |(dist/pweibull $ 2.2 1.5 {}) grid)) "pweibull monotone")
  (is (all |(< (math/abs (- (dist/pweibull (dist/qweibull $ 2.2 1.5 {}) 2.2 1.5 {}) $)) 1e-12)
           [0.01 0.2 0.5 0.8 0.99])
      "qweibull/pweibull round trip")
  # lognormal
  (is (mono? (map |(dist/plnorm $ 0.1 0.5 {}) grid)) "plnorm monotone")
  (is (and (zero? (dist/plnorm 0 0.1 0.5 {})) (= 1 (dist/plnorm math/inf 0.1 0.5 {})))
      "plnorm support bounds")
  (is (all |(< (math/abs (- (dist/plnorm (dist/qlnorm $ 0.1 0.5 {}) 0.1 0.5 {}) $)) 1e-9)
           [0.01 0.2 0.5 0.8 0.99])
      "qlnorm/plnorm round trip")
  # nbinom (REAL size supported)
  (is (mono? (map |(dist/pnbinom $ 5 0.3 {}) (range 0 40))) "pnbinom monotone")
  (is (and (zero? (dist/pnbinom (- 1) 5 0.3 {})) (= 1 (dist/pnbinom math/inf 5 0.3 {})))
      "pnbinom support bounds")
  (is (= (dist/qnbinom (dist/pnbinom 7 5 0.3 {}) 5 0.3 {}) 7) "qnbinom(pnbinom) = k")
  (is (close 1 (dist/pnbinom 1000 0.5 0.3 {})) "pnbinom far support -> 1")
  # vector == scalar, both impls (the both-impl rule)
  (def vx @[0.5 1.2 3.0])
  (each o [{:impl :pure} {:impl :num}]
    (is (deep= (map |(dist/dexp $ 0.7 o) vx) (dist/dexp vx 0.7 o))
        "dexp vector == scalar (both impls)")
    (is (deep= (map |(dist/dlnorm $ 0.1 0.5 o) vx) (dist/dlnorm vx 0.1 0.5 o))
        "dlnorm vector == scalar (both impls)"))
  # the x==0 && shape<1 edge -> +inf (R verbatim; suite-only - the
  # fixture dialect is finite-only)
  (is (= math/inf (dist/dweibull 0 0.5 1 {})) "dweibull(0, shape<1) = +inf"))

(deftest cdf-q-pole-conformance
  # CDF poles + q-endpoints match R
  # (verified live vs R 4.5.0, 2026-10-08). Scalar rows run both
  # impls (the scalar policy: identical pure formulas).
  (each o [{:impl :pure} {}]
    # CDF poles: +inf -> 1, negative/out-of-support -> 0 (no throw)
    (is (= 1 (dist/pchisq math/inf 5 o)) "pchisq(+inf) = 1")
    (is (zero? (dist/pchisq (- math/inf) 5 o)) "pchisq(-inf) = 0")
    (is (= 1 (dist/pgamma math/inf 2 1 o)) "pgamma(+inf) = 1")
    (is (zero? (dist/pgamma (- 1) 2 1 o)) "pgamma(-1) = 0")
    (is (= 1 (dist/ppois math/inf 3 o)) "ppois(+inf) = 1")
    (is (zero? (dist/ppois (- 1) 3 o)) "ppois(-1) = 0")
    (is (= 1 (dist/pbinom math/inf 10 0.5 o)) "pbinom(+inf) = 1")
    (is (zero? (dist/pbinom (- 1) 10 0.5 o)) "pbinom(-1) = 0")
    (is (= 1 (dist/pbinom 12 10 0.5 o)) "pbinom k>n saturates at 1 (R)")
    (is (= 1 (dist/pbeta math/inf 2 3 o)) "pbeta(+inf) = 1")
    (is (zero? (dist/pbeta (- 1) 2 3 o)) "pbeta(-1) = 0")
    (is (= 1 (dist/pbeta 1 2 3 o)) "pbeta(1) = 1 (no throw)")
    # q-endpoints: no bracket saturation on legal p
    (is (= (- math/inf) (dist/qnorm 0 0 1 o)) "qnorm(0) = -inf")
    (is (= math/inf (dist/qnorm 1 0 1 o)) "qnorm(1) = +inf")
    (is (= (- math/inf) (dist/qt 0 3 o)) "qt(0) = -inf")
    (is (= math/inf (dist/qt 1 3 o)) "qt(1) = +inf")
    (is (zero? (dist/qchisq 0 5 o)) "qchisq(0) = 0")
    (is (= math/inf (dist/qchisq 1 5 o)) "qchisq(1) = +inf")
    (is (zero? (dist/qgamma 0 2 1 o)) "qgamma(0) = 0")
    (is (= math/inf (dist/qgamma 1 2 1 o)) "qgamma(1) = +inf")
    (is (zero? (dist/qbeta 0 2 3 o)) "qbeta(0) = 0")
    (is (= 1 (dist/qbeta 1 2 3 o)) "qbeta(1) = 1")
    # pure vector lane: vapply composes the same scalar poles
    (when (= o {:impl :pure})
      (is (deep= @[1 0] (dist/pchisq @[math/inf (- math/inf)] 5 o))
          "pchisq vector poles (pure)"))))

(deftest dbinom-dpois-support-conformance
  # the R support
  # convention - out-of-support k -> 0, non-integer k -> NaN (R: NaN +
  # warning; janet has no warning channel). The old throws withdrawn.
  (each o [{:impl :pure} {:impl :num}]
    (is (zero? (dist/dbinom (- 1) 10 0.5 o)) "dbinom(-1) = 0 (R)")
    (is (zero? (dist/dbinom 12 10 0.5 o)) "dbinom k>n = 0 (R)")
    (is (= (- math/inf) (dist/dbinom (- 1) 10 0.5 (merge {:log true} o)))
        "dbinom(-1, log) = -inf")
    #  reversal: beyond-tolerance non-integer -> 0 (R's
    # own value + warning; the NaN was a mis-decode, reversed)
    (is (zero? (dist/dbinom 2.5 10 0.5 o))
        "dbinom non-integer = 0 (R convention, reversed)")
    (is (zero? (dist/dpois (- 1) 3 o)) "dpois(-1) = 0 (R)")
    (is (zero? (dist/dpois 1.5 3 o))
        "dpois non-integer = 0 (R convention, reversed)")
    (is (close 0.205078125 (dist/dbinom 6 10 0.5 o))
        "dbinom(6,10,.5) unchanged (R spot)")
    (is (close 0.22404180765538775 (dist/dpois 3 3 o))
        "dpois(3,3) unchanged (R spot)"))
  # kernel vector lane: good arrays stay kernel; bad-element arrays
  # fall to the elementwise pure path with R values elementwise
  # (parity bar 1e-9 - the kernel lane's lgamma association diff
  # ~1e-15 is inside the ABI bar, not bitwise)
  (def gv (dist/dbinom @[6 7] 10 0.5 {:impl :num}))
  (is (and (< (math/abs (- 0.205078125 (gv 0))) 1e-9)
           (< (math/abs (- 0.1171875 (gv 1))) 1e-9))
      "dbinom good vector (kernel lane, ABI bar)")
  (def mix (dist/dbinom @[6 2.5 (- 1)] 10 0.5 {:impl :num}))
  (is (and (< (math/abs (- 0.205078125 (mix 0))) 1e-9)
           (zero? (mix 1))
           (zero? (mix 2)))
      "dbinom mixed vector: elementwise R values (fallback lane)"))

(deftest density-cdf-pole-rows
  # every
  # value verified live vs R 4.5.0 (2026-10-09)
  (each o [{:impl :pure} {}]
    # remaining density/CDF poles -> R's 0
    (is (zero? (dist/pf (- math/inf) 3 10 o)) "pf(-inf) = 0 (R)")
    (is (zero? (dist/pf (- 1) 3 10 o)) "pf(-1) = 0 (R)")
    (is (zero? (dist/df-f math/inf 3 10 o)) "df-f(inf) = 0 (R)")
    (is (zero? (dist/df-f (- math/inf) 3 10 o)) "df-f(-inf) = 0 (R)")
    (is (zero? (dist/dbeta math/inf 2 3 o)) "dbeta(inf) = 0 (R)")
    (is (zero? (dist/dbeta (- math/inf) 2 3 o)) "dbeta(-inf) = 0 (R)")
    (is (zero? (dist/dgamma math/inf 2 1 o)) "dgamma(inf) = 0 (R)")
    (is (zero? (dist/dgamma (- math/inf) 2 1 o)) "dgamma(-inf) = 0 (R)")
    (is (zero? (dist/dchisq math/inf 3 o)) "dchisq(inf) = 0 (R)")
    (is (zero? (dist/dchisq (- math/inf) 3 o)) "dchisq(-inf) = 0 (R)")
    (is (zero? (dist/dexp (- 1) 1 o)) "dexp(-1) = 0 (R)")
    (is (zero? (dist/dpois math/inf 2 o)) "dpois(inf) = 0 (R)")
    # dpois lambda=0 point mass
    (is (= 1 (dist/dpois 0 0 o)) "dpois(0,0) = 1 (R)")
    (is (zero? (dist/dpois 3 0 o)) "dpois(k>0,0) = 0 (R)")
    (is (= 0 (dist/dpois 0 0 (merge {:log true} o))) "dpois(0,0,log) = 0")
    # qf endpoints
    (is (zero? (dist/qf 0 3 10 o)) "qf(0) = 0 (R)")
    (is (= math/inf (dist/qf 1 3 10 o)) "qf(1) = +inf (R)")
    # floor(q + 1e-7) coercion (nmath), both lanes
    (is (close 0.3827827864 (dist/pbinom 2.5 10 0.3 o))
        "pbinom(2.5) = P(<=2) (R floor)")
    (is (close 0.6496107184 (dist/pbinom 2.9999999 10 0.3 o))
        "pbinom(2.9999999) = P(<=3) (floor + 1e-7)")
    (is (close 0.6766764161830634 (dist/ppois 2.5 2 o))
        "ppois(2.5) = P(<=2) (R floor)")
    # dnbinom support convention, family-wide
    (is (zero? (dist/dnbinom math/inf 5 0.5 o)) "dnbinom(inf) = 0 (R)")
    (is (zero? (dist/dnbinom (- 1) 5 0.5 o)) "dnbinom(-1) = 0 (R)")
    (is (zero? (dist/dnbinom 2.5 5 0.5 o))
        "dnbinom non-integer = 0 (R convention, family-wide)"))
  # kernel-lane vector rows: non-integer/inf arrays fall to the
  # elementwise floor path; integer arrays stay kernel
  (def pf25 (dist/ppois @[2.5 2] 2 {:impl :num}))
  (is (and (<= (math/abs (- 0.6766764161830634 (get pf25 0))) 1e-9)
           (<= (math/abs (- 0.6766764161830634 (get pf25 1))) 1e-9))
      "ppois vector floors elementwise (num lane)")
  (def pv (dist/ppois @[2 3] 2 {:impl :num}))
  (is (and (<= (math/abs (- 0.6766764161830634 (get pv 0))) 1e-9)
           (<= (math/abs (- 0.857123460498547 (get pv 1))) 1e-9))
      "ppois integer vector stays kernel (num lane)")
  (is (deep= @[1 0 0] (dist/dpois @[0 1 2] 0 {:impl :num}))
      "dpois lambda=0 vector point mass (num lane, fallback)")
  # the m-slot trap guards actionably
  (var caught nil)
  (try (descriptive/variance @[1 2 3] {:impl :num}) ([e] (set caught e)))
  (is (and caught (string/find "precomputed mean" caught))
      "variance opts-into-m trap -> the actionable guard ")
  (set caught nil)
  (try (descriptive/skewness @[1 2 3 4] {:impl :num}) ([e] (set caught e)))
  (is (and caught (string/find "precomputed moments" caught))
      "skewness m/s slot trap -> the actionable guard "))

(deftest dlnorm-tolerance-rows
  # every value verified live
  # vs R 4.5.0 (2026-10-09); S5's warning rows are subprocess probes
  # (the in-process warning fires once per process by contract)
  (each o [{:impl :pure} {}]
    # dlnorm pole (the survivor - its doc comment was right)
    (is (zero? (dist/dlnorm -1 0 1 o)) "dlnorm(-1) = 0 (R)")
    (is (zero? (dist/dlnorm 0 0 1 o)) "dlnorm(0) = 0 (R)")
    (is (= (- math/inf) (dist/dlnorm 0 0 1 (merge {:log true} o)))
        "dlnorm(0, log) = -inf (R)")
    (is (close 0.62749607711592397 (dist/dlnorm 0.5 0 1 o))
        "dlnorm(0.5) unchanged (R spot)")
    # R's relative 1e-9 non-integer tolerance (boundary-measured:
    # x=2 -> 2e-9, x=7 -> 7e-9, x=200 -> 2e-7)
    (is (close 0.3125 (dist/dbinom (+ 2 1.9e-9) 5 0.5 o))
        "dbinom(2+1.9e-9) = dbinom(2) (R tolerance)")
    (is (zero? (dist/dbinom (+ 2 2e-9) 5 0.5 o))
        "dbinom(2+2e-9) = 0 (beyond tolerance, = R's 0)")
    (is (close 0.18044704431548301 (dist/dpois (+ 3 2e-9) 2 o))
        "dpois(3+2e-9) = dpois(3) (R tolerance)")
    (is (close 0.080566406250000006 (dist/dnbinom (+ 7 6.9e-9) 5 0.5 o))
        "dnbinom(7+6.9e-9) = dnbinom(7) (R tolerance)")
    (is (close 6.22301527786114e-61 (dist/dbinom (+ 200 1.9e-7) 200 0.5 o))
        "dbinom(200+1.9e-7) = dbinom(200) (relative at scale)")
    # vector lane: tolerance-accepted elements floor elementwise
    (def tv (dist/dbinom @[(+ 2 1.9e-9) 2] 5 0.5 (merge {:impl :num} o)))
    (is (and (< (math/abs (- 0.3125 (get tv 0))) 1e-9)
             (< (math/abs (- 0.3125 (get tv 1))) 1e-9))
        "dbinom tolerance vector floors elementwise (num lane, ABI bar)"))
  # guard type coverage - every non-number m/s gets the D4 message
  (each bad [@{:impl :num} [1] "m" :kw]
    (var caught nil)
    (try (descriptive/variance @[1 2 3] bad {}) ([e] (set caught e)))
    (is (and caught (string/find "precomputed mean" caught))
        (string "variance m-slot " (type bad) " -> D4 guard"))
    (set caught nil)
    (try (descriptive/kurtosis @[1 2 3 4 5] nil bad {}) ([e] (set caught e)))
    (is (and caught (string/find "precomputed moments" caught))
        (string "kurtosis s-slot " (type bad) " -> D4 guard")))
  # nil m and numeric m still work (the fast path + the auto path)
  (is (= 1 (descriptive/variance @[1 2 3] 2 {})) "numeric m unchanged")
  (is (= 1 (descriptive/variance @[1 2 3] nil {})) "nil m auto-computes"))

(deftest dbinom-tolerance-guard-rows
  # every value
  # verified live vs R 4.5.0 from BOTH sides (2026-10-09)
  (each o [{:impl :pure} {}]
    # BELOW-side tolerance band = the integer mass (R rounds)
    (is (close 0.0439453125 (dist/dbinom (- 2 2e-9) 10 0.5 o))
        "dbinom(2-2e-9) = dbinom(2) (R rounds - the below side)")
    (is (close 0.0439453125 (dist/dbinom (+ 2 1.9e-9) 10 0.5 o))
        "dbinom(2+1.9e-9) = dbinom(2) (the above side stays)")
    (is (close 6.22301527786114e-61 (dist/dbinom (- 200 1.9e-7) 200 0.5 o))
        "dbinom(200-1.9e-7) = dbinom(200) (below at scale)")
    (is (close 0.180447044315483 (dist/dpois (- 3 2.9e-9) 2 o))
        "dpois(3-2.9e-9) = dpois(3) (below)")
    (is (close 0.08056640625 (dist/dnbinom (- 7 6.9e-9) 5 0.5 o))
        "dnbinom(7-6.9e-9) = dnbinom(7) (below)")
    # reversal: beyond-tolerance non-integer -> 0 (R's value)
    (is (zero? (dist/dbinom 2.5 10 0.5 o)) "dbinom(2.5) = 0 (R; the NaN reversed)")
    (is (zero? (dist/dbinom (+ 2 2e-9) 10 0.5 o)) "dbinom(2+2e-9) = 0 (R)")
    (is (zero? (dist/dpois 2.5 2 o)) "dpois(2.5) = 0 (R)")
    (is (zero? (dist/dnbinom 2.5 5 0.5 o)) "dnbinom(2.5) = 0 (R)")
    (is (= (- math/inf) (dist/dbinom 2.5 10 0.5 (merge {:log true} o)))
        "dbinom(2.5, log) = -inf (R)")
    # the tiny-negative p zone -> 0 (R; not P(X<=0))
    (each neg [-1e-10 -1e-8 (- 0.5e-7) -1.5e-7 -1]
      (is (zero? (dist/pbinom neg 10 0.5 o))
          (string "pbinom(" neg ") = 0 (tiny-negative zone, R)"))
      (is (zero? (dist/ppois neg 3 o))
          (string "ppois(" neg ") = 0 (tiny-negative zone, R)")))
    # the floor coercions stay correct
    (is (close 0.3827827864 (dist/pbinom 2.5 10 0.3 o)) "pbinom(2.5) floors (stays)")
    # q-lane poles NaN elementwise, both out-of-domain sides
    (each bad [-1 1.5]
      (is (not= (dist/qnorm bad 0 1 o) (dist/qnorm bad 0 1 o))
          (string "qnorm(" bad ") = NaN (R)")))
    (each [qf a2 a3] [[dist/qbinom 10 0.5] [dist/qweibull 2 1]
                       [dist/qlnorm 0 1] [dist/qnbinom 5 0.5]
                       [dist/qbeta 2 3] [dist/qgamma 2 1]]
      (is (not= (qf 1.5 a2 a3 o) (qf 1.5 a2 a3 o))
          "q family out-of-domain p = NaN (R; no vector-abort throws)")
      (is (not= (qf (- 1) a2 a3 o) (qf (- 1) a2 a3 o))
          "q family negative p = NaN (R)"))
    (each [qf a2] [[dist/qpois 3] [dist/qexp 1]]
      (is (not= (qf 1.5 a2 o) (qf 1.5 a2 o))
          "q family out-of-domain p = NaN (R; no vector-abort throws)")
      (is (not= (qf (- 1) a2 o) (qf (- 1) a2 o))
          "q family negative p = NaN (R)"))
    (def qw (dist/qweibull @[0.25 1 1.5 (- 0.1)] 2 1 o))
    (is (and (close 0.53636002130265183 (get qw 0))
             (= math/inf (get qw 1))
             (not= (get qw 2) (get qw 2))
             (not= (get qw 3) (get qw 3)))
        "qweibull vector: elementwise NaN, the rest computes (R)")
    # mad's positional center honored
    (is (close 1.4826 (descriptive/mad @[1 2 3 4] 2 o))
        "mad(xs, 2) honors the positional center (R's mad(x, 2))"))
  # guards: the D4-class messages
  (var caught nil)
  (try (tests/wilcoxon-sr @[1 2 3] @{:impl :num}) ([e] (set caught e)))
  (is (and caught (string/find "null mean" caught))
      "wilcoxon-sr mu-slot trap -> actionable guard")
  (set caught nil)
  (try (tests/poisson-test 5 "junk") ([e] (set caught e)))
  (is (and caught (string/find "time base" caught))
      "poisson-test t-slot trap -> actionable guard")
  # handles accepted by dlnorm/plnorm kernel lanes
  (each o [{:impl :num}]
    (def h (num/stage @[0.5 1.0] (num/kernels)))
    (def hd (dist/dlnorm h 0 1 o))
    (is (and (<= (math/abs (- 0.6274960771159244 (get hd 0))) 1e-9)
             (<= (math/abs (- 0.3989422804014327 (get hd 1))) 1e-9))
        "dlnorm handle: kernel lane consumes it (gate)")))

(deftest q-domain-rows
  # NaN-p and out-of-domain p
  # return NaN elementwise across ALL q families (R's ISNAN check;
  # janet comparisons are false on NaN so the old guards let it
  # through returning solver garbage). Interiors intact.
  (def nan (/ 0.0 0.0))
  (each o [{:impl :pure} {}]
    # one NaN-p row per family (the twelve + qtukey)
    (each [qf args] [[dist/qnorm [0 1]] [dist/qbinom [10 0.5]]
                     [dist/qpois [3]] [dist/qbeta [2 3]]
                     [dist/qgamma [2 1]] [dist/qexp [1]]
                     [dist/qweibull [2 1]] [dist/qlnorm [0 1]]
                     [dist/qnbinom [5 0.5]] [dist/qt [5]]
                     [dist/qchisq [3]] [dist/qf [4 7]]]
      (is (not= (qf nan ;args o) (qf nan ;args o))
          "q family NaN p = NaN (no solver garbage)"))
    (is (not= (dist/qtukey nan 4 12 o) (dist/qtukey nan 4 12 o))
        "qtukey NaN p = NaN")
    # out-of-domain finite p: the three converted  + qtukey 
    (is (not= (dist/qt 1.5 5 o) (dist/qt 1.5 5 o)) "qt(1.5) = NaN (was bracket garbage)")
    (is (not= (dist/qchisq (- 0.1) 3 o) (dist/qchisq (- 0.1) 3 o))
        "qchisq(-0.1) = NaN (was a whole-vector throw)")
    (is (not= (dist/qf 1.5 4 7 o) (dist/qf 1.5 4 7 o))
        "qf(1.5) = NaN (was bracket garbage)")
    (is (not= (dist/qtukey (- 0.1) 4 12 o) (dist/qtukey (- 0.1) 4 12 o))
        "qtukey(-0.1) = NaN (was a throw)")
    # elementwise vectors: the NaN slot NaNs, the rest computes
    # (deep= cannot compare NaNs - NaN != NaN - so per-element checks)
    (def qv (dist/qnorm @[0.5 nan 0.9] 0 1 o))
    (is (and (zero? (get qv 0))
             (not= (get qv 1) (get qv 1))
             (close 1.2815515655446009 (get qv 2)))
        "qnorm mixed vector: NaN slot NaN, rest computes (R)")
    (def cv (dist/qchisq @[0.05 1.5 (- 0.2)] 3 o))
    (is (and (close 0.351846317749272 (get cv 0))
             (not= (get cv 1) (get cv 1))
             (not= (get cv 2) (get cv 2)))
        "qchisq mixed vector: elementwise, no abort (R)")
    # interiors intact (R spots)
    (is (close 1.959963984540054 (dist/qnorm 0.975 0 1 o)) "qnorm(.975) intact")
    (is (close 2.57058183563632 (dist/qt 0.975 5 o)) "qt(.975,5) intact (R)")
    (is (close 7.81472790325119 (dist/qchisq 0.95 3 o)) "qchisq(.95,3) intact (R)")
    (is (close 4.12031172689764 (dist/qf 0.95 4 7 o)) "qf(.95,4,7) intact (R)")
    (is (close 4.19866023125160 (dist/qtukey 0.95 4 12 o)) "qtukey(.95) intact (R)")))

(deftest far-tail-rows
  # the far-tail rows BOTH SIDES
  # per solver family - the missing  instrumentation; every
  # value verified live vs R 4.5.0 (2026-10-10); the 60-cell crash
  # grid (shapes 5-50 x p band x rates) green, 0 fails
  (each o [{:impl :pure} {}]
    # qgamma: the bar cells (were THROWS) + the band + both
    # sides (near-0 mirrors) + the pinned  cells
    (is (close 49.227431945838042 (dist/qgamma (- 1 1e-16) 5 1 o))
        "qgamma(1-1e-16,5,1) = R (was a THROW)")
    (is (close 60.2931901883112 (dist/qgamma (- 1 2e-16) 10 1 o))
        "qgamma(1-2e-16,10,1) = R (was a THROW)")
    (is (close 79.100907623178358 (dist/qgamma (- 1 5e-16) 20 1 o))
        "qgamma(1-5e-16,20,1) = R")
    (is (close 29.948824800835411 (dist/qgamma (- 1 1e-14) 0.5 1 o))
        "qgamma(1-1e-14,.5) = R (the bar cell)")
    (is (close 34.381626105834208 (dist/qgamma (- 1 1e-16) 0.5 1 o))
        "qgamma(1-1e-16,.5) = R (the bar cell)")
    (is (close 75.287098905440644 (dist/qgamma (- 1 1e-14) 20 1 o))
        "qgamma(1-1e-14,20) = R (the bar cell)")
    (is (close 31.099896028389072 (dist/qgamma (- 1 1e-12) 2 1 o))
        "qgamma(1-1e-12,2,1) = R EXACT (the pinned cell)")
    (is (close 39.235847839113838 (dist/qgamma (- 1 1e-12) 5 1 o))
        "qgamma(1-1e-12,5,1) = R")
    (is (close 1.4142135690397595e-08 (dist/qgamma 1e-16 2 1 o))
        "qgamma(1e-16,2,1) = R (the near-0 side)")
    (is (close 4.4721360216662435e-08 (dist/qgamma 1e-15 2 1 o))
        "qgamma(1e-15,2,1) = R (the near-0 side)")
    # qchisq: the bar cells + the contradiction cell
    (is (close 98.454863891676084 (dist/qchisq (- 1 1e-16) 10 o))
        "qchisq(1-1e-16,10) = R (was a THROW; S8: the constant was a digit transposition)")
    (is (close 162.40041295406323 (dist/qchisq (- 1 1e-16) 40 o))
        "qchisq(1-1e-16,40) = R (was a THROW)")
    (is (close 77.396315490620879 (dist/qchisq (- 1 1e-16) 3 o))
        "qchisq(1-1e-16,3) = R (the 3.4e-3 cell)")
    (is (close 58.91980066488555 (dist/qchisq (- 1 1e-12) 3 o))
        "qchisq(1-1e-12,3) = R (the contradiction cell, was 1.9e-6 off)")
    (is (close 1.1223330973050010e-05 (dist/qchisq 1e-8 3 o))
        "qchisq(1e-8,3) = R (near-0)")
    (is (close 1.5707963267948893e-16 (dist/qchisq 1e-8 1 o))
        "qchisq(1e-8,1) = R (near-0)")
    # interiors unchanged
    (is (close 4.743864518390577 (dist/qgamma 0.95 2 1 o))
        "qgamma(.95,2,1) interior intact")
    (is (close 7.8147279032511836 (dist/qchisq 0.95 3 o))
        "qchisq(.95,3) interior intact"))
  # the mixed vector: elementwise, no abort 
  (def mv (dist/qchisq @[0.05 (- 1 1e-16) 0.5] 10 {}))
  (is (and (close 3.94029913611906 (get mv 0))
           (close 98.4548638916761 (get mv 1))
           (close 9.34181776559196 (get mv 2)))
      "qchisq mixed near-1 vector: elementwise, no abort "))

(deftest qgamma-complement-rows
  # the
  # upper lane (:lower-tail false) answers the RAW upper p, not the
  # materialized double-rounded complement - rows keyed to the input
  # CLASS (upper-lane tiny-p ladder per complement family, shape
  # ladder, both qf fudge orders, the raw-vs-complement control pair,
  # the swallow boundary), every value verified live vs R 4.5.0
  # (rprobe.R c1_ cells, 2026-10-13). Bar: rel <= 1e-12.
  (each o [{:impl :pure} {}]
    # the shape ladder at the swallow boundary p (was rel 2.3e-3)
    (is (close 49.341017193401548 (dist/qgamma 1e-16 5 1 {:lower-tail false}))
        "qgamma(1e-16,5,1,up) = R (was the complement-rounded 49.2274)")
    (is (and (< (math/abs (dist/qgamma 5e-17 5 1 {:lower-tail false})) math/inf)
            (close 50.093435862224609 (dist/qgamma 5e-17 5 1 {:lower-tail false})))
        "qgamma(5e-17,5,1,up) finite = R (was inf, the 1-p==1.0 swallow)")
    (is (close 34.484730479258282 (dist/qgamma 1e-16 0.5 1 {:lower-tail false}))
        "qgamma(1e-16,.5,1,up) = R")
    (is (close 40.568709189188226 (dist/qgamma 1e-16 2 1 {:lower-tail false}))
        "qgamma(1e-16,2,1,up) = R")
    (is (close 81.336032763908946 (dist/qgamma 1e-16 20 1 {:lower-tail false}))
        "qgamma(1e-16,20,1,up) = R")
    (is (close 131.98426500768386 (dist/qgamma 1e-16 50 1 {:lower-tail false}))
        "qgamma(1e-16,50,1,up) = R")
    # the p ladder (bias grows as p shrinks - 6.25e-7 at 1e-12)
    (is (close 39.235823281422618 (dist/qgamma 1e-12 5 1 {:lower-tail false}))
        "qgamma(1e-12,5,1,up) = R (was the complement-rounded 39.23585)")
    (is (close 44.315762597493901 (dist/qgamma 1e-14 5 1 {:lower-tail false}))
        "qgamma(1e-14,5,1,up) = the TRUE quantile (DEVIATION, cells-carrying:
         R prints 44.315762600777347 but R's own pgamma there is off by
         dlog -2.99e-9 (materialized-p noise); ours has dlog 0
         BY R's pgamma - verified 2026-10-13; same class at 1e-13/5e-13
         (R dlog -2.2e-9/-3.9e-11); at 1e-12 and below R converges)")
    (is (close 23.431523423392193 (dist/qgamma 1e-6 5 1 {:lower-tail false}))
        "qgamma(1e-6,5,1,up) = R")
    # qchisq inherits (composes qgamma)
    (is (close 84.411683100633269 (dist/qchisq 1e-16 5 {:lower-tail false}))
        "qchisq(1e-16,5,up) = R (was rel 2.57e-3)")
    (is (close 98.682034386803096 (dist/qchisq 1e-16 10 {:lower-tail false}))
        "qchisq(1e-16,10,up) = R")
    (is (close 65.238636213368267 (dist/qchisq 1e-12 5 {:lower-tail false}))
        "qchisq(1e-12,5,up) = R")
    (is (close 263.96853001536772 (dist/qchisq 1e-16 100 {:lower-tail false}))
        "qchisq(1e-16,100,up) = R")
    # qf fudge lane, BOTH df orders at the far tail
    (is (close 16.882336620126654 (dist/qf 1e-16 5 400001 {:lower-tail false}))
        "qf(1e-16,5,400001,up) = R (was rel 2.57e-3)")
    (is (close 3884270.3117717253 (dist/qf 1e-16 400001 5 {:lower-tail false}))
        "qf(1e-16,400001,5,up) = R (the other fudge order)")
    (is (close 13.047727242673654 (dist/qf 1e-12 5 400001 {:lower-tail false}))
        "qf(1e-12,5,400001,up) = R")
    (is (close 1.0186144513439341 (dist/qf 1e-16 400001 400000 {:lower-tail false}))
        "qf(1e-16,400001,400000,up) = R (the symmetric fudge cell)")
    # the control pair: the upper lane must NOT reproduce the lower
    # lane at the materialized complement (R itself differs - rel
    # 6.25e-7 at 1e-12; the pre-fix code returned EXACTLY this)
    (is (not (close 39.235847839113838 (dist/qgamma 1e-12 5 1 {:lower-tail false})))
        "qgamma(1e-12,5,1,up) != qgamma(1-1e-12,5,1) (the control pair)")
    (is (not (close 49.227431945838042 (dist/qgamma 1e-16 5 1 {:lower-tail false})))
        "qgamma(1e-16,5,1,up) != qgamma(1-1e-16,5,1) (the control pair)")
    # upper-lane interiors hold (large raw upper p)
    (is (close 1.6234863901184209 (dist/qgamma 0.975 5 1 {:lower-tail false}))
        "qgamma(.975,5,1,up) = R (interior control)")
    (is (close 5.5225943453085504 (dist/qf 0.025 4 7 {:lower-tail false}))
        "qf(.025,4,7,up) = R (interior control)")
    # the swallow boundary: 1e-17 finite (the old boundary fired at
    # ~1.1e-16); the subnormal tail below the cf's linear floor stays
    # a documented boundary (see api-contract (d)). finite? guards:
    # close() degenerates on an inf actual ((<= inf inf) via the rel
    # clause) - the rows must pin FINITENESS explicitly
    (is (and (< (math/abs (dist/qgamma 1e-17 5 1 {:lower-tail false})) math/inf)
            (close 51.836936947384352 (dist/qgamma 1e-17 5 1 {:lower-tail false})))
        "qgamma(1e-17,5,1,up) finite = R (16 orders below the old boundary)")
    (is (and (< (math/abs (dist/qgamma 1e-30 5 1 {:lower-tail false})) math/inf)
            (close 83.654670315294041 (dist/qgamma 1e-30 5 1 {:lower-tail false})))
        "qgamma(1e-30,5,1,up) finite = R")
    (is (and (< (math/abs (dist/qgamma 1e-300 5 1 {:lower-tail false})) math/inf)
            (close 713.88597806494431 (dist/qgamma 1e-300 5 1 {:lower-tail false})))
        "qgamma(1e-300,5,1,up) finite = R (R-exact through 1e-320)")))
  # the pinned  lower-lane cells must hold byte-identically
  # (the fix touches only the upper lane; guarded by far-tail-rows)

(deftest qgamma-median-rows
  # 's
  # log-scale Newton lower branch degraded the p=0.5 median at large
  # shape (post rel 7.01e-5..3.26e-4 where pre-CR held 7.79e-13..
  # 7.83e-8); the gate holds the log-scale route to the far-tail
  # domain it was ported for (the 1e-12..1e-16 band),
  # the median rides the pre-CR linear branch. Bars = the measured
  # pre-CR (012ac16) precision class per shape (the lgamma/series
  # noise grows with shape - pre-existing, S2's family); values
  # verified live vs R 4.5.0 (rprobe.R c3_ cells, 2026-10-13).
  (defn rel [x r] (math/abs (/ (- x r) r)))
  (each o [{:impl :pure} {}]
    # the recorded pre-fix values as the bar
    (is (< (rel (dist/qgamma 0.5 2e5 0.5 o) 399999.33333353087) 1e-12)
        "qgamma(.5,2e5,.5) at pre-CR precision (was rel 7.01e-5)")
    (is (< (rel (dist/qgamma 0.5 5e5 0.5 o) 999999.3333334123) 5e-10)
        "qgamma(.5,5e5,.5) at pre-CR precision (was rel 2.55e-4)")
    (is (< (rel (dist/qgamma 0.5 1e6 0.5 o) 1999999.3333333728) 1e-7)
        "qgamma(.5,1e6,.5) at pre-CR precision (was rel 3.26e-4)")
    (is (< (rel (dist/qchisq 0.5 4e5 o) 399999.33333353087) 1e-12)
        "qchisq(.5,4e5) at pre-CR precision (was rel 7.01e-5)")
    (is (< (rel (dist/qchisq 0.5 1e6 o) 999999.3333334123) 5e-10)
        "qchisq(.5,1e6) at pre-CR precision (was rel 2.55e-4)")
    # the p=0.5 +/- 1ulp ladder (the border sits exactly at 0.5)
    (is (< (rel (dist/qgamma (+ 0.5 2e-16) 2e5 0.5 o) 399999.33333353087) 1e-12)
        "qgamma(.5+1ulp,2e5,.5) at pre-CR precision")
    (is (< (rel (dist/qgamma (- 0.5 2e-16) 2e5 0.5 o) 399999.33333353087) 1e-12)
        "qgamma(.5-1ulp,2e5,.5) at pre-CR precision")
    (is (< (rel (dist/qgamma 0.4999999999 2e5 0.5 o) 399999.33333330665) 1e-11)
        "qgamma(.5-1e-10,2e5,.5) at pre-CR precision")
    (is (< (rel (dist/qgamma 0.5000000001 2e5 0.5 o) 399999.33333375509) 1e-11)
        "qgamma(.5+1e-10,2e5,.5) at pre-CR precision")
    (is (< (rel (dist/qgamma 0.5000000001 1e6 0.5 o) 1999999.3333338741) 1e-7)
        "qgamma(.5+1e-10,1e6,.5) at pre-CR precision")
    # the 0.55 control (upper branch, never regressed - the gate must
    # not move it)
    (is (< (rel (dist/qgamma 0.55 2e5 0.5 o) 400111.73867717071) 5e-12)
        "qgamma(.55,2e5,.5) control (upper branch)")
    (is (< (rel (dist/qgamma 0.55 1e6 0.5 o) 2000250.6665055074) 5e-12)
        "qgamma(.55,1e6,.5) control (upper branch)")
    # small-shape median controls (the linear branch's home domain)
    (is (close 9.3418177655919692 (dist/qgamma 0.5 5 0.5 o))
        "qgamma(.5,5,.5) interior control")
    (is (close 3.368600385977321 (dist/qgamma 0.25 5 1 o))
        "qgamma(.25,5,1) interior control")
    (is (< (rel (dist/qgamma 0.5 2e4 0.5 o) 39999.333335308678) 1e-12)
        "qgamma(.5,2e4,.5) at pre-CR precision")
    # the qf 4e5+1 fudge cells, both orders (rel 1.12e-6
    # at df1>df2;  improved 1000x vs pre, residual recorded)
    (is (< (rel (dist/qf 0.95 400001 4e5 o) 1.0036887163794868) 2e-6)
        "qf(.95,400001,4e5) <= the recorded 1.12e-6 class (boundary)")
    (is (< (rel (dist/qf 0.975 400001 400000 o) 1.0043971250813357) 1e-6)
        "qf(.975,400001,400000) <= the recorded 5.16e-7 class")
    (is (< (rel (dist/qf 0.95 4e5 400001 o) 1.0036808449008021) 1e-12)
        "qf(.95,4e5,400001) the other order (R-precision)")
    (is (< (rel (dist/qf 0.5 4e5 400001 o) 0.99999833333382715) 1e-6)
        "qf(.5,4e5,400001) the p=0.5 fudge cell (was rel 7.01e-5)")))

(deftest qgamma-huge-shape-rows
  # the
  # huge-shape family. (a) the qgamma final clamp no longer launders
  # a diverged solver's NaN into a confident 0 - it errors LOUD with
  # the scoped message (qchisq(0.5,8e6)+, qf(0.025,1e9,1e9) were 0
  # where R finite); (b) gamma-p's series truncation (R routes
  # x~shape to ppois_asymp/bd0 - pgamma.c decoded; our series needs
  # > 1000 iterations from shape ~2e4 at x=shape) is the DOCUMENTED
  # boundary with measured cells; (c) the finite imprecision band is
  # pinned per-cell. Values verified live vs R 4.5.0 (rprobe.R c2_
  # cells, 2026-10-13); the jagged 2e5-vs-3e5 onset is the
  # seed-rescue's shape dependence (recorded in api-contract (b)).
  (defn rel [x r] (math/abs (/ (- x r) r)))
  (each o [{:impl :pure} {}]
    # the finite imprecision band (interior p, huge df)
    (is (< (rel (dist/qchisq 0.5 1e5 o) 99999.333334123468) 1e-7)
        "qchisq(.5,1e5) imprecision boundary (rel 4.56e-8, documented)")
    (is (< (rel (dist/qchisq 0.5 2e5 o) 199999.3333337284) 1e-5)
        "qchisq(.5,2e5) imprecision boundary (rel 6.29e-6, documented)")
    (is (< (rel (dist/qchisq 0.5 3e5 o) 299999.3333335967) 1e-12)
        "qchisq(.5,3e5) inside the band (the jagged onset)")
    (is (< (rel (dist/qchisq 0.5 2e6 o) 1999999.3333333728) 1e-7)
        "qchisq(.5,2e6) imprecision boundary (rel 7.83e-8)")
    (is (< (rel (dist/qchisq 0.5 4e6 o) 3999999.333333353) 1e-5)
        "qchisq(.5,4e6) imprecision boundary (rel 7.92e-6, finite)")
    (is (< (rel (dist/qchisq 0.025 1e7 o) 9991236.6690538954) 1e-4)
        "qchisq(.025,1e7) near-0 tail boundary (rel 8.43e-5)")
    # gamma-p's own truncation onset (the p-lane boundary; the q-lane
    # rides the seed-rescue so it does NOT throw here)
    (is (< (rel (dist/pgamma 1e4 1e4 1 o) 0.50132980833995511) 1e-11)
        "pgamma(1e4,1e4) exact-class (series converges, rel 5.1e-12)")
    (is (< (rel (dist/pgamma 1e5 1e5 1 o) 0.50042052211036514) 2e-3)
        "pgamma(1e5,1e5) TRUNCATION boundary (rel 1.58e-3, documented)")
    (is (< (rel (dist/pgamma 1e6 1e6 1 o) 0.50013298076087254) 0.35)
        "pgamma(1e6,1e6) TRUNCATION boundary (rel 0.32, documented)")
    # the qf fudge lanes that stay finite
    (is (< (rel (dist/qf 0.95 1e9 1e9 o) 1.0000735612274694) 1e-10)
        "qf(.95,1e9,1e9) finite (rel 3.5e-11)")
    (is (< (rel (dist/qf 0.5 1e6 1e6 o) 0.99999933333341229) 5e-10)
        "qf(.5,1e6,1e6) finite (rel 1.32e-10)"))
  # the LOUD class: the diverged solver errors, never a silent 0
  (each o [{:impl :pure} {}]
    (each [th lbl] [[(fn [] (dist/qchisq 0.5 8e6 o)) "qchisq(.5,8e6)"]
                      [(fn [] (dist/qchisq 0.5 1e7 o)) "qchisq(.5,1e7)"]
                      [(fn [] (dist/qchisq 0.5 1e10 o)) "qchisq(.5,1e10)"]
                      [(fn [] (dist/qgamma 0.5 5e6 0.5 o)) "qgamma(.5,5e6,.5)"]
                      [(fn [] (dist/qf 0.025 1e9 1e9 o)) "qf(.025,1e9,1e9)"]]
      (var caught nil)
      (try (th) ([e] (set caught e)))
      (is (and caught (string/find "huge-shape" caught))
          (string lbl " errors loud (was a confident 0)")))))

(deftest qgamma-huge-shape-loud-rows
  # the qgamma
  # huge-shape silent imprecision band - the contract pinned qchisq
  # cells only; these are the q-lane cells. Measured live vs R 4.5.0
  # 2026-10-13: gamma-p's series truncation begins ~shape 2e4 (>1000
  # iterations at x~shape), the error climbs through the band to
  # rel ~4e-4 at 2e6-3.5e6, and the loud guard fires from shape 4e6
  # (qchisq df 8e6 = shape 4e6, the documented onset). Both sides of
  # the onset pinned; the guard itself untouched.
  (defn rel [x r] (math/abs (/ (- x r) r)))
  (each o [{:impl :pure} {}]
    # the truncation onset (the series' first visible cells)
    (is (< (rel (dist/qgamma 0.5 2e4 1 o) 19999.666667654339) 1e-12)
        "qgamma(.5,2e4) onset (rel 1.2e-13, series just converging)")
    (is (< (rel (dist/qgamma 0.5 3e4 1 o) 29999.666667325109) 1e-9)
        "qgamma(.5,3e4) onset (rel 6.6e-11, first >1e-12 cell)")
    (is (< (rel (dist/qgamma 0.5 5e4 1 o) 49999.666667061734) 1e-7)
        "qgamma(.5,5e4) in-band (rel 4.6e-8)")
    # the band proper: silent finite, error growing with shape
    (is (< (rel (dist/qgamma 0.5 1e5 1 o) 99999.666666864199) 1e-5)
        "qgamma(.5,1e5) in-band (rel 6.3e-6)")
    (is (< (rel (dist/qgamma 0.5 2e5 1 o) 199999.66666676544) 1e-4)
        "qgamma(.5,2e5) in-band (rel 7.1e-5)")
    (is (< (rel (dist/qgamma 0.5 5e5 1 o) 499999.66666670615) 4e-4)
        "qgamma(.5,5e5) in-band (rel 2.8e-4)")
    (is (< (rel (dist/qgamma 0.5 1e6 1 o) 999999.66666668642) 4.5e-4)
        "qgamma(.5,1e6) in-band (rel 4.0e-4, the jitter class)")
    (is (< (rel (dist/qgamma 0.5 2e6 1 o) 1999999.6666666765) 4.5e-4)
        "qgamma(.5,2e6) in-band (rel 4.2e-4)")
    (is (< (rel (dist/qgamma 0.5 3e6 1 o) 2999999.666666673) 4.5e-4)
        "qgamma(.5,3e6) in-band (rel 4.0e-4)")
    (is (< (rel (dist/qgamma 0.5 3.5e6 1 o) 3499999.6666666721) 4.5e-4)
        "qgamma(.5,3.5e6) the last silent cell (rel 4.0e-4)")
    # the loud side of the onset (the guard fires, never a laundered 0)
    (each [th lbl] [[(fn [] (dist/qgamma 0.5 4e6 1 o)) "qgamma(.5,4e6)"]
                    [(fn [] (dist/qgamma 0.5 1e7 1 o)) "qgamma(.5,1e7)"]]
      (var caught nil)
      (try (th) ([e] (set caught e)))
      (is (and caught (string/find "huge-shape" caught))
          (string lbl " errors loud at the onset (was a confident 0)")))))


(deftest huge-shape-ladder-rows
  # the (b)
  # wording "interior p diverges loud" overstated the loud region.
  # Measured p-ladders at huge shape (live R 4.5.0, 2026-10-13): the
  # loud set is a NARROW p-band centred near p=0.5 that WIDENS with
  # shape (shape 1e7: p in [0.46,0.51] throws; 2e7: [0.44,0.50] throws)
  # - cells BELOW the band stay silent with in-band imprecision, and
  # interior p at/above the band's upper edge stays FINITE (rel
  # ~1e-9 at the band edge decaying to 1e-11..1e-13 by p=0.55).
  # Rows pin both sides at both shapes.
  (defn rel [x r] (math/abs (/ (- x r) r)))
  (each o [{:impl :pure} {}]
    (is (< (rel (dist/qgamma 0.44 1e7 1 o) 9999522.2676956169) 1e-3)
        "qgamma(.44,1e7) below-band silent (rel 2.5e-4, in-band class)")
    (each [th lbl] [[(fn [] (dist/qgamma 0.46 1e7 1 o)) "qgamma(.46,1e7)"]
                    [(fn [] (dist/qgamma 0.48 1e7 1 o)) "qgamma(.48,1e7)"]
                    [(fn [] (dist/qgamma 0.5 1e7 1 o)) "qgamma(.5,1e7)"]
                    [(fn [] (dist/qgamma 0.44 2e7 1 o)) "qgamma(.44,2e7)"]
                    [(fn [] (dist/qgamma 0.5 2e7 1 o)) "qgamma(.5,2e7)"]]
      (var caught nil)
      (try (th) ([e] (set caught e)))
      (is (and caught (string/find "huge-shape" caught))
          (string lbl " in the loud band, errors")))
    # upper interior p: finite past the band edge, error decaying away
    (is (< (rel (dist/qgamma 0.51 1e7 1 o) 10000078.941723162) 1e-9)
        "qgamma(.51,1e7) upper interior finite (rel 2.6e-10, band-edge)")
    (is (< (rel (dist/qgamma 0.52 1e7 1 o) 10000158.267058615) 1e-10)
        "qgamma(.52,1e7) upper interior R-exact (rel 6.0e-12)")
    (is (< (rel (dist/qgamma 0.55 1e7 1 o) 10000397.047992457) 1e-10)
        "qgamma(.55,1e7) upper interior R-exact (rel 8.0e-12)")
    (is (< (rel (dist/qgamma 0.51 2e7 1 o) 20000111.778441038) 1e-8)
        "qgamma(.51,2e7) upper interior finite (rel 9.2e-9, band-edge)")
    (is (< (rel (dist/qgamma 0.52 2e7 1 o) 20000223.961146832) 1e-8)
        "qgamma(.52,2e7) upper interior finite (rel 1.9e-9, band-edge)")
    (is (< (rel (dist/qgamma 0.55 2e7 1 o) 20000561.646552231) 1e-10)
        "qgamma(.55,2e7) upper interior R-exact (rel 1.8e-13)")))


(deftest qbeta-floor-gap-rows
  # the S4
  # floor's documented a-domain ("the a < ~1e-4 medians") misses the
  # boundary band around a=b=0.05, measured live vs R 4.5.0
  # (2026-10-13, p=1e-16 lane): above (a=b=0.1) the log-x lane is
  # R-exact (rel 2.9e-14); AT a=b=0.05 ours underflows the subnormal
  # root (9.71e-315) where R keeps 5.5626846462680035e-309 (rel 1.0 -
  # the seed underflows between grid and root); below (a=b<=0.04)
  # ours returns the bracket floor 8.88e-16 where R underflows to 0
  # (both unreferenced no-answers). Accept-knowingly, pinned per
  # cell; no code change (the honest measured boundary).
  (each o [{:impl :pure} {}]
    (is (< (math/abs (/ (- (dist/qbeta 1e-16 0.1 0.1 o) 8.8692806555504202e-158)
                       8.8692806555504202e-158)) 1e-12)
        "qbeta(1e-16,.1,.1) R-exact above the gap (rel 2.9e-14)")
    (is (close 9.7118236031466185e-315 (dist/qbeta 1e-16 0.05 0.05 o))
        "qbeta(1e-16,.05,.05) the gap cell accept-knowingly (ours subnormal 9.71e-315, R 5.56e-309)")
    (is (close 8.8817841970012523e-16 (dist/qbeta 1e-16 0.04 0.04 o))
        "qbeta(1e-16,.04,.04) below the gap: the bracket floor where R underflows to 0")
    (is (close 8.8817841970012523e-16 (dist/qbeta 1e-16 0.02 0.02 o))
        "qbeta(1e-16,.02,.02) the floor class stable")))


(deftest logp-rejection
  # :
  # :log-p was silently ignored by every q family - the p stayed raw,
  # so a legitimate R log.p call (qgamma(log(0.5),2,1,log.p=TRUE))
  # corrupted to an out-of-domain input and NaN-poisoned with no
  # error. Loud rejection now (the  sibling-option class):
  # the VALIDATION message names the option and the fix.
  (each o [{:impl :pure} {}]
    (each [th lbl] [[(fn [] (dist/qnorm (math/log 0.5) 0 1 (merge o {:log-p true}))) "qnorm"]
                    [(fn [] (dist/qlnorm (math/log 0.5) 0 1 (merge o {:log-p true}))) "qlnorm"]
                    [(fn [] (dist/qexp (math/log 0.5) 2 (merge o {:log-p true}))) "qexp"]
                    [(fn [] (dist/qweibull (math/log 0.5) 1.5 1 (merge o {:log-p true}))) "qweibull"]
                    [(fn [] (dist/qgamma (math/log 0.5) 2 1 (merge o {:log-p true}))) "qgamma"]
                    [(fn [] (dist/qbeta (math/log 0.5) 2 3 (merge o {:log-p true}))) "qbeta"]
                    [(fn [] (dist/qchisq (math/log 0.5) 4 (merge o {:log-p true}))) "qchisq"]
                    [(fn [] (dist/qt (math/log 0.5) 7 (merge o {:log-p true}))) "qt"]
                    [(fn [] (dist/qf (math/log 0.5) 3 5 (merge o {:log-p true}))) "qf"]
                    [(fn [] (dist/qbinom (math/log 0.5) 48 0.5 (merge o {:log-p true}))) "qbinom"]
                    [(fn [] (dist/qpois (math/log 0.5) 10 (merge o {:log-p true}))) "qpois"]
                    [(fn [] (dist/qnbinom (math/log 0.5) 8 0.3 (merge o {:log-p true}))) "qnbinom"]
                    [(fn [] (dist/qtukey (math/log 0.5) 4 12 (merge o {:log-p true}))) "qtukey"]]
      (var caught nil)
      (try (th) ([e] (set caught e)))
      (is (and caught (string/find ":log-p is not supported" caught))
          (string lbl " rejects :log-p loud (was a silent NaN)")))
    # the deep-cell NaN class (raw log(1e-10) = -23, deep out of
    # domain): both lanes, the upper lane too
    (each [th lbl] [[(fn [] (dist/qgamma (math/log 1e-10) 2 1 (merge o {:log-p true})))
                     "qgamma(log 1e-10) deep cell"]
                    [(fn [] (dist/qgamma (math/log 0.5) 2 1 (merge o {:log-p true :lower-tail false})))
                     "qgamma(log .5) upper lane"]
                    [(fn [] (dist/qbinom (math/log 0.5) 48 0.5 (merge o {:log-p true :lower-tail false})))
                     "qbinom(log .5) upper lane"]]
      (var caught nil)
      (try (th) ([e] (set caught e)))
      (is (and caught (string/find ":log-p is not supported" caught))
          (string lbl " rejects loud, never a silent NaN")))
    # the unaffected neighbour:lower-tail alone still decodes
    (is (< (math/abs (/ (- (dist/qgamma (- 1 1e-15) 2 1 (merge o {:lower-tail false}))
                           4.4703483581543002e-08) 4.4703483581543002e-08)) 1e-12)
        "qgamma(1-1e-15,up) neighbour unaffected")))

(deftest qbeta-far-tail-rows
  # :
  # solve-monotone's 1e-15 ABSOLUTE bracket width strands every
  # direct-lane root below ~1e-15 (qbeta(1e-12,0.5,0.5) was 8.88e-16
  # vs R 2.47e-24, rel 3.6e8); the fix ports R qbeta.c's use_log_x
  # Newton (lines 444-509): u = log(x), y = log I_x(a,b) computed
  # DIRECTLY from beta-inc's log-space combination (- never a
  # linear-scale I that underflows), R's w step formula, engaged
  # when the leading-order seed x0 = (tgt*a*B(a,b))^(1/a) lands
  # below 1e-13 (R's own cut is seed < e^-5; ours is conservative -
  # every verified mid cell keeps the bracket path). Values
  # verified live vs R 4.5.0 (rprobe.R c4_ cells, 2026-10-13).
  (defn rel [x r] (math/abs (/ (- x r) r)))
  (each o [{:impl :pure} {}]
    # the bar cells: the arcsine far tail, R-exact
    (is (< (rel (dist/qbeta 1e-12 0.5 0.5 o) 2.4674011002723402e-24) 1e-12)
        "qbeta(1e-12,.5,.5) = R (was the 8.88e-16 floor, rel 3.6e8)")
    (is (< (rel (dist/qbeta 1e-16 0.5 0.5 o) 2.4674011002723404e-32) 1e-12)
        "qbeta(1e-16,.5,.5) = R (was rel 3.6e16)")
    (is (< (rel (dist/qbeta 1e-8 0.5 0.5 o) 2.4674011002723406e-16) 1e-12)
        "qbeta(1e-8,.5,.5) = R (was rel 1.34e-3)")
    (is (< (rel (dist/qbeta 1e-14 0.5 0.5 o) 2.4674011002723408e-28) 1e-12)
        "qbeta(1e-14,.5,.5) = R")
    # the a<1 class (x ~ p^{1/a}, 20 orders of floor)
    (is (< (rel (dist/qbeta 1e-6 0.05 0.05 o) 9.7118236027201415e-115) 1e-12)
        "qbeta(1e-6,.05,.05) = R (was 8.88e-16, rel 9e98)")
    # qf's beta composition rides the fix finite garbage
    (is (and (< (math/abs (dist/qf 5e-324 4 7 {:lower-tail false})) math/inf)
            (< (rel (dist/qf 5e-324 4 7 (merge {:lower-tail false} o)) 6.3514747985027126e+92) 1e-12))
        "qf(5e-324,4,7,up) = R 6.35e92 (was the finite garbage 7.57e15)")
    # the pinned  cells must hold (bracket path byte-stable)
    (is (close 2.5819889063604971e-09 (dist/qbeta 1e-16 2 5 o))
        "qbeta(1e-16,2,5) near-0 side intact (the pinned cell)")
    (is (close 0.99999999643651671 (dist/qbeta 1e-16 3.5 2 {:lower-tail false}))
        "qbeta(1e-16,3.5,2,up) mirror intact (the  bar cell)")
    (is (close 0.9995497181989903 (dist/qbeta (- 1 1e-16) 2 5 o))
        "qbeta(1-1e-16,2,5) mirror intact")
    # interiors on the bracket path
    (is (close 0.64123457899767478 (dist/qbeta 0.975 2 5 o))
        "qbeta(.975,2,5) interior intact")
    # the DOCUMENTED FLOOR (R-warns extreme class): a=1e-8, b=1 puts
    # the median at x = 0.5^(1e8) ~ e^-6.9e7, below the double grid;
    # R floors at 2^-1024 and warns "full precision may not have been
    # achieved" - OUR floor is the bracket's 8.88e-16 (the seed
    # underflows to 0, the log-x lane disengages); both unreferenced
    (is (close 8.8818596130072009e-16 (dist/qbeta 0.5 1e-8 1 o))
        "qbeta(.5,1e-8,1) the documented floor, ours 8.88e-16 (R floors at 2^-1024+warns)")
    # ACCEPT-KNOWINGLY : the arcsine median at
    # extreme (a,b) - OURS is the symmetric-correct 0.5; R's 1.87e-41
    # comes with R's own full-precision warning (reference divergence)
    (is (< (rel (dist/qbeta 0.5 1e-8 1e-8 o) 0.5) 1e-7)
        "qbeta(.5,1e-8,1e-8) = the arcsine-correct median 0.5 (accept-knowingly; R 1.87e-41 WARNS)")))

(deftest discrete-search-rows
  # the discrete lane decoded from qDiscrete_search.h -
  # the 8eps left-continuity fuzz on the RAW p (lower shrinks p,
  # upper grows p only when 1-p > 32eps), the Cornish-Fisher seed
  # on the flag-aware qnorm, and do_search on the FLAG-THREADED CDF
  # (lower walks P(K<=y), upper walks P(K>y) DIRECTLY - the old
  # complement walk's threshold 1-v+8eps saturated above 1 at
  # v < ~2e-15 so the walk never advanced, and the upper lane's
  # materialized target rounded to 1.0 at subnormal p and walked
  # to inf). pbinom's upward recurrence rounds ABOVE 1 in the
  # saturated upper region (an impossible CDF value) - recomputed
  # from the backward tail there. Every value verified live vs
  # R 4.5.0 (rprobe.R S5 cells, 2026-10-13).
  (each o [{:impl :pure} {}]
    # (a) the qpois lower walk (was the normal-approx seed)
    (is (= 30 (dist/qpois 1e-16 100 o))
        "qpois(1e-16,100) = 30 (was the seed 18)")
    (is (= 751 (dist/qpois 1e-16 1000 o))
        "qpois(1e-16,1000) = 751 (was 740)")
    (is (= 9189 (dist/qpois 1e-16 10000 o))
        "qpois(1e-16,10000) = 9189 (was 9178)")
    (is (= 32 (dist/qpois 2e-15 100 o))
        "qpois(2e-15,100) = 32 (was 30)")
    # (b) the upper-lane walks take the target at the right index
    #    (pbinom's impossible >1 values were the root)
    (is (= 968 (dist/qbinom 1e-16 1000 0.9 (merge {:lower-tail false} o)))
        "qbinom(1e-16,1000,0.9,up) = 968 (was 960)")
    (is (= 89 (dist/qbinom 1e-16 100 0.5 (merge {:lower-tail false} o)))
        "qbinom(1e-16,100,0.5,up) = 89 (was 88)")
    (is (= 47 (dist/qbinom 1e-16 50 0.4 (merge {:lower-tail false} o)))
        "qbinom(1e-16,50,0.4,up) = 47 (was 46)")
    (is (= 529 (dist/qnbinom 1e-16 100 0.3 (merge {:lower-tail false} o)))
        "qnbinom(1e-16,100,0.3,up) = 529 (was 528)")
    # the pbinom CDF row (the impossible 1.0000000000001048)
    (is (close 0.99999999999961542 (dist/pbinom 960 1000 0.9 o))
        "pbinom(960,1000,0.9) = R (was 1.0000000000001048, impossible)")
    (is (= 1 (dist/pbinom 999 1000 0.9 o))
        "pbinom(999,1000,0.9) = 1 (the saturated pole intact)")
    # (c) the upper-lane termination at tiny p (was inf)
    (is (= 238 (dist/qpois 1e-300 5 (merge {:lower-tail false} o)))
        "qpois(1e-300,5,up) = 238 (was inf)")
    (is (= 2063 (dist/qnbinom 1e-300 10 0.3 (merge {:lower-tail false} o)))
        "qnbinom(1e-300,10,0.3,up) = 2063 (was inf)")
    # (d) the fuzz constant 16eps -> 32eps (_xf_=4, source-verified
    #     qDiscrete_search.h:184,202) - THE FUZZ COVERAGE ROW: its
    #     removal mutation must fail here
    (is (= 1 (dist/qbinom (- 1 5e-15) 48 0.5 (merge {:lower-tail false} o)))
        "qbinom(1-5e-15,48,0.5,up) = 1 (the 32eps window; was 0 at 16eps)")
    # (e) the flipped-lane cells (these
    #     now match R via the flag-aware walk)
    (is (= 33 (dist/qpois 1e-16 5 (merge {:lower-tail false} o)))
        "qpois(1e-16,5,up) = 33 (the far-tail cell, was 32)")
    (is (= 167 (dist/qnbinom 1e-16 10 0.3 (merge {:lower-tail false} o)))
        "qnbinom(1e-16,10,0.3,up) = 167 (the far-tail cell, was 166)")
    (is (= 30 (dist/qpois (- 1 1e-16) 100 (merge {:lower-tail false} o)))
        "qpois(1-1e-16,100,up) = 30 (the near-1 upper cell, was 29)")
    # interiors - both lanes, all three families (byte-stable class)
    (is (= 10 (dist/qpois 0.5 10 o)) "qpois(.5,10) interior intact")
    (is (= 120 (dist/qpois 0.975 100 o)) "qpois(.975,100) interior intact")
    (is (= 120 (dist/qpois 0.025 100 (merge {:lower-tail false} o)))
        "qpois(.025,100,up) interior intact (the flip mirror)")
    (is (= 5 (dist/qbinom 0.5 10 0.5 o)) "qbinom(.5,10,.5) interior intact")
    (is (= 95 (dist/qbinom 0.975 100 0.9 o)) "qbinom(.975,100,.9) interior intact")
    (is (= 22 (dist/qnbinom 0.5 10 0.3 o)) "qnbinom(.5,10,.3) interior intact")))



(deftest qlnorm-endpoint-rows
  # the (= v 0)/(= v 1)
  # endpoint branches were NOT flag-aware - the upper lane returned
  # the LOWER lane's endpoints (0/inf swapped). R: qlnorm = exp(qnorm)
  # maps the flag-aware -inf/inf to 0/Inf on both lanes. Endpoint
  # class verified live vs R 4.5.0: qlnorm(0,0,1)=0, qlnorm(1,0,1)=Inf,
  # qlnorm(0,0,1,lt=F)=Inf, qlnorm(1,0,1,lt=F)=0 - same at any
  # meanlog/sdlog (qnorm endpoints are parameter-free). Exact = rows,
  # not close - inf/0 pinning cannot pass vacuously (S8's close()
  # degeneracy class).
  (each o [{:impl :pure} {}]
    (is (= 0 (dist/qlnorm 0 0 1 o)) "qlnorm(0,0,1) = 0 (lower, intact)")
    (is (= math/inf (dist/qlnorm 1 0 1 o)) "qlnorm(1,0,1) = inf (lower, intact)")
    (is (= math/inf (dist/qlnorm 0 0 1 (merge {:lower-tail false} o)))
       "qlnorm(0,0,1,up) = inf (was 0 - inverted)")
    (is (= 0 (dist/qlnorm 1 0 1 (merge {:lower-tail false} o)))
       "qlnorm(1,0,1,up) = 0 (was inf - inverted)")
    # the endpoint class is parameter-free - non-default meanlog/sdlog
    (is (= math/inf (dist/qlnorm 0 0.1 0.5 (merge {:lower-tail false} o)))
       "qlnorm(0,.1,.5,up) = inf (endpoint independent of params)")
    (is (= 0 (dist/qlnorm 1 0.1 0.5 (merge {:lower-tail false} o)))
       "qlnorm(1,.1,.5,up) = 0 (endpoint independent of params)")
    # interiors byte-stable on both lanes (the flip delegates to
    # qnorm - only the endpoints bypassed it)
    (is (close 7.0990713842313333 (dist/qlnorm 0.025 0 1 (merge {:lower-tail false} o)))
       "qlnorm(.025,0,1,up) interior intact")
    (is (close 0.00026865508415775056 (dist/qlnorm 1e-16 0 1 o))
       "qlnorm(1e-16,0,1) interior intact")))



(deftest handle-gate-rows
  # the FULL-marker gate (dlnorm's
  # shape, ) wired into the 13 kernel-backed d/p families -
  # under the default (num) impl a partial-marker struct fell into
  # num/stage and died with janet-num internals noise ("can not
  # convert nil nil to 64 bit signed integer", "could not find
  # method :* for 8 or :r* for nil") - a staging crash, not an
  # actionable error. Bar rows: the two reproduced cells
  # (dpois {:n 2}, dexp {:a 1}); the class row sweeps the full
  # partial-marker set over every gated family, both impls.
  # Audit found 13: dnorm pnorm dt dchisq pchisq df dbeta pbeta dgamma
  # pgamma dpois ppois dexp (composition families + dlnorm/plnorm
  # were already actionable).
  (each o [{:impl :pure} {}]
    # the two reproduced error cells (verbatim from the recorded report)
    (var caught nil)
    (try (dist/dpois {:n 2} 5 o) ([e] (set caught e)))
    (is (and caught (string/find "caller misuse" (string caught)))
       "dpois({:n 2}) -> actionable misuse error (was 'can not convert nil nil to 64 bit signed integer')")
    (set caught nil)
    (try (dist/dexp {:a 1} 1 o) ([e] (set caught e)))
    (is (and caught (string/find "caller misuse" (string caught)))
       "dexp({:a 1}) -> actionable misuse error (was 'could not find method :* for 8 or :r* for nil')")
    # the class row: every gated family, every partial marker
    (let [fams [[:dnorm (fn [x] (dist/dnorm x 0 1 o))]
               [:pnorm (fn [x] (dist/pnorm x 0 1 o))]
               [:dt (fn [x] (dist/dt x 5 o))]
               [:dchisq (fn [x] (dist/dchisq x 5 o))]
               [:pchisq (fn [x] (dist/pchisq x 5 o))]
               [:df (fn [x] (dist/df-f x 4 7 o))]
               [:dbeta (fn [x] (dist/dbeta x 2 3 o))]
               [:pbeta (fn [x] (dist/pbeta x 2 3 o))]
               [:dgamma (fn [x] (dist/dgamma x 2 1 o))]
               [:pgamma (fn [x] (dist/pgamma x 2 1 o))]
               [:dpois (fn [x] (dist/dpois x 5 o))]
               [:ppois (fn [x] (dist/ppois x 5 o))]
                [:dexp (fn [x] (dist/dexp x 1 o))]]]
      (each [nm f] fams
       (each bad [{:n 2} {:bytes 8} {:a 1} {:n 2 :bytes 8 :data "junk"}]
         (set caught nil)
         (try (f bad) ([e] (set caught e)))
         (is (and caught (string/find "caller misuse" (string caught)))
             (string/format "%s %q -> actionable misuse error" nm bad))))))
  # full packed handles stay INTENDED kernel inputs - the gate must
  # not touch them (the  ppois/dpois gate shape)
  (let [h (((num/kernels) 'pack-f64) @[1 2 3])]
    (is (= 3 (length (dist/ppois h 5 {})))
       "ppois on a real packed handle still stages (gate does not reject handles)")))



(deftest qt-as396-rows
  # the AS 396 port replaces the t^2=F(1,df)
  # composition; the large-df ladder (rel < 1e-10 bar) + interiors +
  # far tails, every value verified live vs R 4.5.0 (2026-10-10).
  # 1e8/1e9 sit below the polish guard (our dt/beta-inc precision
  # boundary, recorded in the CR) and ride the seed alone.
  (each o [{:impl :pure} {}]
    # the ladder: was rel 3e-6 at 1e10, 2.1x garbage at 1e16, 215x at 1e20
    (is (close 1.9599640082627663 (dist/qt 0.975 1e8 o))
        "qt(.975,1e8) = R (seed lane)")
    (is (close 1.9599639869123249 (dist/qt 0.975 1e9 o))
        "qt(.975,1e9) = R (seed lane)")
    (is (close 1.9599639847772807 (dist/qt 0.975 1e10 o))
        "qt(.975,1e10) = R EXACT (bar cell, was rel 3e-6)")
    (is (close 1.9599639845424257 (dist/qt 0.975 1e12 o))
        "qt(.975,1e12) = R EXACT")
    (is (close 1.9599639845400538 (dist/qt 0.975 1e16 o))
        "qt(.975,1e16) = R (was 4.2147, 2.1x garbage)")
    (is (close 1.9599639845400536 (dist/qt 0.975 1e18 o))
        "qt(.975,1e18) = R (was 42.15, 21x)")
    (is (close 1.9599639845400536 (dist/qt 0.975 1e20 o))
        "qt(.975,1e20) = R (was 421.47, 215x)")
    # interiors: the exact branches + the AS 396 middle/P>P0 lanes
    (is (close 2.5705818356363137 (dist/qt 0.975 5 o))
        "qt(.975,5) = R")
    (is (close 6.0166631044279288 (dist/qt 0.975 1.5 o))
        "qt(.975,1.5) = R")
    (is (close 2.1788128296672284 (dist/qt 0.975 12 o))
        "qt(.975,12) = R EXACT")
    (is (close 4.7852896286383331 (dist/qt 0.999 7 o))
        "qt(.999,7) = R")
    # far tails (the middle branch, y0 <= 0.05+a): both the pinned
    #  cell and the band below it
    (is (close 2434.2914303699285 (dist/qt (- 1 1e-16) 5 o))
        "qt(1-1e-16,5) = R (pinned  cell)")
    (is (close 393.95870264999013 (dist/qt (- 1 1e-12) 5 o))
        "qt(1-1e-12,5) = R")
    (is (close 62.404506048184757 (dist/qt (- 1 1e-8) 5 o))
        "qt(1-1e-8,5) = R")
    (is (close -62.404506110967276 (dist/qt 1e-8 5 o))
        "qt(1e-8,5) = R (near-0 lane; R's near-1/near-0 doubles differ by rel 1e-9 - path noise)")))

(deftest qf-fudge-rows
  # R qf.c's extreme-df fudge ported
  # (df > 4e5); every value verified live vs R 4.5.0 (2026-10-10).
  # The df1==Inf==df2 cell rides the fudge's 1.0. The 1e-12/1e-16
  # LOWER-tail qf cells (beta-inc complement) are S4's slice.
  (each o [{:impl :pure} {}]
    # the fudge cells (were rel 4.9e-6 / 1.1e-6 off)
    (is (close 2.5665003988060051 (dist/qf 0.975 5 1e6 o))
        "qf(.975,5,1e6) = R (fudge lane, was 4.9e-6)")
    (is (close 0.1662423226973325 (dist/qf 0.025 5 1e6 o))
        "qf(.025,5,1e6) = R (fudge lane)")
    (is (close 6.0153153768227856 (dist/qf 0.975 1e6 5 o))
        "qf(.975,1e6,5) = R (fudge lane, was 1.1e-6)")
    (is (close 0.38963563008414986 (dist/qf 0.025 1e6 5 o))
        "qf(.025,1e6,5) = R (fudge lane, was 4.9e-6)")
    # 4e5-boundary: the qbeta lane holds ~1e-11 rel here (S4's polish
    # tightened it from 6.9e-11; the residual is the (1/x - 1)
    # composition's absolute floor at x ~= 1-3e-5) - rel 1e-9
    (is (<= (math/abs (- (dist/qf 0.975 5 4e5 o) 2.5665319429180045))
            (* 1e-9 2.5665319429180045))
        "qf(.975,5,4e5) = R to rel 1e-9 (fudge boundary, qbeta lane)")
    (is (close 1.0027737014679261 (dist/qf 0.975 1e6 1e6 o))
        "qf(.975,1e6,1e6) = R (both sides huge, df1<=df2 lane)")
    (is (= 1.0 (dist/qf 0.975 math/inf math/inf o))
        "qf(.975,Inf,Inf) = 1 (R's both-Inf cell)")
    # far tails both sides + interiors (the qbeta lane unchanged)
    (is (<= (math/abs (- (dist/qf (- 1 1e-12) 4 7 o) 7212.9945286892053))
            (* 1e-9 7212.9945286892053))
        "qf(1-1e-12,4,7) = R to rel 1e-9 (upper far tail, qbeta lane)")
    (is (<= (math/abs (- (dist/qf (- 1 1e-16) 4 7 o) 97300.454482355402))
            (* 1e-9 97300.454482355402))
        "qf(1-1e-16,4,7) = R to rel 1e-9 (pinned  cell, qbeta lane)")
    (is (close 6.2365030831734725e-05 (dist/qf 1e-8 4 7 o))
        "qf(1e-8,4,7) = R (lower far tail)")
    (is (close 517.11998956158232 (dist/qf (- 1 1e-8) 4 7 o))
        "qf(1-1e-8,4,7) = R")
    (is (close 5.5225943453085495 (dist/qf 0.975 4 7 o))
        "qf(.975,4,7) = R (interior intact)")))

(deftest lower-tail-rows
  # :lower-tail family-wide - all 13
  # families, both tails, every value verified live vs R 4.5.0
  # (2026-10-11). The old code silently computed the lower tail on
  # ten of thirteen (qt(0.025,5,lt-false) was -2.5705818, R +2.5705818)
  (each o [{:impl :pure} {}]
    (def u (merge o {:lower-tail false}))
    (is (close 1.9599639845400538 (dist/qnorm 0.025 0 1 u))
        "qnorm(.025,up) = +1.9599639845")
    (is (close 7.0990713842313333 (dist/qlnorm 0.025 0 1 u))
        "qlnorm(.025,0,1,up) = 7.0990713842 = R")
    (is (close 3.6888794541139363 (dist/qexp 0.025 1 u))
        "qexp(.025,1,up) = 3.6888794541 = R (the bar cell)")
    (is (close 1.9206455826398414 (dist/qweibull 0.025 2 1 u))
        "qweibull(.025,2,1,up) = 1.9206455826 = R (the bar cell)")
    (is (close 3.6123438338619791 (dist/qgamma 0.025 3 2 u))
        "qgamma(.025,3,2,up) = 3.6123438339 = R")
    (is (close 0.64123457899767589 (dist/qbeta 0.025 2 5 u))
        "qbeta(.025,2,5,up) = R")
    (is (close 12.832501994030029 (dist/qchisq 0.025 5 u))
        "qchisq(.025,5,up) = 12.8325019940 = R")
    (is (close 2.5705818356363150 (dist/qt 0.025 5 u))
        "qt(.025,5,up) = +2.5705818356 = R (the bar cell, was -2.5706)")
    (is (close 5.5225943453085291 (dist/qf 0.025 4 7 u))
        "qf(.025,4,7,up) = 5.5225943453 = R (the bar cell, was 0.1102)")
    (is (= 6 (dist/qbinom 0.025 10 0.3 u))
        "qbinom(.025,10,.3,up) = 6 = R")
    (is (= 10 (dist/qpois 0.025 5 u))
        "qpois(.025,5,up) = 10 = R")
    (is (= 43 (dist/qnbinom 0.025 10 0.3 u))
        "qnbinom(.025,10,.3,up) = 43 = R")
    # the flipped upper tails at p=0.975 (the other side)
    (is (close -1.9599639845400536 (dist/qnorm 0.975 0 1 u))
        "qnorm(.975,up) = -1.9599639845 = R")
    (is (close 0.025317807984289897 (dist/qexp 0.975 1 u))
        "qexp(.975,1,up) = R")
    (is (close -2.5705818356363199 (dist/qt 0.975 5 u))
        "qt(.975,5,up) = -2.5705818356 = R")
    (is (close 0.11020326820103399 (dist/qf 0.975 4 7 u))
        "qf(.975,4,7,up) = R")
    (is (= 0 (dist/qbinom 0.975 10 0.3 u))
        "qbinom(.975,10,.3,up) = 0 = R")
    (is (= 1 (dist/qpois 0.975 5 u))
        "qpois(.975,5,up) = 1 = R")
    # qtukey (R qtukey.c:140's R_DT_qIv flip; the 13th family): the
    # flipped lane's band is NOT the lower lane's - the
    # near-1 upper zone sits at the secant's eps 1e-4 floor where
    # BOTH impls are noise (ours ~0.44x R across 1-1e-8..1-1e-16,
    # qtukey(1-1e-16,4,12,up) 1.76e-15 vs R 4.04e-15), and at p =
    # 1-1e-16 the lower lane is a CONVENTION SPLIT (ours 0, R NaN
    # "convergence failed" - both no-answers); rows pin OUR values
    # with the class reason, never an R-agreement claim (cells
    # measured 2026-10-13)
    (is (<= (math/abs (- (dist/qtukey 0.025 4 12 u) 4.761395574956687))
            (* 1e-8 4.761395574956687))
        "qtukey(.025,4,12,up) = 4.7613955750 to rel 1e-8 = R (the 13th family)")
    (is (<= (math/abs (- (dist/qtukey 0.975 4 12 u) 0.58561782116619532))
            (* 1e-8 0.58561782116619532))
        "qtukey(.975,4,12,up) = R to rel 1e-8")
    (is (<= (math/abs (- (dist/qtukey 0.025 4 12) 0.5856178211661951))
            (* 1e-8 0.5856178211661951))
        "qtukey(.025,4,12) lower lane unchanged = R to rel 1e-8")
    # the solver-floor class rows - pin OURS, the class
    # documented (never an R-agreement claim: the R side is NaN or
    # floor noise)
    (is (close 1.7613248589497419e-15 (dist/qtukey (- 1 1e-16) 4 12 u))
        "qtukey(1-1e-16,4,12,up) = 1.76e-15 OURS - the secant-floor class (R 4.04e-15, rel 0.56, both noise; api-contract 5p)")
    (is (= 0 (dist/qtukey (- 1 1e-16) 4 12))
        "qtukey(1-1e-16,4,12) = 0 OURS - the convention split (R NaN, convergence failed; both no-answers)")))

(deftest packed-marker-guards
  # the FULL handle marker (n AND bytes
  # AND ptr since  - was data-or-close; partial markers
  # actionable, never a staging crash) + t-test-one's nil mu default
  # (the guard's own remedy)
  (each o [{:impl :pure} {}]
    (var caught nil)
    (each bad [{:n 2} {:bytes 8} {:n 2 :bytes 8 :data "junk"} {:n 2 :bytes 8}]
      (set caught nil)
      (try (dist/dlnorm bad 0 1 o) ([e] (set caught e)))
      (is (and caught (string/find "caller misuse" (string caught)))
          (string/format "dlnorm %q -> the actionable misuse error (was a staging crash)"
                         bad))
      (set caught nil)
      (try (dist/plnorm bad 0 1 o) ([e] (set caught e)))
      (is (and caught (string/find "caller misuse" (string caught)))
          (string/format "plnorm %q -> the actionable misuse error" bad))))
  # t-test-one: nil mu = 0 (the documented remedy, was a nil-
  # arithmetic crash one line past its own error message)
  (let [r (tests/t-test-one @[1 2 3 4 5] nil {})
        r0 (tests/t-test-one @[1 2 3 4 5] 0 {})]
    (is (and (close (r :statistic) (r0 :statistic))
             (close (r :p-value) (r0 :p-value)))
        "t-test-one(xs, nil) = the mu-0 result (the guard's remedy now works)")))

(deftest subnormal-rows
  # subnormal p on the flipped lane -
  # the 1-p complement swallowed 5e-324 to exactly 1.0 (qnorm
  # returned +Inf, qt -Inf). Values live vs R 4.5.0 (2026-10-11);
  # the qnorm cells carry R's r>27 asymptotic-branch precision
  # (rel ~1.3e-9, the documented class - same as the forward side)
  (each o [{:impl :pure} {}]
    (def u (merge o {:lower-tail false}))
    (is (<= (math/abs (- (dist/qnorm 5e-324 0 1 u) 38.467405617144344))
            (* 1e-8 38.467405617144344))
        "qnorm(5e-324,up) = 38.467405617 to rel 1e-8 (was +Inf) - the bar")
    (is (close -38.467405667754072 (dist/qnorm 5e-324 0 1 o))
        "qnorm(5e-324) forward unchanged (rel 1.3e-9, accepted class)")
    (is (close -7.1894860044921933e+64 (dist/qt 5e-324 5 o))
        "qt(5e-324,5) = -7.1894860045e64 = R EXACT (was -Inf) - the bar")
    (is (close 7.1894860044921933e+64 (dist/qt 5e-324 5 u))
        "qt(5e-324,5,up) = +R EXACT")
    (is (close -1.5683925619292198e+60 (dist/qt 1e-300 5 o))
        "qt(1e-300,5) = R EXACT")
    (is (close 744.44007192138122 (dist/qexp 5e-324 1 u))
        "qexp(5e-324,1,up) = 744.440071921 = R EXACT (raw-log lane)")
    (is (close 744.44007192138122 (dist/qweibull 5e-324 1 1 u))
        "qweibull(5e-324,1,1,up) = R EXACT")
    (is (close 3675.8370396426658 (dist/qlnorm (- 1 1e-16) 0 1 o))
        "qlnorm(1-1e-16) = R EXACT")))

(deftest far-tail-remaining-rows
  # the far-tail rows for the families
  # not covered in S1/S2 - near-1 (1-1e-12, 1-1e-16) + the near-0
  # mirrors, all verified live vs R 4.5.0 (2026-10-11), all EXACT.
  # qbeta's far-tail rows land in S4 (beta-inc's near-1 complement
  # is that slice's fix; recorded in the CR resolution)
  (each o [{:impl :pure} {}]
    (is (close 7.0344869100478356 (dist/qnorm (- 1 1e-12) 0 1 o))
        "qnorm(1-1e-12) = R EXACT")
    (is (close 8.2095361516013856 (dist/qnorm (- 1 1e-16) 0 1 o))
        "qnorm(1-1e-16) = R EXACT (the pinned  cell)")
    (is (close -8.2220822161304348 (dist/qnorm 1e-16 0 1 o))
        "qnorm(1e-16) = R EXACT (near-0 side)")
    (is (close 1135.1123480102162 (dist/qlnorm (- 1 1e-12) 0 1 o))
        "qlnorm(1-1e-12) = R EXACT")
    (is (close 3675.8370396426658 (dist/qlnorm (- 1 1e-16) 0 1 o))
        "qlnorm(1-1e-16) = R EXACT")
    (is (close 0.00026865508415775056 (dist/qlnorm 1e-16 0 1 o))
        "qlnorm(1e-16) = R EXACT (near-0 side)")
    (is (close 27.631043237893358 (dist/qexp (- 1 1e-12) 1 o))
        "qexp(1-1e-12,1) = R EXACT")
    (is (close 36.736800569677101 (dist/qexp (- 1 1e-16) 1 o))
        "qexp(1-1e-16,1) = R EXACT")
    (is (close 1.0000000000005e-12 (dist/qexp 1e-12 1 o))
        "qexp(1e-12,1) = R EXACT (near-0 side)")
    (is (close 15.769571621988982 (dist/qweibull (- 1 1e-12) 2 3 o))
        "qweibull(1-1e-12,2,3) = R EXACT")
    (is (close 18.183267174165756 (dist/qweibull (- 1 1e-16) 2 3 o))
        "qweibull(1-1e-16,2,3) = R EXACT")
    (is (close 3.0000000000000004e-08 (dist/qweibull 1e-16 2 3 o))
        "qweibull(1e-16,2,3) = R EXACT (near-0 side)")
    (is (= 10 (dist/qbinom (- 1 1e-12) 10 0.3 o))
        "qbinom(1-1e-12,10,.3) = 10 = R")
    (is (= 10 (dist/qbinom (- 1 1e-16) 10 0.3 o))
        "qbinom(1-1e-16,10,.3) = 10 = R")
    (is (= 0 (dist/qbinom 1e-16 10 0.3 o))
        "qbinom(1e-16,10,.3) = 0 = R (near-0 side)")
    (is (= 27 (dist/qpois (- 1 1e-12) 5 o))
        "qpois(1-1e-12,5) = 27 = R")
    (is (= 31 (dist/qpois (- 1 1e-16) 5 o))
        "qpois(1-1e-16,5) = 31 = R")
    (is (= 0 (dist/qpois 1e-16 5 o))
        "qpois(1e-16,5) = 0 = R (near-0 side)")
    (is (= 136 (dist/qnbinom (- 1 1e-12) 10 0.3 o))
        "qnbinom(1-1e-12,10,.3) = 136 = R")
    (is (= 157 (dist/qnbinom (- 1 1e-16) 10 0.3 o))
        "qnbinom(1-1e-16,10,.3) = 157 = R")
    (is (= 0 (dist/qnbinom 1e-16 10 0.3 o))
        "qnbinom(1e-16,10,.3) = 0 = R (near-0 side)")))

(deftest qbeta-complement-rows
  # qbeta's mirrored solve (R's
  # swap_tail) + swap-side Newton - the near-1 CDF's double spacing
  # equals the complement's scale, so the old direct solve was
  # staircase-blind (qbeta(1e-16,3.5,2,up) was 25% off). All values
  # verified live vs R 4.5.0 (2026-10-11); qf's cells ride the same
  # fix through the composition. These are the S3-deferred qbeta
  # far-tail rows (reorder recorded in the CR resolution)
  (each o [{:impl :pure} {}]
    (def u (merge o {:lower-tail false}))
    # qbeta far tails, both sides + both lanes
    (is (close 0.99721663905675362 (dist/qbeta (- 1 1e-12) 2 5 o))
        "qbeta(1-1e-12,2,5) = R EXACT (was rel 3.1e-8)")
    (is (close 0.9995497181989903 (dist/qbeta (- 1 1e-16) 2 5 o))
        "qbeta(1-1e-16,2,5) = R EXACT (was rel 5.8e-5)")
    (is (close 2.5819897863610085e-07 (dist/qbeta 1e-12 2 5 o))
        "qbeta(1e-12,2,5) = R (near-0 side, was rel 2.8e-9)")
    (is (close 2.5819889063604971e-09 (dist/qbeta 1e-16 2 5 o))
        "qbeta(1e-16,2,5) = R (near-0 side, was rel 3.7e-7)")
    (is (close 0.99999999643651671 (dist/qbeta 1e-16 3.5 2 u))
        "qbeta(1e-16,3.5,2,up) = R EXACT (the bar cell, was 25%)")
    (is (close 0.99999964365157168 (dist/qbeta 1e-12 3.5 2 u))
        "qbeta(1e-12,3.5,2,up) = R EXACT")
    (is (close 0.64123457899767478 (dist/qbeta 0.975 2 5 o))
        "qbeta(.975,2,5) interior intact")
    # qf's far-tail cells through the composition
    (is (close 6.2360957509000059e-09 (dist/qf 1e-16 4 7 o))
        "qf(1e-16,4,7) = R EXACT (the bar cell, was 25%)")
    (is (close 6.2360997182819844e-07 (dist/qf 1e-12 4 7 o))
        "qf(1e-12,4,7) = R EXACT (was rel 3.9e-5)")
    (is (close 7212.9945286892053 (dist/qf (- 1 1e-12) 4 7 o))
        "qf(1-1e-12,4,7) = R (tightened from the S2 rel-1e-9 row)")
    (is (close 97300.454482355402 (dist/qf (- 1 1e-16) 4 7 o))
        "qf(1-1e-16,4,7) = R (tightened from the S2 rel-1e-9 row)")
    # THE DOCUMENTED DEVIATION (user decision on measured evidence):
    # qf(1-1e-12,1,1) solves I_x(.5,.5) = 1e-12 -> x ~= 2.47e-24, eight
    # orders below the linear-scale solver's 1e-15 bracket floor and
    # unreachable by the sqrt-tail Newton (linear steps); R reaches it
    # via qbeta.c's log-x scale (use_log_x) - the deep-scope class,
    # recorded as a boundary, not fixed ( resolution). R itself
    # is tail-limited in this cell class per the v6 finding.
    # RESOLVED - the direct lane's log-x
    # Newton (this slice's port) reaches the root; the plateau's
    # constant 1.126e15 across 1-1e-11..1-1e-16 is gone, every band
    # cell tracks R at rel <= 2e-13 (cr014-s4 rows carry the band)
    (is (close 4.0530266635530702e+23 (dist/qf (- 1 1e-12) 1 1 o))
        "qf(1-1e-12,1,1) = R (was the 1.126e15 plateau -  resolved the deviation)")))

(deftest mechanical-rows
  # p-adjust method validation fires
  # BEFORE the n<=1 short-circuit (R errors at any length); tukey-hsd
  # validates the required fit keys (no hang on malformed structs)
  (var caught nil)
  (try (tests/p-adjust @[0.1] {:method :banana}) ([e] (set caught e)))
  (is (and caught (string/find "unknown method" caught))
      "p-adjust n=1 bad method -> the VALIDATION shape (was silent)")
  (set caught nil)
  (try (tests/p-adjust @[0.1 0.2] {:method :banana}) ([e] (set caught e)))
  (is (and caught (string/find "unknown method" caught))
      "p-adjust n=2 bad method -> same shape (unchanged)")
  (set caught nil)
  (try (regsurv/tukey-hsd {:mse 1 :df-residual 1 :means @[1 2]} {})
    ([e] (set caught e)))
  (is (and caught (string/find "missing :k" caught))
      "tukey-hsd malformed fit -> VALIDATION, not a hang"))

(deftest arity-rows
  # opts is ALWAYS optional + trailing - the
  # api-contract §1 row is now true as written. The old required-opts
  # tier-A/B/C signatures and the variance mid-signature are gone.
  (is (= 2 (descriptive/mean @[1 2 3])) "mean with NO opts works")
  (is (close 0.9331927987311419 (dist/pnorm 1.5 0 1)) "pnorm no-opts works (R spot)")
  (is (= 1 (dist/pchisq math/inf 5)) "pole row through no-opts lane")
  (is (close 1 (descriptive/variance @[1 2 3])) "variance auto-mean (m optional)")
  (is (= 1 (descriptive/sd @[1 2 3])) "sd auto-mean (sample sd, R spot)")
  (is (zero? (descriptive/skewness @[1 2 3 4 5])) "skewness auto m/s (symmetric = 0)")
  (is (= 1 (descriptive/variance @[1 2 3] 2 {:impl :pure}))
      "variance precomputed-mean fast path intact (new order)")
  (is (= 1 (descriptive/variance @[1 2 3] 2 {:impl :num}))
      "same via :num (scalar policy identical)"))


(deftest exact-tests-contracts
  # p.adjust monotone clamps + ordering; the exact tests' p in
  # [0,1]; CI sanity; sign-test zero-drop conservation
  (def pv @[0.01 0.04 0.03 0.005 0.04 0.2])
  (each method [:bonferroni :holm :hochberg :bh :by]
    (def adj (tests/p-adjust pv {:method method}))
    (is (= (length adj) (length pv)) "length preserved")
    (is (all |(and (>= $ 0) (<= $ 1)) adj) "adjusted p in [0,1]")
    (is (all (fn [i] (<= (pv i) (adj i))) (range 0 6))
        "adjusted >= raw (the monotone clamps)"))
  (each o [{:impl :pure} {:impl :num}]
    (is (deep= (tests/p-adjust pv (merge o {:method :holm}))
               (tests/p-adjust pv (merge o {:method :holm})))
        "p-adjust stable within impl")
    (is (deep= (tests/p-adjust pv (merge {:method :holm} o))
               (tests/p-adjust pv {:method :holm}))
        "p-adjust num == pure"))
  # the error taxonomy
  (var caught nil)
  (try (tests/p-adjust @[-0.1 0.5] {}) ([e] (set caught e)))
  (is (and caught (string/find "p in [0,1]" caught))
      "negative p -> the VALIDATION shape")
  (set caught nil)
  (try (tests/p-adjust @[0.5] {:n 0}) ([e] (set caught e)))
  (is (and caught (string/find "n >= the number" caught))
      "n < lp -> the VALIDATION shape")
  # the exact tests
  (def b (tests/binom-test 7 20 0.3 {}))
  (is (and (>= (b :p-value) 0) (<= (b :p-value) 1)) "binom p in [0,1]")
  (is (and (<= (get (b :conf-interval) 0) (b :estimate))
           (<= (b :estimate) (get (b :conf-interval) 1)))
      "CP CI brackets the estimate")
  (def pt (tests/poisson-test 12 2.5 1.4 {}))
  (is (and (<= (get (pt :conf-interval) 0) (pt :estimate))
           (<= (pt :estimate) (get (pt :conf-interval) 1)))
      "rate CI brackets the estimate")
  (def st (tests/sign-test @[0.5 -0.3 0.2 0.9 -0.1 0.01 -0.02] {}))
  (is (= 0 (st :n-zero)) "no zeros dropped")
  (def st2 (tests/sign-test @[0.3 0.0 -0.2 0.0 0.1 -0.05 0.4] {}))
  (is (= 2 (st2 :n-zero)) "zeros dropped (the R convention)")
  (is (= 5 (st2 :parameter)) "n = the nonzero count after the drop"))

(deftest prop-test-contracts
  # the YATES cap algebra; the null-DEPENDENT cc CI; the k=2
  # given-p lane (df k, forced two-sided, nil CI); both-impl parity;
  # the error taxonomy
  (def r0 (tests/prop-test 5 10 0.5 {}))
  (is (= 0 (r0 :statistic)) "x == n*p: the capped correction is 0")
  (is (= 1 (r0 :p-value)) "exact-null p = 1")
  (def r1 (tests/prop-test 5 10 0.6 {}))
  (is (close (r1 :statistic) 0.1041666666666667)
      "|x - np| = 1: the full 0.5 correction applies")
  (is (not (close (get (r1 :conf-interval) 0) (get (r0 :conf-interval) 0)))
      "the cc CI is null-DEPENDENT (the cap carries n*p)")
  (def r2 (tests/prop-test 0 1 0.5 {}))
  (is (and (= 0 (r2 :statistic)) (= 1 (r2 :p-value)))
      "degenerate n=1, x=0: stat 0, p 1")
  (def r3 (tests/prop-test @[5 5] @[10 10] nil {}))
  (is (= 0 (r3 :statistic)) "k=2 DELTA=0: the correction caps to 0")
  (is (nil? (r3 :null-value)) "the equality lane has no null value")
  (is (deep= @[0.5 0.5] (r3 :estimate)) "k=2 estimate = the pair")
  (def r4 (tests/prop-test @[7 3] @[12 15] @[0.5 0.4] {:alternative :less}))
  (is (= 2 (r4 :parameter)) "k=2 given p: df = k")
  (is (= :two-sided (r4 :alternative)) "k=2 given p: alternative forced two-sided")
  (is (nil? (r4 :conf-interval)) "k=2 given p: no CI (R verbatim)")
  # asymmetric sanity: the CI brackets the estimate / DELTA
  (def r5 (tests/prop-test 7 20 0.3 {}))
  (is (and (<= (get (r5 :conf-interval) 0) (r5 :estimate))
           (<= (r5 :estimate) (get (r5 :conf-interval) 1)))
      "k=1 CI brackets p-hat")
  (is (= 0.35 (r5 :estimate)) "estimate = x/n")
  (def r6 (tests/prop-test @[7 3] @[12 15] nil {}))
  (is (and (<= (get (r6 :conf-interval) 0) (- 0.5833333333333334 0.2))
           (<= (- 0.5833333333333334 0.2) (get (r6 :conf-interval) 1)))
      "k=2 CI brackets DELTA")
  # one-sided CI shapes: less @[0 hi], greater @[lo 1]
  (def rl (tests/prop-test 9 12 0.5 {:alternative :less}))
  (is (and (= 0 (get (rl :conf-interval) 0)) (< (get (rl :conf-interval) 1) 1))
      "less CI = [0, hi]")
  (is (close (rl :p-value) 0.92554266341061719)
      "one-sided less p (sign*sqrt, the erfc route)")
  (def rg (tests/prop-test 9 12 0.5 {:alternative :greater}))
  (is (and (> (get (rg :conf-interval) 0) 0) (= 1 (get (rg :conf-interval) 1)))
      "greater CI = [lo, 1]")
  (is (close (rg :p-value) 0.074457336589382833)
      "one-sided greater p (the erfc far-tail route)")
  # both impls: identical results (scalar composition, the
  # no-spurious-:num-rejection property)
  (each o [{:impl :pure} {:impl :num}]
    (is (deep= (tests/prop-test 7 20 0.3 o)
               (tests/prop-test 7 20 0.3 {}))
        "prop-test identical across impls (k=1)")
    (is (deep= (tests/prop-test @[7 3] @[12 15] nil o)
               (tests/prop-test @[7 3] @[12 15] nil {}))
        "prop-test identical across impls (k=2)"))
  # the error taxonomy (VALIDATION shapes + the COMPOSITION guard)
  (var caught nil)
  (try (tests/prop-test 7 0 0.3 {}) ([e] (set caught e)))
  (is (and caught (string/find "n must be positive" caught))
      "n <= 0 -> the VALIDATION shape")
  (set caught nil)
  (try (tests/prop-test 13 10 0.3 {}) ([e] (set caught e)))
  (is (and caught (string/find "not be greater than n" caught))
      "x > n -> the VALIDATION shape")
  (set caught nil)
  (try (tests/prop-test 5 10 0 {}) ([e] (set caught e)))
  (is (and caught (string/find "p must be in (0,1)" caught))
      "p outside (0,1) -> the VALIDATION shape")
  (set caught nil)
  (try (tests/prop-test @[1 2 3] @[4 5 6] nil {}) ([e] (set caught e)))
  (is (and caught (string/find "1 or 2 proportions" caught))
      "k > 2 -> the charter VALIDATION shape")
  (set caught nil)
  (try (tests/prop-test 5 10 0.3 {:alternative :both}) ([e] (set caught e)))
  (is (and caught (string/find ":two-sided/:less/:greater" caught))
      "bad alternative -> the VALIDATION shape")
  (set caught nil)
  (try (tests/prop-test {:x 5} 10 0.3 {}) ([e] (set caught e)))
  (is (and caught (string/find "unpack packed handles first" caught))
      "a struct x -> the COMPOSITION shape"))

(deftest fixture-dialect-parser-contract
  # fixture-csv-dialect.md (retrofit R6): the JOIN rule - a record is
  # 3 logical fields but args MAY contain commas; expected is the LAST
  # field; header lines (containing `fn,`) and lines with `"` skip;
  # empty lines filter. The twice-shipped bug class, pinned.
  (def text (string "fn,args,expected\n"
                   "foo,1 2,3 4,7\n"
                   "\"quoted\",header,9\n"
                   "\n"
                   "bar,4.5,6"))
  (def rows @[])
  (each line (string/split "\n" text)
    (when (and (not (empty? line)) (not (string/find "fn," line))
               (not (string/find "\"" line)))
      (array/push rows line)))
  (is (= 2 (length rows)) "header, quoted-header, and empty lines skipped")
  (def fields (string/split "," (rows 0)))
  (is (and (= "foo" (fields 0))
           (= "1 2,3 4" (string/join (slice fields 1 -2) ","))
           (= 7 (scan-number (last fields))))
      "the JOIN rule: args reconstructed with commas, expected = LAST")
  (def fields2 (string/split "," (rows 1)))
  (is (and (= "bar" (fields2 0)) (= 6 (scan-number (last fields2))))
      "the simple record parses; last-field rule holds"))

(deftest cran-oracle-contracts
  # epi-2x2 NaN triples + the case-control surface; cohens-d lane
  # semantics + the exact J algebra; the one-sided ncp CI limits;
  # both-impl parity (incl. the resident moment path); error taxonomy
  # epi-2x2: the zero-cell pole NaNs the RR/OR triples (epiR's guard)
  (def z (tests/epi-2x2 @[@[0 5] @[7 20]] {}))
  (each m [(z :rr) (z :or)]
    (is (not= (m :estimate) (m :estimate)) "zero cell: est is NaN")
    (is (all |(not= $ $) (m :conf-interval)) "zero cell: NaN CI triple"))
  (is (= -0.25925925925925924 ((z :rd) :estimate))
      "zero cell: RD stays finite (no cell divisions)")
  # cohort vs case-control surface
  (def cc (tests/epi-2x2 @[@[12 5] @[7 20]] {:method :case-control}))
  (is (and (nil? (cc :rr)) (nil? (cc :rd)) (nil? (cc :rd-score)))
      "case-control: no cohort measures")
  (is (nil? ((tests/epi-2x2 @[@[12 5] @[7 20]] {}) :missing))
      "cohort-count: the cohort measures present")
  # the proportion-scale contract (epiR displays x100; ours are raw)
  (def r-epi (tests/epi-2x2 @[@[12 5] @[7 20]] {}))
  (is (= 0.44662309368191727 ((r-epi :rd) :estimate))
      "RD = a/N1 - c/N0 (proportions, not the x100 display scale)")
  # cohens-d lane semantics (df per the decoded .effect_size_difference)
  (def d-pool (tests/cohens-d @[1 2 3 4 5 6] @[3 4 5 6 7 8] {}))
  (is (= 10 (d-pool :df)) "pooled: df = n1+n2-2")
  (is (close 10 ((tests/cohens-d @[1 2 3 4 5 6] @[3 4 5 6 7 8]
                            {:pooled-sd false}) :df))
      "equal spreads: Welch df coincides with pooled here (the
      float sum lands an ulp off the integer - close, not =)")
  (def d-unp (tests/cohens-d @[2.5 3.1 2.8 3.4 2.9 3.2] @[5.1 4.8 5.4 5.0 4.9]
                             {:pooled-sd false}))
  (is (not= 9 (d-unp :df)) "unequal spreads: Welch df differs from pooled")
  (is (= 4 ((tests/cohens-d @[2 4 6 8 10] nil {}) :df))
      "one-sample: df = n-1")
  (is (= 4 ((tests/cohens-d @[1 2 3 4 5] @[2 3 3 5 6] {:paired true}) :df))
      "paired: df = pairs-1")
  # hedges-g: the EXACT J algebra (g = J*d, bounds scaled the same)
  (def g-pool (tests/hedges-g @[1 2 3 4 5 6] @[3 4 5 6 7 8] {}))
  (def J10 (math/exp (- (utils/lgamma 5) (math/log (math/sqrt 5)) (utils/lgamma 4.5))))
  (is (close (* J10 (d-pool :estimate)) (g-pool :estimate))
      "g = J(d) * d (the exact gamma form)")
  (is (close (* J10 (get (d-pool :conf-interval) 0))
             (get (g-pool :conf-interval) 0))
      "the ncp CI bounds scale by the same J")
  # the one-sided CI limits (.adjust_ci + .limit_ci)
  (def d-less (tests/cohens-d @[1 2 3 4 5 6] @[3 4 5 6 7 8]
                              {:alternative :less}))
  (is (and (= (- math/inf) (get (d-less :conf-interval) 0))
           (< (get (d-less :conf-interval) 1) 0))
          "less: [-inf, hi] at the 2*conf-1 level")
  (def d-gr (tests/cohens-d @[1 2 3 4 5 6] @[3 4 5 6 7 8]
                            {:alternative :greater}))
  (is (and (<= (math/abs (- (get (d-gr :conf-interval) 0) -2.0708442740711948))
               1e-6)
           (= math/inf (get (d-gr :conf-interval) 1)))
          "greater: [lo, +inf] (lo vs the effectsize row at its
          optim-noise tolerance; ours is the exact solve)")
  # the pnt route (dist): monotone in ncp; spot-checked vs R's pt
  (is (< (dist/pt 1 10 {:ncp 2}) (dist/pt 1 10 {:ncp 1}))
      "pt :ncp is monotone decreasing in ncp")
  (is (close (dist/pt 1 10 {:ncp 2}) 0.15851505953749501)
      "pt :ncp spot value (R pt(1, 10, 2))")
  # both impls: identical results (the resident moment path is
  # bitwise; the ncp series is scalar)
  (each o [{:impl :pure} {:impl :num}]
    (is (deep= (tests/cohens-d @[1 2 3 4 5 6] @[3 4 5 6 7 8] o)
               (tests/cohens-d @[1 2 3 4 5 6] @[3 4 5 6 7 8] {}))
        "cohens-d identical across impls")
    (is (deep= (tests/cohens-d @[1 2 3 4 5] @[2 3 3 5 6] (merge {:paired true} o))
               (tests/cohens-d @[1 2 3 4 5] @[2 3 3 5 6] {:paired true}))
        "the paired vmap2 diff path is bitwise across impls")
    (is (deep= (tests/epi-2x2 @[@[12 5] @[7 20]] o)
               (tests/epi-2x2 @[@[12 5] @[7 20]] {}))
        "epi-2x2 identical across impls"))
  # the error taxonomy
  (var caught nil)
  (try (tests/cohens-d @[1 1 1] @[2 2 2] {}) ([e] (set caught e)))
  (is (and caught (string/find "nonzero sample sd" caught))
      "constant input -> the VALIDATION shape")
  (set caught nil)
  (try (tests/cohens-d @[1] @[2] {}) ([e] (set caught e)))
  (is (and caught (string/find ">= 2 observations" caught))
      "n < 2 -> the VALIDATION shape")
  (set caught nil)
  (try (tests/cohens-d @[1 2 3] @[2 3 4] {:alternative :both}) ([e] (set caught e)))
  (is (and caught (string/find ":two-sided/:less/:greater" caught))
      "bad alternative -> the VALIDATION shape")
  (set caught nil)
  (try (tests/cohens-d @[1 2 3] nil {:paired true}) ([e] (set caught e)))
  (is (and caught (string/find "needs a second sample" caught))
      "paired without ys -> the VALIDATION shape")
  (set caught nil)
  (try (tests/epi-2x2 @[@[12 5] @[7 20]] {:method :cohort-time}) ([e] (set caught e)))
  (is (and caught (string/find ":cohort-count/:case-control" caught))
      "bad method -> the VALIDATION shape")
  (set caught nil)
  (try (tests/epi-2x2 @[@[12 5] @[7 -20]] {}) ([e] (set caught e)))
  (is (and caught (string/find "nonnegative" caught))
      "negative cell -> the VALIDATION shape"))

(deftest cor-test-contracts
  # the exact/approx switch decode; the |r|=1 and
  # NaN poles; n<4 lanes; both-impl parity (pearson resident);
  # the error taxonomy. The DP spot values are R-bitwise.
  (def a @[2.1 3.4 1.8 2.9 2.5 3.9 1.1 2.2])
  (def b @[3.1 6.2 2.8 8.4 4.9 7.7 1.2 5.5])
  # the exact DPs: bitwise vs R's C_pKendall / C_pRho
  (is (= 0.86388888888888893 (ranksum/kendall-p 10 6))
      "kendall-p spot (R pkendall(10, 6))")
  (is (= 0.041666666666666664 (ranksum/spearman-p 2 4 true))
      "spearman-p spot (R prho(2, 4, lower))")
  (is (close 0.0009463183422
             ((tests/cor-test (range 1 11) @[2 1 4 3 6 5 8 7 10 9]
                              {:method :kendall}) :p-value))
      "the ckendall DP through cor-test (R exact)")
  # the switch decode: exact vs approx differ observably
  (is (not (close ((tests/cor-test a b {:method :kendall}) :p-value)
                  ((tests/cor-test a b {:method :kendall :exact false}) :p-value)))
      "kendall exact != S-approx (the switch is real)")
  (is (not (close ((tests/cor-test a b {:method :spearman}) :p-value)
                  ((tests/cor-test a b {:method :spearman :exact false}) :p-value)))
      "spearman brute != t-approx (the switch is real)")
  # tied input: the auto-approx fires (z finite, not the raw T)
  (def tz ((tests/cor-test @[1 1 2 2 3 3 4 5 6] @[2 3 1 4 2 5 3 6 5]
                           {:method :kendall}) :statistic))
  (is (close 1.7228220865861292 tz)
      "tied kendall z matches R's S-approx (the auto-approx fired)")
  # the |r| = 1 pole: t = inf, p = 0, CI degenerate (n > 3 for the CI)
  (def r1 (tests/cor-test @[1 2 3 4 5] @[2 4 6 8 10] {}))
  (is (= math/inf (r1 :statistic)) "|r|=1: t = +inf")
  (is (= 0 (r1 :p-value)) "|r|=1: two-sided p = 0")
  (is (deep= @[1 1] (r1 :conf-interval)) "|r|=1: the Fisher-z CI collapses")
  # constant input mirrors R's NA (NaN triple)
  (def rc (tests/cor-test @[1 1 2 2] @[3 3 3 3] {}))
  (each k [:estimate :statistic :p-value]
    (is (not= (rc k) (rc k)) (string "constant y: NaN " k)))
  # n=2 kendall tied: R's 0/0 NaN lane
  (def r2t (tests/cor-test @[1 1] @[1 2] {:method :kendall}))
  (is (not= (r2t :p-value) (r2t :p-value))
      "n=2 tied kendall: NaN (R's v2 term 0/0)")
  # n=3 pearson: no CI; n=2 pearson: error; n=2 kendall: fine
  (is (nil? ((tests/cor-test @[1 2 9] @[2 4 5] {}) :conf-interval))
      "n=3 pearson: no CI (R computes it only for n > 3)")
  (var caught nil)
  (try (tests/cor-test @[1 2] @[2 4] {}) ([e] (set caught e)))
  (is (and caught (string/find ">= 3 observations" caught))
      "pearson n < 3 -> the VALIDATION shape")
  (set caught nil)
  (try (tests/cor-test a b {:method :rho}) ([e] (set caught e)))
  (is (and caught (string/find ":pearson/:spearman/:kendall" caught))
      "bad method -> the VALIDATION shape")
  (set caught nil)
  (try (tests/cor-test a b {:conf-level 1.5}) ([e] (set caught e)))
  (is (and caught (string/find "in [0,1]" caught))
      "conf-level out of range -> the VALIDATION shape")
  (set caught nil)
  (try (tests/cor-test a @[1 2] {}) ([e] (set caught e)))
  (is (and caught (string/find "same length" caught))
      "length mismatch -> the VALIDATION shape")
  (set caught nil)
  (try (tests/cor-test @[1 2 math/nan] @[1 2 3] {}) ([e] (set caught e)))
  (is (and caught (string/find "finite observations" caught))
      "non-finite input -> the VALIDATION shape (ours; R drops NAs)")
  (set caught nil)
  (try (tests/cor-test (num/stage @[1 2 3 4 5 6 7 8] (num/kernels)) b
                       {:method :spearman})
    ([e] (set caught e)))
  (is (and caught (string/find "unpack packed handles first" caught))
      "handle into the sort-bound lane -> the COMPOSITION shape")
  # both impls: pearson is RESIDENT (handles chain); identical results
  (each o [{:impl :pure} {:impl :num}]
    (is (deep= (tests/cor-test a b o) (tests/cor-test a b {}))
        "cor-test identical across impls (pearson)"))
  (def r-pure (tests/cor-test a b {:impl :pure}))
  (def r-num (tests/cor-test a b {:impl :num}))
  (is (close (r-pure :p-value) (r-num :p-value))
      "pearson p identical across impls (scalar tails)"))

(deftest mcnemar-mh-contracts
  # the cc gates (mcnemar's b!=c condition;
  # mantelhaen's |DELTA| >= 0.5 cap); the NaN poles; the homog
  # surface; both-impl parity; the error taxonomy
  # mcnemar: the correction applies ONLY when the discordant cells
  # differ - equal margins get stat 0 even with :correct
  (def eq (tests/mcnemar-test @[@[15 5] @[5 2]] {}))
  (is (and (= 0 (eq :statistic)) (= 1 (eq :p-value))
           (= 0 ((tests/mcnemar-test @[@[15 5] @[5 2]] {:correct false})
                 :statistic)))
      "n12 == n21: stat 0 both cc settings (the cc gate)")
  (is (= 0.8421052631578947
         ((tests/mcnemar-test @[@[15 12] @[7 2]] {}) :statistic))
      "b != c: the cc subtracts 1 (R's convention)")
  (def zd (tests/mcnemar-test @[@[15 0] @[0 2]] {}))
  (each k [:statistic :p-value]
    (is (not= (zd k) (zd k)) (string "zero discordant pairs: NaN " k)))
  # mantelhaen: the |DELTA| >= 0.5 cap - DELTA 0 gets NO correction
  (def zt (tests/mantelhaen-test @[@[@[5 5] @[5 5]]] {}))
  (is (= 0 (zt :statistic))
      "DELTA = 0: the cap does not fire (stat 0 with :correct on)")
  (is (and (nil? (zt :homogeneity)) (not (nil? (zt :estimate))))
      "K = 1: the homog keys are nil, the base surface complete")
  # the sparse pole: all bc-zero strata -> OR inf (R's x/0)
  (def sp (tests/mantelhaen-test @[@[@[3 0] @[0 3]] @[@[2 0] @[0 2]]] {}))
  (is (= math/inf (sp :estimate)) "all-BC-zero: the common OR is inf")
  # both impls: identical results (scalar compositions)
  (def t3 @[@[@[12 5] @[7 20]] @[@[8 3] @[6 10]] @[@[20 9] @[14 30]]])
  (each o [{:impl :pure} {:impl :num}]
    (is (deep= (tests/mantelhaen-test t3 o) (tests/mantelhaen-test t3 {}))
        "mantelhaen identical across impls")
    (is (deep= (tests/mcnemar-test @[@[15 12] @[7 2]] o)
               (tests/mcnemar-test @[@[15 12] @[7 2]] {}))
        "mcnemar identical across impls"))
  # the error taxonomy
  (var caught nil)
  (try (tests/mcnemar-test @[@[1 2 3] @[4 5]] {}) ([e] (set caught e)))
  (is (and caught (string/find "2x2 table" caught))
      "non-square mcnemar input -> the VALIDATION shape")
  (set caught nil)
  (try (tests/mcnemar-test @[@[1 -2] @[3 4]] {}) ([e] (set caught e)))
  (is (and caught (string/find "nonnegative" caught))
      "negative cell -> the VALIDATION shape")
  (set caught nil)
  (try (tests/mantelhaen-test @[@[@[1 0] @[0 0]]] {}) ([e] (set caught e)))
  (is (and caught (string/find "> 1" caught))
      "stratum n < 2 -> the VALIDATION shape")
  (set caught nil)
  (try (tests/mantelhaen-test t3 {:alternative :both}) ([e] (set caught e)))
  (is (and caught (string/find ":two-sided/:less/:greater" caught))
      "bad alternative -> the VALIDATION shape"))

(deftest ptukey-pure-contracts
  # the pure studentized range - R's
  # ptukey.c/qtukey.c verbatim (Copenhaver-Holland + AS 70). The
  # R-pinned spot values; the roundtrip; monotonicity; the
  # nranges effect; the poles; both-impl parity + the vector shape
  # (the C kernel + :num vector lane land at )
  (each [q k df ev] [[2.5 4 10 0.65804329973306763]
                     [3.9 5 24 0.92570794527737921]
                     [1.5 8 300 0.035746771957133584]
                     [4.6 3 2000 0.99666941071606585]
                     [0.5 2 9000 0.27631810893711167]
                     [2.5 2 30000 0.9229001282564584]]
    (is (close (dist/ptukey q k df {}) ev)
        (string "ptukey(" q "," k "," df ") matches R")))
  (each [p k df ev] [[0.05 3 4 0.43420371178581874]
                     [0.95 5 10 4.6542930755433387]
                     [0.99 2 120 3.7015924818950512]]
    (is (<= (math/abs (- (dist/qtukey p k df {}) ev)) 1e-6)
        (string "qtukey(" p "," k "," df ") matches R (secant eps 1e-4)")))
  # the roundtrip (bounded by the secant's eps 1e-4: R's own qtukey
  # stops there - the CDF residual ~ slope * residual quantile)
  (is (<= (math/abs (- (dist/ptukey (dist/qtukey 0.95 4 20 {}) 4 20 {}) 0.95))
          1e-5)
      "ptukey(qtukey(p)) = p (within the secant's eps)")
  # monotone in q; nranges shrinks the CDF (more ranges -> smaller)
  (is (< (dist/ptukey 2 4 10 {}) (dist/ptukey 2.5 4 10 {}))
      "ptukey monotone in q")
  (is (< (dist/ptukey 2.5 3 10 {:nranges 2}) (dist/ptukey 2.5 3 10 {}))
      "nranges > 1: more ranges -> smaller CDF")
  # the poles + boundaries
  (is (= 0 (dist/ptukey 0 3 10 {})) "q = 0 -> 0")
  (is (= 0 (dist/qtukey 0 3 10 {})) "p = 0 -> 0")
  (is (= math/inf (dist/qtukey 1 3 10 {})) "p = 1 -> +inf")
  # the vector shape: the DEFAULT lane now routes the ptukey-v
  # kernel - elementwise agreement within the ABI bar
  (def qv @[2.5 3.9 4.6])
  (def qv-num (dist/ptukey qv 4 20 {}))
  (def qv-pure (seq [q :in qv] (dist/ptukey q 4 20 {})))
  (var qv-worst 0)
  (for i 0 3 (set qv-worst (max qv-worst (math/abs (- (qv-num i) (qv-pure i))))))
  (is (<= qv-worst 1e-9)
      "the vector q lane: default :num == elementwise pure (ABI bar)")
  (each o [{:impl :pure} {:impl :num}]
    (is (close (dist/ptukey 2.5 4 10 o) (dist/ptukey 2.5 4 10 {}))
        "ptukey identical across impls (identical code)"))
  (var caught nil)
  (try (dist/ptukey 2.5 1 10 {}) ([e] (set caught e)))
  (is (and caught (string/find "nmeans >= 2" caught))
      "nmeans < 2 -> the VALIDATION shape (ours; R NaNs)")
  (set caught nil)
  (try (dist/ptukey 2.5 4 1 {}) ([e] (set caught e)))
  (is (and caught (string/find "df >= 2" caught))
      "df < 2 -> the VALIDATION shape")
  # out-of-domain p -> NaN (the old
  # throw withdrawn - the q-pole convention, family-wide)
  (is (not= (dist/qtukey 1.5 3 10 {}) (dist/qtukey 1.5 3 10 {}))
      "p out of range -> NaN (was the VALIDATION throw)"))

(deftest ptukey-num-contracts
  # the ptukey-v kernel lane - the ABI
  # parity bar 1e-9; resident handles chain; scalar :num = the pure
  # formula; the NaN pole divergence documented (the pure carries
  # R's ISNAN gate; the kernel's pole change rides its own CR)
  (def qs @[0.3 1.1 2.5 3.9 4.6])
  (def pure (map |(dist/ptukey $ 5 32 {:impl :pure}) qs))
  (def num (dist/ptukey qs 5 32 {:impl :num}))
  (var worst 0)
  (for i 0 (length qs)
    (set worst (max worst (math/abs (- (num i) (pure i))))))
  (is (<= worst 1e-9)
      (string "vector :num == pure within the ABI bar (worst " worst ")"))
  (is (all |(number? $) num) "the kernel lane unpacks (no handle leaks)")
  (is (deep= (dist/ptukey (num/stage qs (num/kernels)) 5 32 {:impl :num})
             num)
      "packed handles chain (RESIDENT)")
  (is (close (dist/ptukey 2.5 4 10 {:impl :num}) (dist/ptukey 2.5 4 10 {}))
      "scalar :num = the pure formula")
  # the NaN pole:  LANDED (janet-num 1.2.0+, the 4-arg
  # ISNAN gate returning the caller's NaN) - pure and kernel now
  # AGREE with R (NaN -> NaN); the pre-CR divergence is unpinned
  (is (not= (dist/ptukey math/nan 4 10 {}) (dist/ptukey math/nan 4 10 {}))
      "pure NaN q -> NaN (R's ISNAN gate - the user decision)")
  (def knan (first (dist/ptukey @[math/nan] 4 10 {:impl :num})))
  (is (not= knan knan)
      "kernel NaN q -> NaN (the  pole gate; pure/kernel agreement)"))

(deftest aov-contracts
  # the ANOVA table
  # contracts, the Welch/classic oneway lanes, the Tukey pair
  # convention + poles (janet division throws - every guard explicit,
  # R's Inf/NaN mirrored), both-impl parity. R VALUE conformance
  # lives in the aov harness rows (294 rows, both impls).
  (def balanced @[@[1 2 3 4 5] @[2 3 4 5 6] @[5 6 7 8 9]])
  (def fit (regsurv/aov balanced {:impl :pure}))
  (def tg (get (get fit :table) :group))
  (def tr (get (get fit :table) :residual))
  # the SS identity: SSB + SSW == TSS; table F == lm's overall F
  (var tss 0)
  (def ys (array ;(mapcat identity balanced)))
  (def ybar (descriptive/mean ys {}))
  (each g balanced (each v g (+= tss (* (- v ybar) (- v ybar)))))
  (is (close (+ (tg :ss) (tr :ss)) tss) "SSB + SSW == TSS (the orthogonal decomposition)")
  (is (close (tg :f) (fit :f)) "table F == lm's overall F (the one-way identity)")
  (is (and (>= (tg :f) 0) (<= 0 (tg :p) 1)) "F >= 0, p in [0,1]")
  (is (= (tg :df) 2) "group df = k-1")
  (is (= (tr :df) 12) "residual df = n-k")
  # the Welch/classic divergence (unequal spreads move the lanes)
  (def spread @[@[1 2 3 4 5] @[10 20 30 40 50] @[0.5 1 1.5 2 2.5]])
  (def w (tests/oneway-test spread))
  (def c (tests/oneway-test spread {:var-equal true}))
  (is (and (not= (w :statistic) (c :statistic))
           (not= (get (w :parameter) 1) (get (c :parameter) 1)))
      "welch F/df2 differ from the classic lane on unequal spreads")
  (is (= (w :var-equal) false) "the DEFAULT is the welch lane (R's default)")
  # the pair convention: column-major lower-tri, later-minus-earlier
  (def t (regsurv/tukey-hsd fit {:impl :pure}))
  (is (deep= (map |(get $ :label) (t :pairs)) @["2-1" "3-1" "3-2"])
      "pair order: 2-1, 3-1, 3-2 (R's column-major lower-tri)")
  (is (close ((first (t :pairs)) :diff) 1) "diff = later minus earlier level")
  # conf shifts the CI width, never p
  (def t90 (regsurv/tukey-hsd fit {:conf-level 0.9 :impl :pure}))
  (is (= ((first (t90 :pairs)) :p) ((first (t :pairs)) :p))
      "conf.level moves the width only (p is conf-independent)")
  (is (< ((first (t :pairs)) :lwr) ((first (t90 :pairs)) :lwr))
      "the 0.95 CI is wider than the 0.9 CI")
  # :ordered sorts by means; :i/:j stay the ORIGINAL indices
  (def two @[@[5 6 7] @[1 2 3]])
  (def to (regsurv/tukey-hsd (regsurv/aov two) {:ordered true}))
  (is (= ((first (to :pairs)) :label) "1-2")
      "ordered: the later-sorted level minus the earlier (label 1-2)")
  (is (and (= ((first (to :pairs)) :i) 0) (= ((first (to :pairs)) :j) 1))
      "ordered pairs carry ORIGINAL indices")
  # both-impl parity: the lm P1.3 parity carries the table; the tukey
  # p's route the ptukey vector lane (num) vs the scalar loop (pure)
  (def fit-n (regsurv/aov balanced {:impl :num}))
  (is (close ((get (get fit-n :table) :group) :f) (tg :f))
      "aov num == pure (the lm P1.3 parity carries the table)")
  (def t-n (regsurv/tukey-hsd fit-n {:impl :num}))
  (var worst 0)
  (each [pn pp] (map (fn [a b] [a b]) (t-n :pairs) (t :pairs))
    (set worst (max worst (math/abs (- (pn :p) (pp :p))))))
  (is (<= worst 1e-9)
      (string "tukey p num == pure within the ABI bar (worst " worst ")"))
  # poles: k = 1, n_i = 1 (oneway), df.res = 1 (tukey), non-finite
  (var caught nil)
  (try (regsurv/aov @[@[1 2 3]]) ([e] (set caught e)))
  (is caught "aov k = 1 -> VALIDATION (R: not enough groups)")
  (set caught nil)
  (try (tests/oneway-test @[@[1 2 3] @[4]]) ([e] (set caught e)))
  (is caught "oneway n_i = 1 -> VALIDATION (R: not enough observations)")
  (set caught nil)
  (try (tests/oneway-test @[@[1 2 3] @[4 math/nan]]) ([e] (set caught e)))
  (is caught "non-finite input -> VALIDATION")
  (set caught nil)
  (try (regsurv/aov @[@[1 2] @[3 4] "ab"]) ([e] (set caught e)))
  (is caught "non-indexed group -> VALIDATION")
  # singleton groups are LEGAL in aov (unlike oneway.test); but
  # n = k+1 leaves df.res = 1 -> tukey-hsd throws ptukey's df >= 2
  (def tiny (regsurv/aov @[@[1.0] @[2 3]]))
  (is (= (tiny :k) 2) "singleton group legal in aov")
  (set caught nil)
  (try (regsurv/tukey-hsd tiny {}) ([e] (set caught e)))
  (is caught "df.res = 1 -> tukey-hsd errors (ptukey's df >= 2 VALIDATION)")
  # MSE = 0, unequal means: F = inf, p = 0, est = ±inf -> p 0
  (def const-fit (regsurv/aov @[@[2 2 2] @[5 5 5]] {:impl :pure}))
  (def ctg (get (get const-fit :table) :group))
  (is (= math/inf (ctg :f)) "constant groups, unequal means: F = inf (R's Inf)")
  (is (= 0 (ctg :p)) "pf(inf) upper -> p 0")
  (is (= 0 ((first ((regsurv/tukey-hsd const-fit {}) :pairs)) :p))
      "tukey est = ±inf -> p 0 (both impls)")
  # all-identical y: 0/0 - F NaN, and the ptukey NaN pole rides
  # through as AGREEMENT post-p NaN in BOTH lanes (the
  # kernel's ISNAN gate matches the pure's; R's own value)
  (def same-fit (regsurv/aov @[@[3 3 3] @[3 3 3]] {:impl :pure}))
  (is (not= ((get (get same-fit :table) :group) :f)
            ((get (get same-fit :table) :group) :f))
      "all-identical: F = NaN (R's 0/0)")
  (def same-pure ((first ((regsurv/tukey-hsd same-fit {:impl :pure}) :pairs)) :p))
  (is (not= same-pure same-pure) "pure: NaN est -> NaN p (rides R)")
  (def same-num ((first ((regsurv/tukey-hsd same-fit {:impl :num}) :pairs)) :p))
  (is (not= same-num same-num)
      "kernel: NaN est -> NaN p (the  gate; both lanes agree)")
  # a zero-variance group in the welch lane: w = inf -> R's NaN
  (def zw (tests/oneway-test @[@[1 1 1] @[2 4 6]]))
  (is (and (not= (zw :statistic) (zw :statistic))
           (not= (zw :p-value) (zw :p-value)))
      "welch zero-variance group: F/p NaN (R's Inf-weighted NaN)"))

(deftest packed-gate-rows
  # the widened gate
  # (:ptr REQUIRED, :close clause deleted). The crash-class cell is
  # asserted from OUTSIDE the process - a janet child runs the
  # close-only marker through dnorm (pre-fix: SEGFAULT exit 139,
  # reproduced twice; in-process would kill the suite) and
  # the harness checks exit != 139 AND the scoped caller-misuse error.
  # The junk-:ptr cell is the DEFINED behaviour (the contract row):
  # caller-supplied bytes unpack as data - the return IS dnorm of the
  # bytes' IEEE interpretation (0x4141414141414141 = 2261634.5098039214;
  # its density underflows to exactly 0, the log lane pins the decode).
  (let [cell "(import stats/dist) (try (do (printf \"GOT:%q\" (dist/dnorm (struct :n 1 :bytes 8 :close (fn [_] nil)) 0 1)) (print \"DONE\")) ([e] (print \"CAUGHT:\") (print e)))"
        out (string "/tmp/cr016-sub-" (gensym) "-" (math/floor (* 1000 (os/clock))))
        _ (sh (string "janet -e " (sh-quote cell) " > " (sh-quote out) " 2>&1; echo \"EXIT:$?\" >> " (sh-quote out)))
        txt (string (slurp-string out))
        _ (os/rm out)
        lines (filter |(not (empty? $)) (string/split "\n" txt))
        tail (last lines)
        exit (scan-number (string/replace "EXIT:" "" tail))]
    (is (not= 139 exit)
        "close-only marker: child exit != 139 (pre-fix: SEGFAULT)")
    (is (string/find "caller misuse" txt)
        "close-only marker: the scoped caller-misuse error (subprocess, both impls share the gate)"))
  # in-process rows on the :pure impl only - the default (:num) lane
  # is the subprocess cell above (pre-fix it segfaults in-process)
  (var caught nil)
  (try (dist/dnorm (struct :n 1 :bytes 8 :close (fn [_] nil)) 0 1 {:impl :pure})
    ([e] (set caught e)))
  (is (and caught (string/find "caller misuse" (string caught)))
      "close-only marker, :pure impl -> the same scoped misuse error")
  # junk :ptr bytes = data: 8x 0x41 decodes to 2261634.5098039214
  (let [h (struct :n 1 :bytes 8 :ptr @"AAAAAAAA")
        junk (dist/dnorm h 0 1)
        junk-log (dist/dnorm h 0 1 {:log true})
        ieee 2261634.5098039214]
    (is (= 1 (length junk)) "junk :ptr: one element (:n 1)")
    (is (= 0 (first junk))
        "junk :ptr: dnorm(2.26e6) underflows to exactly 0 (the pinned cell)")
    (is (close (first junk-log) (dist/dnorm ieee 0 1 {:log true}))
        "junk :ptr: log lane = dnorm of the bytes' IEEE interpretation (2261634.5098039214)"))
  # correct bytes through a real handle: unchanged, already green
  # (the  row above: ppois stages a real packed handle)
  (let [h (((num/kernels) 'pack-f64) @[0.5 1.0])
        r (dist/dnorm h 0 1)]
    (is (and (= 2 (length r))
             (close (get r 0) (dist/dnorm 0.5 0 1))
             (close (get r 1) (dist/dnorm 1.0 0 1)))
        "real handle (correct bytes) -> correct numbers (gate does not reject)")))

(deftest buffer-sufficiency-rows
  # shape is necessary
  # but NOT sufficient - the kernel reads :n x :bytes unconditionally,
  # so a short :ptr is an OOB read (silent junk for small overreads,
  # SEGFAULT exit 139 for the large one). All three rows are
  # subprocess-asserted: exit != 139 AND the scoped caller-misuse
  # error. The green row pins zero false positives: a REAL pack-f64
  # handle (:bytes == (length :ptr) by construction) still passes.
  (defn sub [cell]
    (let [out (string "/tmp/cr017-sub-" (gensym) "-" (math/floor (* 1000 (os/clock))))
          _ (sh (string "janet -e " (sh-quote cell) " > " (sh-quote out) " 2>&1; echo \"EXIT:$?\" >> " (sh-quote out)))
          txt (string (slurp-string out))
          _ (os/rm out)
          lines (filter |(not (empty? $)) (string/split "\n" txt))
          exit (scan-number (string/replace "EXIT:" "" (last lines)))]
      [exit txt]))
  (let [[exit txt] (sub "(import stats/dist) (try (do (printf \"GOT:%q\" (dist/dnorm (struct :n 1000000 :bytes 8000000 :ptr @\"AAAAAAAA\") 0 1)) (print \"DONE\")) ([e] (print \"CAUGHT:\") (print e)))")]
    (is (not= 139 exit)
        "n 1000000/bytes 8M/8-byte ptr: child exit != 139 (the pre-fix blocking cell: SEGFAULT)")
    (is (string/find "caller misuse" txt)
        "n 1000000 short ptr: the scoped caller-misuse error (subprocess)"))
  (let [[exit txt] (sub "(import stats/dist) (try (do (printf \"GOT:%q\" (dist/dnorm (struct :n 2 :bytes 16 :ptr @\"AAAAAAAA\") 0 1)) (print \"DONE\")) ([e] (print \"CAUGHT:\") (print e)))")]
    (is (not= 139 exit) "n 2/bytes 16/8-byte ptr: child exit != 139")
    (is (string/find "caller misuse" txt)
        "n 2/bytes 16 short ptr: the scoped caller-misuse error (pre-fix: silent 2-value overread, exit 0)"))
  (let [[exit txt] (sub "(import stats/dist) (try (do (printf \"GOT:%q\" (dist/dnorm (struct :n 1 :bytes 8 :ptr @\"AAA\") 0 1)) (print \"DONE\")) ([e] (print \"CAUGHT:\") (print e)))")]
    (is (not= 139 exit) "n 1/bytes 8/3-byte ptr: child exit != 139")
    (is (string/find "caller misuse" txt)
        "n 1/bytes 8 short ptr: the scoped caller-misuse error (pre-fix: silent 1-value overread, exit 0)"))
  # the gate widened to TABLE-shaped
  # markers - dict? at the guard sites; a table carrying handle keys
  # gets the SAME scoped error (pre-fix: unscoped kernel error, exit 0)
  (var tbl-caught nil)
  (try (dist/dnorm @{:n 1 :bytes 8 :ptr @"AAAAAAAA"} 0 1)
    ([e] (set tbl-caught e)))
  (is (and tbl-caught (string/find "caller misuse" (string tbl-caught)))
      "table-shaped marker: the scoped caller-misuse error (not the unscoped kernel error)")
  (let [h (((num/kernels) 'pack-f64) @[0.5 1.0])
        r (dist/dnorm h 0 1)]
    (is (and (= 2 (length r))
             (close (get r 0) (dist/dnorm 0.5 0 1))
             (close (get r 1) (dist/dnorm 1.0 0 1)))
        "real pack-f64 handle (sufficient by construction) -> still passes, zero false positives")))

(def reports (run-tests! :silent? true :no-exit? true))
(var failures 0)
(each r reports
  (each f (r :failures)
    (++ failures)
    (printf "FAIL %s: %s expect=%q actual=%q"
            (r :test) (f :note) (f :expect) (f :actual))))
(printf "stats suite: %d tests, %d failures" (length reports) failures)
(when (pos? failures) (os/exit 1))
(print "STATS SUITE GREEN")
