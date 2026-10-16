# janet-stats - public API contract (non-Quint, v0.1.0)

The janet-stats analogue of janet-num's kernel-abi.md. Canonical
signatures live in the modules' docstrings; this document states what
signatures cannot carry: the calling conventions, the input-dispatch
taxonomy, and the ERROR TAXONOMY (four classes - the fuzz and smokes
assert message shapes against it).

## 1. Call shape

Every public fn: `(fn data... &opt opts)` - data first, a trailing
opts TABLE last (never a trailing option value). opts is always
optional with documented defaults; every key is a keyword.
UNIFIED (2026-10-08): the older surfaces
that took opts as a REQUIRED positional arg and the variance/sd/
skewness/kurtosis mid-signature `[xs opts m(,s)]` are gone - all
public fns now spell this row as written; the moment fns take the
precomputed moments as optional leading params with auto-computed
defaults (variance `[xs &opt m opts]`, R's var shape). The
resolution key `:impl` applies to ALL public fns (dual-impl, no
exceptions): `{:impl :pure|:num|:impl}` - the default `:impl` is
machine-dynamic (`:num` with libjn, `:pure` without; FORMAL SPEC:
the janet-num repository's specs/impl-resolution.qnt - probe-once,
strict `:num`, `:pure` never loads).

### The opts table (cross-fn keys)

| key | fns | default | meaning |
|---|---|---|---|
| `:impl` | ALL | `:impl` | implementation resolution |
| `:lower-tail` | the thirteen q families (qnorm, qlnorm, qexp, qweibull, qgamma, qbeta, qchisq, qt, qf, qbinom, qpois, qnbinom, qtukey) | true | R's `lower.tail`: false computes the UPPER-tail quantile. Family-wide (was silently ignored by ten of thirteen). Decoded per R: dpq.h's `0.5-p+0.5` complement / raw-side `-log(p)` lanes (qexp, qweibull, qlnorm, qnorm) and qt.c's neg/P two-sided mass - subnormal p survives on those five + the discrete walks' `p*(1+8eps)` fuzz (qDiscrete_search.h). The qgamma/qchisq/qf's qgamma-composed lanes answer the RAW upper p (the qbeta CIv treatment, dpq.h; R's qgamma.c END-block log-p Newton) - the old swallow fired at p < ~1.1e-16 (the materialized 1-p == 1.0, ~16 orders above the documented subnormal boundary); measured: upper-lane p from 1 down to 1e-320 R-exact (1e-300 exact, 1e-320 rel 2e-13); the remaining boundary is denormal p < ~1e-321 where gamma-q's linear-scale continued fraction underflows (R's pgamma computes the log scale internally) - 5e-324 returns 0 vs R 767.8. The log-x lane fixed the qbeta/qf subnormal family too (measured 2026-10-13: qbeta(5e-324,3,7,up)=1 = R exact; qf(5e-324,5,5,up)=4.1350967057808038e+129 = R exact; qf(5e-324,1,1) and qf(1e-320,2,2) unflagged) - only the qgamma/qchisq 5e-324 → 0 boundary above remains. The old text ("p < 2.2e-308 swallows to 1.0") was measured-wrong on the swallow threshold AND is now superseded on the qgamma families. DEVIATION (cells-carrying, measured 2026-10-13): in the upper-lane band p ~ (1e-14, 1e-12), R's own qgamma answers carry materialized-p solver noise (R's pgamma at R's own qgamma(1e-14,5,1,up) is off by dlog -2.99e-9; -2.2e-9 at 1e-13, -3.9e-11 at 5e-13, -3e-12 at 1e-12, converged at 1e-16); our raw-v Newton converges to the TRUE quantile (dlog 0 by R's own pgamma), so ours differs from R's printed value by rel ≤ 1.5e-9 there - truer than the reference, pinned at rel 1e-12 of the true value in the far-tail suite rows |
| `:log-p` | the thirteen q families | REJECTED | R's `log.p` has NO v1 lane - every q family reads `:lower-tail` only, and `:log-p` was silently ignored (since v1), so `qgamma(log 0.5, 2, 1, {:log-p true})` corrupted to an out-of-domain raw p and NaN-poisoned with no error. Now VALIDATION, loud: `stats: <fn>: :log-p is not supported in v1 - pass the raw p in [0,1]` (the sibling-option rejection class). A raw p is the only supported form; the d/p families' log lane is `:log`, a different key |
| tail precision (no key - a documented domain claim) | the q families | - | The q lane is R-exact over the legal domain EXCEPT the recorded boundaries: (a) `qf(1-1e-12,1,1)`-class cells - RESOLVED (the qbeta direct lane's log-x Newton, R qbeta.c's use_log_x port): the old constant 1.126e15 plateau across p 1-1e-11..1-1e-16 is gone, every band cell tracks R at rel ≤ 2e-13 (R 4.05e21..3.29e31; the retired deviation, superseded by the qbeta-complement suite rows); the remaining extreme-(a,b) floor: seeds that UNDERFLOW at x0 = (tgt·a·B(a,b))^{1/a} (the a < ~1e-4 medians, true root below the double grid - e.g. qbeta(0.5,1e-8,1): I_x(a,1) = x^a puts the median at e^-6.9e7) disengage the log-x lane and return the bracket floor 8.88e-16, a DOCUMENTED floor (R floors at 2^-1024 and warns "full precision may not have been achieved" - both unreferenced); the boundary band around a=b=0.05 is measured at the p=1e-16 lane - a=b=0.1 stays R-exact (rel 2.9e-14), a=b=0.05 underflows the subnormal root (ours 9.71e-315 where R keeps 5.56e-309, accept-knowingly), a=b≤0.04 returns the floor 8.88e-16 where R underflows to 0 (both unreferenced no-answers) - rows in the qbeta-floor-gap suite rows; the extreme-(a,b) interiors where R itself warns are ACCEPT-KNOWINGLY: at (0.5,1e-8,1e-8) ours is the arcsine-correct symmetric median 0.5 where R prints 1.87e-41 with its own warning (reference divergence); (b) the HUGE-SHAPE family (rewritten from an earlier wrong mechanism text - measured 2026-10-13): the mechanism is gamma-p's SERIES TRUNCATION at x≈shape (the series needs >1000 iterations from shape ~2e4 at x=shape; R routes the x~alph region to ppois_asymp/bd0 asymptotics, pgamma.c decoded - our series-only port is the boundary; the deep bd0/stirlerr port is recorded future work). Domain, measured: pgamma itself degrades from shape ~1e5 (pgamma(1e5,1e5) rel 1.58e-3, (1e6,1e6) rel 0.32 - boundary rows in the huge-shape suite rows); the q lane rides the seed + the Newton guard-rescue, so qchisq stays FINITE with per-cell measured imprecision over df 2e5..4e6 (qchisq(0.5,2e5) rel 6.3e-6, (0.5,3e5) 4e-13 - a jagged, seed-dependent onset, (0.5,1e6) 1.3e-10, (0.5,4e6) 7.9e-6, (0.025,1e7) 8.4e-5); the qgamma q-lane band is pinned cell-by-cell (truncation onset ~shape 2e4, rel climbing 4.6e-8 at 5e4 / 6.3e-6 at 1e5 / ~4e-4 at 2e6..3.5e6 - rows in the huge-shape-loud suite rows); at p=0.5 the solver DIVERGES from shape ≳ 4e6 (df 8e6) and qgamma/qchisq/qf ERROR LOUD with the scoped huge-shape message (the old NaN→0 clamp laundered a confident 0 - qchisq(0.5,8e6)=0, qf(0.025,1e9,1e9)=0, vs finite R); the loud set is NOT all interior p - it is a NARROW p-band centred near p=0.5 that WIDENS with shape (p-ladders live-measured: shape 1e7: p∈[0.46,0.51] throws; 2e7: [0.44,0.50] throws); below-band cells stay silent with in-band imprecision (1e7 p=0.44 rel 2.5e-4) and interior p at/above the band's upper edge stays FINITE (rel ~1e-9 at the edge decaying to 1e-11..1e-13 by p=0.55) - rows in the huge-shape-ladder suite rows; (c) qt's AS 396 polish skipped at df ≥ 1e7 - our dt/beta-inc precision boundary; the seed alone tracks R to <1e-10 at the 1e10..1e20 bar cells; the [1e6,1e7) band BELOW the guard carries a ~2e-9 residual (measured 2026-10-13: qt(0.975,9.5e6) rel 2.18e-9, (0.975,9.99e6) 1.89e-9 - same order as the earlier 3.34e-9/2.89e-10, NOT a regression; at/above 1e7 the guard is rel 0 EXACT, better than before); (d) a KNIFE EDGE (accepted knowingly): one qpois fuzz boundary cell (p ≈ 1-7e-16, lambda 0.1) may differ from R by one count - R itself is ulp-arbitrary there (its own q_DISCRETE fuzz rounds the other way); no row claims coverage |
| `:conf-level` | t-tests, lm, glm, bayes-lm, conj posteriors, binom-test, poisson-test, sign-test, prop-test, epi-2x2, cohens-d, hedges-g, cor-test (pearson) | 0.95 | interval coverage |
| `:alternative` | t-tests, perm, ks-test, cor tests, binom-test, poisson-test, sign-test, prop-test, cohens-d, hedges-g, cor-test | `:two-sided` | `:two-sided/:greater/:less` |
| `:correct` | prop-test, mcnemar-test, mantelhaen-test | true | the continuity correction - THREE gates: prop-test's YATES cap (min(0.5, deviation)); mcnemar subtracts 1 only when the discordant cells differ; mantelhaen 0.5 only when |DELTA| >= 0.5 |
| `:exact` | ks-test, cor-test (kendall/spearman) | nil (auto) | force the p lane: ks two-samp nx·ny < 10000; kendall n < 50 AND no ties; spearman n ≤ 1290 AND no ties (ties → the approx lane silently; R warns - the advisory convention) |
| `:continuity` | cor-test | false | kendall: S = sign(S)(|S|−1); spearman: the t-approx den + 1 |
| `:method` | p-adjust (`:bonferroni :holm :hochberg :bh :by :none`, +`:fdr`), epi-2x2 (`:cohort-count` default / `:case-control`), cor-test (`:pearson` default / `:spearman` / `:kendall`) | `:holm` / `:cohort-count` / `:pearson` | TRIPLE semantics: adjustment family vs study parameterisation vs correlation method - the three fns never share a call |
| `:n` | p-adjust | `length(ps)` | the comparison count (n >= lp) |
| `:mu` | sign-test (the null center), cohens-d, hedges-g (the null SHIFT: est = (mean-diff − mu)/s - R's mu) | 0 | per-fn meaning noted |
| `:paired` | cohens-d, hedges-g | false | d = mean(x−y)/sd(x−y); resident: the diffs are vmap2 :sub of staged inputs |
| `:pooled-sd` | cohens-d, hedges-g | true | pooled sp vs the sqrt((s1²+s2²)/2) variant (df n1+n2−2 vs Welch-Satterthwaite) |
| `:ncp` | pt | nil (central) | the noncentral-t CDF route (R's pt(q, df, ncp)); AS 243 twin series - composition, comp-only |
| `:nranges` | ptukey, qtukey | 1 | R's rr - the number of studentized ranges (the two-way extension); nmeans stays the per-range mean count |
| `:mean/:sd`, `:df`, `:df1/:df2`, `:shape1/:shape2`, `:shape/:rate` | ks-test one-sample (the named-CDF params) | `:pnorm` 0/1; others required | `y` = `:pnorm :pt :pchisq :pf :pbeta :pgamma`, or an arbitrary CDF fn |
| `:var-equal` | t-test-two | false | pooled vs Welch |
| `:ties` | coxph | `:efron` | `:efron/:breslow` |
| `:exact` | perm-test-two | auto (≤20000 splits) | force exact/MC lane |
| `:B`, `:seed` | bootstrap, perm MC, samplers | 2000/999, "stats" | replicates, deterministic seed |
| `:max-it`, `:epsilon` | glm, coxph, clogit, NM | 25/200, 1e-8..1e-10 | iteration bound, tolerance |
| `:steps`, `:burn`, `:adaptive`, `:prop-sd` | mh | 5000/1000/true/1 | chain shape (burn discipline: SPEC mc-chain.qnt) |
| `:m0/:k0/:a0/:b0/:l0` | gibbs, bayes-lm | flat/Jeffreys | priors |

## 2. Input dispatch (what a fn accepts)

- **Resident (kernel-backed)**: packed handles chain through without
  unpacking (stats/num's `stage`); janet arrays stage in (one boundary
  copy). Vector+`:num` = kernels. The full-marker gate
  validates SHAPE, not provenance, and bounds shape INCLUDING buffer
  sufficiency: a marker with any handle keys present must carry a :ptr
  buffer of at least :bytes bytes (length :ptr >= :bytes - the
  kernel validates :bytes == :n x 8 (element width) BEFORE reading,
  then reads :n x element-width; a short :ptr with a width-consistent
  :bytes is an OOB read (measured: silent junk small,
  segfault large) - the gate's length :ptr >= :bytes closes it, and
  with the kernel's width check the pair is jointly sound: no
  gate-passing marker can cause an out-of-bounds read); struct OR
  table markers (the dist guards cover tables too
  - a table carrying handle keys gets the same
  scoped error, not the unscoped kernel error; a table with ALL
  handle keys and a sufficient :ptr is ALSO rejected scoped, because
  real pack-f64 handles are structs and a table is never valid
  input) that miss any handle key, lack :ptr, or carry a short :ptr
  reject with the scoped caller-misuse error. A sufficient struct
  marker with caller-supplied :ptr bytes proceeds: the num/stage
  path returns the bytes' IEEE-double interpretation - provenance and
  content of a sufficient caller-supplied buffer are the caller's
  responsibility. The gate bounds shape (now including sufficiency),
  never content: malformed marker shapes and insufficient buffers
  error with the scoped caller-misuse error (measured:
  the pre-fix short-:ptr segfault, exit 139, and the silent small
  overreads are closed); width-inconsistent :bytes on an otherwise
  sufficient marker surfaces the kernel's element-width error
  (unscoped, pre-existing); caller-supplied buffer content is
  unpacked as data.
- **Scalar inputs**: run the IDENTICAL pure formulas under
  `:num`/`:impl` (honest scalar policy - kernels are a vector play).
- **Composition (documented)**: R-exact recurrences (pbinom),
  per-element sign branches (pt/pf), monotone inversions (q-family),
  sorts (quantiles/ranks), scatter-adds (kruskal), user-fn loops
  (bootstrap/jackknife), matched sets (clogit) - `:num` runs the
  identical janet code. Packed handles are REJECTED with the
  unpack-first error (composition needs janet-side values).
- Scalar-or-array APIs (dist d/p/q): scalars under `:num` = the pure
  formula; NaN-in-NaN-out on vector+`:num` (no validation pass),
  with the DOCUMENTED exception: the dbinom/dpois/dnbinom support
  convention: near-integers within R's
  relative tolerance |x - m| <= 1e-9*max(1,|x|), m = floor(x+0.5)
  (R's nmath R_NON_INT_TOL - ROUND semantics, boundary-measured
  from BOTH sides: x=2: 2±2e-9, x=7: 7±6.9e-9, x=200: 200±1.9e-7)
  decode to the integer; beyond the tolerance AND out-of-support
  k -> 0 (log: -inf; R's own value + warning - an early NaN here was
  a mis-decode, reversed 2026-10-09 on live evidence). Vector
  lanes: exact-integer arrays run the kernel; any other element
  falls to the elementwise pure path (kernel-ok gate).
- **OBSERVABILITY**: R emits
  separation/support warnings (glm "fitted probabilities numerically
  0 or 1", dbinom non-integer); janet has no warning channel -
  VALUES are the contract, silence is deliberate. The ONE exception:
  a stderr WARNING (once per process, at the first dynamic `:impl`
  call) fires when the janet-num module was not findable AND libjn
  is detectably on disk ($JANET_NUM_LIB path or the pinned
  ~/.local/lib/janet-num/libjn.so) - the half-installed state where
  the documented remedies apply. Genuinely-absent-libjn pure-only
  machines stay silent (the documented deterministic default).
  Explicit `{:impl :pure}` never warns.
- **CDF POLES / q-ENDPOINTS (live
  R-verified) + THE TINY-NEGATIVE ZONE +
  Q-LANE POLES:**
  R-verified 2026-10-08**: scalars follow R at the poles -
  pchisq/pgamma/ppois/pbinom/pbeta at +inf → 1; negative or
  out-of-support scalars → 0 (no throw - the old range throws are
  withdrawn; pbinom k>n saturates at 1, R's convention); the raw
  v<0 pole runs BEFORE any coercion (the -1e-10..-1.5e-7 zone is
  0, R - not P(X<=0)); ALL q families return NaN elementwise for
  p outside [0,1] (R's NaN+warning convention, both sides - the
  whole-vector throws are withdrawn); qnorm/qt
  0/1 → ∓inf, qchisq/qgamma 0/1 → 0/+inf, qbeta 0/1 → 0/1 - no
  bracket-saturation values on legal p. Vector+`:num` keeps the
  NaN-trusts contract at ±inf elements (kernel domain); the pure
  vector lane composes the scalar poles via vapply.

## 3. Error taxonomy (four classes - asserted by the fuzz/smokes)

| class | source | examples | message shape |
|---|---|---|---|
| VALIDATION | caller misuse | empty arrays, `n <= p`, p out of [0,1], non-integer k, ragged matrices | `stats: <fn> needs ...` - states the requirement |
| ALGORITHM | the estimator's own failure | glm/coxph/clogit non-convergence, non-finite deviance (boundary), not-SPD chol, trim removes everything | names the cause + the knob (`:max-it`, boundary) |
| ENVIRONMENT | the machine | strict `{:impl :num}` without libjn | names the installer (`libs/install-janet-num.sh` / `$JANET_NUM_LIB`) - SPEC impl-resolution.qnt |
| COMPOSITION | wrong input type for the lane | packed handle into a composition/sort-bound fn | "unpack packed handles first" |

Every message is actionable: it names the offending value or the fix.
No silent fallthrough between impls (one-path discipline).

## 4. Outcome families (iterative estimators)

- **Error-terminating** (SPEC: iterative-solver.qnt): glm, coxph,
  clogit - results ONLY when converged; failure has exactly two
  causes (non-finite objective / iteration bound). Nelder-Mead is
  NOT in this family: it returns `{:converged bool}` and never errors
  on non-convergence (the flag family).
- **Chain samplers** (SPEC: mc-chain.qnt): mh/gibbs return draws of
  length `steps - burn` (min 0); adaptation only during burn-in;
  determinism per `:seed` (RNG paths are never cross-compared -
  moments/convergence only).

## 5. Return shapes

Result structs carry the R-convention names (`:coefficients`, `:se`,
`:p`, `:conf-interval`, ...) plus `:iters`/`:converged` where the
algorithm iterates. Simulation entries carry their draws. Vector
outputs are janet arrays; handles never leak to the user (unpack at
the pipeline's end).

Additions (the CRAN-oracle pair): `epi-2x2`
(table opts; rows = EXPOSURE, cols = OUTCOME - epiR's as.columns
layout) → cohort-count: `{:rr :or :rd :rd-score}` each
`{:estimate :conf-interval}` - RR log/Katz, OR Woolf, RD Wald +
the Miettinen-Nurminen score interval (the halving bisection
mirrored verbatim); zero-cell estimates NaN the RR/OR triples
(epiR's guard; RD stays finite); estimates are PROPORTIONS (epiR's
units=100 is display scale); `:case-control` → `:or` ONLY (the
cohort measures are undefined - keys absent). `cohens-d` /
`hedges-g` (xs ys opts) → `{:estimate :conf-interval :alternative
:conf-level :mu :paired :pooled-sd :df :se}` (+ `:J` on hedges-g):
the exact noncentral-t CI (the .get_ncp_t equation via exact
monotone inversion; effectsize's optim carries ~1e-8 solver noise -
the fixture CI rows tolerance ABS); g = J·d with the exact gamma J,
the CI bounds scaled by the same J; one-sided = the 2·conf−1 level
solves with the unbounded side ±inf; CONSTANT input (zero sd) →
VALIDATION (R silently yields ±Inf/NA - our convention).
`pt {:ncp d}` - R's pnt.c verbatim (AS 243): the reflection, the
del>37.63/df>4e5 normal escape, the twin series (errbd 1e-12
gate); verified vs R's pt to 4e-13 (the reflected tails - C's
LDOUBLE headroom; direct-branch rows bitwise).

Additions: `cor-test` (xs ys opts) → pearson:
`{:estimate :statistic :parameter (df) :p-value :conf-interval
:null-value}` - the Fisher-z CI (tanh'd bounds; one-sided
[−1, hi] / [lo, 1]: R tanh's the WHOLE cint, tanh(±inf) = ±1),
n > 3 only, the |r| = 1 pole gives t ±inf / p 0 / CI [±1, ±1];
kendall/spearman: `{:estimate :statistic :p-value :null-value}`
- the exact DPs (kendall: the ckendall DP; spearman: brute
histogram n ≤ 9 / the AS-89 EDGEWORTH series n ≥ 10 - R's
"exact" IS Edgeworth there, the decode finding) with the auto
S-approx / t-approx on ties or `:exact false` (kendall n ≥ 50
defaults to the approx too); CONSTANT input mirrors R's NA as
NaN triples (incl. the kendall n=2-tied v2 0/0 → NaN); pearson
is kernel-RESIDENT (packed handles chain); spearman/kendall are
sort-bound (the unpack-first COMPOSITION guard). The DPs live in
ranksum (pure, no opts): `kendall-p` (R's pkendall: P(T ≤ q)),
`spearman-p` (R's prho: P(S < is) lower / P(S ≥ is) upper).

Additions: `mcnemar-test` (table opts) →
`{:statistic :parameter (df=1) :p-value}` - the paired-proportions
2×2: cc subtracts 1 ONLY when n12 ≠ n21 (R's any(x − t(x) ≠ 0)
gate); zero discordant pairs → NaN (R's 0/0; our explicit guard).
`mantelhaen-test` (tables opts; K 2×2 strata, rows-major
@[a b]/@[c d]) → `{:statistic :parameter :p-value :estimate
(conf-interval) :alternative :null-value 1}` - the MH chi-square
(the 0.5 cc caps at |DELTA| ≥ 0.5, R's pattern), common OR +
the Robins-Breslow-Greenland CI, one-sided CIs [0, hi]/[lo, inf]
- PLUS the epiR-provenance charter keys: `:rr`
(Greenland-Robins pooled RR + log CI) and `:homogeneity`
({:rr {:statistic :p-value} :or {...}} - the inverse-variance
weighted log-scale Q tests, df K−1; the OR lane via the
HALDANE-ANSCOMBE +0.5 correction; nil at K = 1). Both fns are
scalar compositions; both-impl rows assert identical results.

Additions: `ptukey` (q nmeans df opts) and
`qtukey` (p nmeans df opts) - the studentized range, R's
ptukey.c/qtukey.c VERBATIM (Copenhaver-Holland quadrature; the
AS 70 qinv + the eps 1e-4 secant - qtukey rows tolerance the
secant's stopping noise, ~1e-7 observed; the roundtrip is bounded
by the same eps). `:nranges` (R's rr, default 1) exposes the
two-way lane. VALIDATION ours (R NaNs with a warning): df >= 2,
nmeans >= 2, nranges >= 1, qtukey p in [0,1] (p = 0 → 0, p = 1 →
+inf mirror R).

Measured live vs R 4.5.0 (2026-10-13): the near-1 zone is R's OWN
no-answer / solver-floor zone - documented divergences, not chased.
(1) CONVENTION SPLIT: at p = 1-1e-16 (4,12, lower lane) R returns
NaN with "convergence failed in 'qtukey'"; ours returns the endpoint
0 - both are no-answers (R's secant gives up inside its own eps
1e-4 tolerance, ours returns the p=1 pole); documented divergence,
no row claims R-agreement there. (2) THE FLIPPED LANE'S SOLVER-FLOOR
BAND: the upper lane at near-1 p sits at the secant's floor, where
BOTH answers are noise - qtukey(1-1e-16,4,12,up): ours
1.7613248589497419e-15 vs R 4.0425319937476875e-15 (rel 0.56;
(1-1e-12,up): ours 1.586425300445987e-11 vs R 3.6411085667468211e-11;
(1-1e-8,up): 1.5864603032054689e-07 vs 3.64118891721706e-07 - the
same ~0.44x offset class, the AS 70 qinv + eps 1e-4 secant floor);
the interior flipped lane carries the ordinary secant band
(qtukey(1-1e-4,4,12,up) rel 3e-6 vs the tails' ~5e-9); the suite's
flipped-lane rows pin OUR values with the class reason, never an
R-agreement claim.

**Step 4 (janet-num 1.1.0)**: `ptukey` is now
kernel-backed on vectors - the pnorm-v dispatch class (RESIDENT:
packed handles chain, arrays stage once, outputs unpacked; the
`ptukey-v` kernel mirrors the R-verified pure, parity bar ≤1e-9
rel, measured ~1e-13/1.5e-11-worst). Scalar q under `:num` = the
identical pure formula (the honest scalar policy). The NaN pole
DIVERGED pre-1.2.0 (pure NaN - R's ISNAN gate; kernel 1, the
mirrored pre-fix decode); the kernel's 4-arg ISNAN
gate landed (janet-num 1.2.0+, returning the caller's NaN) - both lanes
now NaN, R's own value, agreement suite-pinned. `qtukey` stays
scalar composition (the secant is caller-side loop over ptukey -
no vector surface in R either).

aov additions: `aov` (groups opts) → the lm fit
PLUS `:k :means :group-n :df-residual :mse :table` - the table is
`{:group {:df :ss :ms :f :p} :residual {:df :ss :ms}}` (the
orthogonal decomposition; R's anova(aov) row values, fixture
rows pin them); rides the lm machinery (the fit routes
:impl; the table F == lm's overall F, suite-pinned). `oneway-test`
(groups opts) → `{:statistic :parameter [num-df den-df] :p-value
:var-equal}` - `:var-equal` DEFAULTS FALSE (R's default is the
WELCH lane, like t-test-two); classic lane = the same F as the aov
table; VALIDATION before the branch, both lanes: k >= 2, n_i >= 2
(R's messages). `tukey-hsd` (aov-fit opts) → `{:pairs
(:i :j :label :diff :lwr :upr :p) :conf-level}` - pair order =
R's lower-triangle COLUMN-MAJOR (k=3: 2-1, 3-1, 3-2; later level
minus earlier, level ORDER unless `:ordered`); p adj = ptukey
UPPER tail = 1 − ptukey(|est|, k, df.res) (R's lower.tail=FALSE
flip - the decode detail); ONE vector ptukey call for the whole
table (`:num` routes ptukey-v; parity bar 1e-9); qtukey stays
scalar (all pairs share the critical value) - the num lane's
remaining cost (bench.md). POLES (suite-pinned): mse = 0 with
unequal means → est ±inf → p 0 (R); all-identical y → est NaN →
p NaN in BOTH lanes (the ptukey NaN pole - kernel
agreement post-1.2.0); janet division THROWS where R
yields Inf/NaN - the explicit guards (aov F, oneway F, the welch
w_i = n_i/v_i zero-variance group → R's Inf → NaN F/p); df.res <
2 → tukey-hsd throws ptukey's VALIDATION; singleton groups LEGAL
in aov (unlike oneway.test). The pf +inf pole fixed to R's 1
(was NaN - found by the mse=0 probe).

ecdf-ks additions: `ks-test` → `{:statistic :p-value
:alternative :exact :ties}` (R-convention keys; `:exact` = the R
exact flag; `:ties` is the documented extension carrying R's tie
warning). `ecdf` returns a CALLABLE (R's ecdf): scalar or indexed
input → scalar / array of proportions, right-continuous, ties
share one step value.

survival additions: `survfit-km` → `{:n :time :surv :n-risk
:n-event :n-censor}` - aligned arrays over ALL distinct times (R
survfit's output contract; censor-only times carry S forward; the
risk sets INCLUDE same-time censorings - the km-logrank.qnt
convention). `survdiff-logrank` → `{:statistic :p-value :observed
:expected :variance :n :df}` (two-group, rho = 0; chi via
(O1−E1)²/V, p via pchisq). No new cross-fn opts keys.

exact-test additions: `p-adjust` → array of adjusted p (length preserved,
R's stable-tie ordering); `binom-test` (x n p opts) → `{:statistic
:parameter :p-value :conf-interval :estimate :null-value
:alternative}` (CP CI; the relErr two-sided convention) - the
sign-test's base, prop-test completes it; `poisson-test` (x t r opts;
x/t scalar or 2-array) → the same shape (estimate = rate or rate
ratio; the SMR lane is the one-sample call with t = expected);
`sign-test` (xs opts) → the binom-test shape + `:n-zero`.

proportion-test additions: `prop-test` (x n p opts; x/n scalar or 1-2 array, p
scalar or 2-array) → the same htest shape with `:parameter` = df
(k−1 pooled, k given - k=1 is ALWAYS 1); `:estimate` = p-hat or the
pair; `:null-value` = nil on the k=2 equality lane; `:conf-interval`
= @[lo hi] - NIL for k=2 with p given (R gives no CI there);
`:alternative` = the EFFECTIVE alternative (k=2 given p is forced
two-sided, R verbatim). The YATES cap: the correction is capped by
the observed deviation from the null (k=1 |x−np|; k=2 pooled
|DELTA|/(1/n1+1/n2)) - the cc CI is null-DEPENDENT; k=2 given p is
UNCAPPED. Non-integer x allowed (no rounding - R verbatim, unlike
binom-test); p in (0,1) EXCLUSIVE; two-sided p via the far-tail
accurate gamma-q upper chi-square, one-sided via the exact erfc
normal route.
