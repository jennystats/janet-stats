# janet-stats

Statistics library for Janet: distribution functions, hypothesis
tests, regression, survival analysis, resampling and Bayesian
estimation - 146 functions, pure Janet, R-convention semantics.

Built for statistical computing without leaving the Janet process.
The core library is self-contained: the `typed` module it uses is
vendored in-repo - no install step. Part of [Jenny
Stats](https://github.com/jennystats), a home for Janet data and
statistics packages. Jenny is just Janet.

Optional [janet-num](https://github.com/jennystats/janet-num)
integration gives packed-vector kernels (64/128/256-bit lanes) for
large inputs; without it, everything runs in pure Janet (packed-
handle inputs on kernel paths error with guidance instead of kernel
acceleration, and large-input paths are slower).

## Status

51-test suite, green. Precision boundaries are
measured, not guessed - see [Precision](#precision) below and the
contract in `specs/api-contract.md`.

## Install and running

janet-stats runs on Janet 1.42+. The core library is pure Janet and
self-contained: the `typed` module it uses is vendored in
`jennystats/typed.janet` (matching
[janet-kit](https://github.com/jennystats/janet-kit) 0.1.2+, which
carries the same docstring-hoisting behaviour - see `VENDORED.md`).
No install step:

```bash
# from the repo root
JANET_PATH=. janet -e '(import jennystats/stats :prefix "")'
```

Optional: the `janet-num` companion library (built C kernels) enables
packed-vector lanes and large-input paths. `stats` detects it
automatically. Without it, everything runs in pure Janet - packed-
handle inputs on kernel paths error with guidance instead of kernel
acceleration, and large-input paths are slower.

All jenny modules install under a `jennystats/` namespace directory:
`(import jennystats/stats)` - collision-proof against bare-name
shadowing of core symbols or other packages, by construction. Prefer
a shorter call prefix? Alias at import (e.g.
`(import jennystats/stats :as stats)`), or use the `:prefix ""` form
shown in the quick start below.

## Quick start

```janet
(import jennystats/stats :prefix "")
```

### Distributions (d/p/q families, R convention)

```janet
(dnorm 0 0 1)        # -> 0.398942280401433
(pnorm 1.96 0 1)     # -> 0.97500210485178
(qnorm 0.975 0 1)    # -> 1.95996398454005
(qgamma 0.5 2 1)     # -> 1.67834699001666
```

R semantics: `:lower-tail` defaults true. Where a function has no
implemented option, the library rejects loudly rather than returning
a silent wrong value.

### Hypothesis tests

```janet
(t-test-one [2 4 6 8 10] 5)
# -> {:estimate 6 :statistic 0.707... :p-value 0.518...
#      :stderr 1.414... :parameter 4 :conf-interval (...)}

(wilcoxon-sr @[4 6 8])          # signed-rank on differences, mu = 0
# -> {:n 3 :statistic 6 :p-value 0.25}

(binom-test 8 20 0.5)           # -> {:p-value 0.5034... ...}

(qtukey 0.95 4 10)              # studentised range at 95%: 4.32658211600917

(chisq-test observed expected)
(fisher-test-2x2 [[12 5] [7 10]])
(prop-test [45 65] [100 100])
(cor-test x y)
(oneway-test grp1 grp2 grp3)
```

`fisher-test-2x2` returns both `:odds-ratio` (the sample odds
ratio, 3.4286 on the cell above) and `:estimate` (the conditional
maximum-likelihood estimate, 3.2992); the `:estimate` is the
number R's `fisher.test` reports as `$estimate`.

Note `qtukey` takes p as a probability in [0, 1] (R convention), not
the range: `qtukey(2.0, 4, 10)` is NaN in R and NaN here.

### Regression and survival

```janet
(bayes-lm y x)          # Bayesian linear model
(glm y x)               # generalised linear model
(coxph times status covariates)  # Cox proportional hazards
(survfit-km times status) # Kaplan-Meier
(survdiff-logrank ...)  # log-rank test
```

Cox and GLM take matrix-shaped covariates (struct/table form); see
`specs/api-contract.md` for the exact call shapes.

### Resampling and Bayesian estimation

```janet
(def b (bootstrap xs (fn [d idx] (mean (map |(get d $) idx))) {:B 500}))
(boot-ci b)                     # -> percentile interval (tuple)
(perm-test-two a b)
(post-beta-binomial ...)
(post-gamma-poisson ...)
(hpd samples 0.95)
```

Bootstrap takes the data and a two-argument statistic `(fn [data
indices])` - it resamples indices and calls your statistic with the
resample; `boot-ci` then computes the percentile interval from the
`boot-out` record.

### Multivariate

```janet
(pca matrix)
(k-means k data)
(pearson x y)
(spearman a b)
(kendall a b)
```

### Descriptive

```janet
(ecdf xs)
(quantile-t7 xs 0.5)
(iqr xs)
(skewness xs)
(cohens-d a b)
```

## The function surface

146 functions across five families:

- **Distributions** (`d`/`p`/`q` prefixes, R naming): normal, t,
  chi-squared, F, binomial, poisson, beta, gamma, exponential,
  negative binomial, plus Tukey's studentised range (`qtukey`).
  Microbenchmarks: 2.83-63 microseconds per call (scalar to
  packed-lane).
- **Tests**: t (one/two/paired), Wilcoxon signed-rank, chi-squared
  (test + goodness-of-fit), Fisher exact (2x2), proportion, binomial,
  poisson, Pearson/Spearman/Kendall correlation, one-way ANOVA,
  Tukey HSD, Kruskal-Wallis, sign test.
- **Regression and survival**: Bayesian linear model, GLM, Cox
  proportional hazards, Kaplan-Meier, log-rank.
- **Resampling and Bayes**: bootstrap, bootstrap CI, jackknife,
  permutation tests, Gibbs sampling, Beta-binomial / Gamma-Poisson /
  Normal-Normal posteriors, HPD intervals.
- **Multivariate and linear algebra**: PCA, k-means, Pearson,
  Jacobi eigenvalues, QR and Cholesky decompositions, solvers.

The full contract - every function, option and measured precision
boundary - is `specs/api-contract.md`.

## Precision

Precision boundaries are measured against R 4.5.0 cell by cell, and
the contract documents every one that is not exact. Three examples
from the measured set:

- q-lane (quantile) functions are R-exact over the legal domain, with
  documented exceptions at huge shape parameters (gamma family series
  truncation from shape ~2e4; loud rejection with a scoped error
  beyond shape 4e6 at p = 0.5, where the honest answer is "the series
  does not converge" rather than a silent wrong value).
- qbeta has a documented subnormal floor near a = b = 0.05, and a
  documented bracket floor below a = b = 0.04, where R underflows to
  zero.
- Every documented deviation carries a measured relative error, the
  R reference value, and the regression rows that pin it.

The philosophy: no silent wrong values. Where the library cannot be
exact, it says so in the error or in the contract.

## Tests

The full oracle suite (R-reference regression rows, distribution
cells, survival and resampling checks) runs from the repo root:

```bash
./run_tests.sh        # RESULT: pass
```

The suite skips gracefully when the optional `janet-num` companion
is absent; the pure-Janet paths are always tested.

## License

AGPLv3 - see `LICENSE`.

## AI assistance

The code in this repository was written with the assistance of a
large language model and reviewed by its maintainer.
