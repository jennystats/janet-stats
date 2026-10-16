# verification probe (R side). Rscript rprobe.R
f <- function(name, v) cat(sprintf("%s\t%.17g\n", name, v))
UP <- list(lower.tail=FALSE)
f("qgamma_1e-16_5_1_up", qgamma(1e-16, 5, 1, lower.tail=FALSE))
f("qgamma_1e-12_5_1_up", qgamma(1e-12, 5, 1, lower.tail=FALSE))
f("qgamma_5e-17_5_1_up", qgamma(5e-17, 5, 1, lower.tail=FALSE))
f("qchisq_1e-16_5_up", qchisq(1e-16, 5, lower.tail=FALSE))
f("qchisq_1e-16_10_up", qchisq(1e-16, 10, lower.tail=FALSE))
f("qf_1e-16_400001_5_up", qf(1e-16, 400001, 5, lower.tail=FALSE))
f("qgamma_1-1e-12_5_1", qgamma(1-1e-12, 5, 1))
f("qchisq_1-1e-16_10", qchisq(1-1e-16, 10))
f("gamma-p_1e5_1e5", pgamma(1e5, 1e5))
f("gamma-p_1e6_1e6", pgamma(1e6, 1e6))
f("gamma-p_1e7_1e7", pgamma(1e7, 1e7))
f("qchisq_0.5_1e5", qchisq(0.5, 1e5))
f("qchisq_0.5_4e5", qchisq(0.5, 4e5))
f("qchisq_0.5_1e6", qchisq(0.5, 1e6))
f("qchisq_0.5_8e6", qchisq(0.5, 8e6))
f("qchisq_0.5_1e7", qchisq(0.5, 1e7))
f("qchisq_0.5_1e10", qchisq(0.5, 1e10))
f("qchisq_0.025_1e7", qchisq(0.025, 1e7))
f("qf_0.025_1e9_1e9", qf(0.025, 1e9, 1e9))
f("qf_0.95_400001_4e5", qf(0.95, 400001, 4e5))
f("qf_0.5_4e5_400001", qf(0.5, 4e5, 400001))
f("qf_0.025_4e5_4e5", qf(0.025, 4e5, 4e5))
f("qgamma_0.5_2e5_0.5", qgamma(0.5, 2e5, 0.5))
f("qgamma_0.5_5e5_0.5", qgamma(0.5, 5e5, 0.5))
f("qgamma_0.5_1e6_0.5", qgamma(0.5, 1e6, 0.5))
f("qchisq_0.5_2e5", qchisq(0.5, 2e5))
f("qgamma_0.55_2e5_0.5", qgamma(0.55, 2e5, 0.5))
f("qbeta_1e-8_0.5_0.5", qbeta(1e-8, 0.5, 0.5))
f("qbeta_1e-12_0.5_0.5", qbeta(1e-12, 0.5, 0.5))
f("qbeta_1e-16_0.5_0.5", qbeta(1e-16, 0.5, 0.5))
f("qbeta_0.5_1e-8_1", qbeta(0.5, 1e-8, 1))
f("qbeta_0.5_1e-8_1e-8", qbeta(0.5, 1e-8, 1e-8))
f("qf_5e-324_4_7_up", qf(5e-324, 4, 7, lower.tail=FALSE))
f("qpois_1e-16_100", qpois(1e-16, 100))
f("qpois_1e-16_1000", qpois(1e-16, 1000))
f("qpois_1e-16_10000", qpois(1e-16, 10000))
f("qpois_2e-15_100", qpois(2e-15, 100))
f("qpois_1-1e-16_100_up", qpois(1-1e-16, 100, lower.tail=FALSE))
f("qbinom_1e-16_1000_0.9_up", qbinom(1e-16, 1000, 0.9, lower.tail=FALSE))
f("qbinom_1e-16_100_0.5_up", qbinom(1e-16, 100, 0.5, lower.tail=FALSE))
f("qbinom_1e-16_50_0.4_up", qbinom(1e-16, 50, 0.4, lower.tail=FALSE))
f("qnbinom_1e-16_100_0.3_up", qnbinom(1e-16, 100, 0.3, lower.tail=FALSE))
f("qpois_1e-300_5_up", qpois(1e-300, 5, lower.tail=FALSE))
f("qnbinom_1e-300_10_0.3_up", qnbinom(1e-300, 10, 0.3, lower.tail=FALSE))
f("pbinom_960_1000_0.9", pbinom(960, 1000, 0.9))
f("qbinom_1-5e-15_48_0.5_up", qbinom(1-5e-15, 48, 0.5, lower.tail=FALSE))
# discrete lane class extension
f("c5_qpois_1e-16_100", qpois(1e-16, 100))
f("c5_qpois_1e-16_1000", qpois(1e-16, 1000))
f("c5_qpois_1e-16_10000", qpois(1e-16, 10000))
f("c5_qpois_2e-15_100", qpois(2e-15, 100))
f("c5_qpois_1e-16_5_up", qpois(1e-16, 5, lower.tail=FALSE))
f("c5_qnbinom_1e-16_10_0.3_up", qnbinom(1e-16, 10, 0.3, lower.tail=FALSE))
f("c5_qbinom_1e-16_1000_0.9_up", qbinom(1e-16, 1000, 0.9, lower.tail=FALSE))
f("c5_pbinom_960_1000_0.9", pbinom(960, 1000, 0.9))
f("c5_pbinom_999_1000_0.9", pbinom(999, 1000, 0.9))
# interior controls (both lanes, all three families)
f("c5_qpois_0.5_10", qpois(0.5, 10))
f("c5_qpois_0.975_100", qpois(0.975, 100))
f("c5_qpois_0.025_100_up", qpois(0.025, 100, lower.tail=FALSE))
f("c5_qbinom_0.5_10_0.5", qbinom(0.5, 10, 0.5))
f("c5_qbinom_0.975_100_0.9", qbinom(0.975, 100, 0.9))
f("c5_qnbinom_0.5_10_0.3", qnbinom(0.5, 10, 0.3))
f("qlnorm_0_0_1_up", qlnorm(0, 0, 1, lower.tail=FALSE))
f("qlnorm_1_0_1_up", qlnorm(1, 0, 1, lower.tail=FALSE))
f("qtukey_1-1e-15_4_12", qtukey(1-1e-15, 4, 12))
f("qtukey_1-1e-16_4_12", qtukey(1-1e-16, 4, 12))
f("qtukey_1-1e-16_4_12_up", qtukey(1-1e-16, 4, 12, lower.tail=FALSE))
f("qtukey_1e-16_5_10_up", qtukey(1e-16, 5, 10, lower.tail=FALSE))
f("qt_0.975_9.5e6", qt(0.975, 9.5e6))
f("qt_0.975_9.99e6", qt(0.975, 9.99e6))
f("qt_0.975_1e7", qt(0.975, 1e7))
f("qt_0.975_1e8", qt(0.975, 1e8))
# rows keyed to input classes
# upper-lane tiny-p ladder per complement family, shape ladder
for (sh in c(0.5, 2, 5, 20, 50)) {
  for (p in c(1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16, 5e-17)) {
    f(sprintf("c1_qgamma_%g_%g_1_up", p, sh), qgamma(p, sh, 1, lower.tail=FALSE))
  }
}
for (df in c(1, 5, 10, 40, 100)) {
  for (p in c(1e-12, 1e-16)) {
    f(sprintf("c1_qchisq_%g_%g_up", p, df), qchisq(p, df, lower.tail=FALSE))
  }
}
# qf fudge lane both orders at far-tail p
for (p in c(1e-12, 1e-16)) {
  f(sprintf("c1_qf_%g_5_400001_up", p), qf(p, 5, 400001, lower.tail=FALSE))
  f(sprintf("c1_qf_%g_400001_5_up", p), qf(p, 400001, 5, lower.tail=FALSE))
  f(sprintf("c1_qf_%g_400001_400000_up", p), qf(p, 400001, 400000, lower.tail=FALSE))
}
# control pair: lower lane near-1 must differ from upper lane tiny-p
f("c1_qgamma_1-1e-12_5_1", qgamma(1-1e-12, 5, 1))
f("c1_qgamma_1-1e-16_5_1", qgamma(1-1e-16, 5, 1))
# subnormal/swallow boundary probes
f("c1_qgamma_1e-17_5_1_up", qgamma(1e-17, 5, 1, lower.tail=FALSE))
f("c1_qgamma_1e-30_5_1_up", qgamma(1e-30, 5, 1, lower.tail=FALSE))
f("c1_qgamma_1e-300_5_1_up", qgamma(1e-300, 5, 1, lower.tail=FALSE))
# upper lane interior control (large raw upper p)
f("c1_qgamma_0.975_5_1_up", qgamma(0.975, 5, 1, lower.tail=FALSE))
f("c1_qf_0.025_4_7_up", qf(0.025, 4, 7, lower.tail=FALSE))
# large-shape family extension
for (p in c(0.5, 0.5+2^-53, 0.5-2^-53, 0.4999999999, 0.5000000001, 0.55)) {
  for (sh in c(2e5, 5e5, 1e6)) {
    f(sprintf("c3_qgamma_%.17g_%g_0.5", p, sh), qgamma(p, sh, 0.5))
  }
}
for (df in c(2e5, 4e5, 1e6)) {
  f(sprintf("c3_qchisq_0.5_%g", df), qchisq(0.5, df))
}
# the qf 4e5+1 fudge cells, both orders
f("c3_qf_0.95_400001_4e5", qf(0.95, 400001, 4e5))
f("c3_qf_0.95_4e5_400001", qf(0.95, 4e5, 400001))
f("c3_qf_0.975_400001_400000", qf(0.975, 400001, 400000))
f("c3_qf_0.5_4e5_400001", qf(0.5, 4e5, 400001))
# small-shape median controls (the branch must not move there)
f("c3_qgamma_0.5_5_0.5", qgamma(0.5, 5, 0.5))
f("c3_qgamma_0.25_5_1", qgamma(0.25, 5, 1))
f("c3_qgamma_0.5_2e4_0.5", qgamma(0.5, 2e4, 0.5))
# huge-shape family extension
# interior p at the huge-df ladder (the imprecision boundary + the
# divergence boundary - the throw cells have no R value, they pin
# the ERROR; these pin the measured imprecision class)
for (df in c(1e5, 2e5, 3e5, 1e6, 2e6, 4e6)) {
  f(sprintf("c2_qchisq_0.5_%g", df), qchisq(0.5, df))
}
f("c2_qchisq_0.025_1e7", qchisq(0.025, 1e7))
# gamma-p truncation onset (x = shape, the p-lane boundary)
f("c2_pg_1e4_1e4", pgamma(1e4, 1e4))
f("c2_pg_1e5_1e5", pgamma(1e5, 1e5))
f("c2_pg_1e6_1e6", pgamma(1e6, 1e6))
# qf fudge-lane huge-df (the silent-0 class)
f("c2_qf_0.025_1e9_1e9", qf(0.025, 1e9, 1e9))
f("c2_qf_0.95_1e9_1e9", qf(0.95, 1e9, 1e9))
f("c2_qf_0.5_1e6_1e6", qf(0.5, 1e6, 1e6))
# qbeta direct-lane floor extension
f("c4_qbeta_1e-6_0.05_0.05", qbeta(1e-6, 0.05, 0.05))
f("c4_qbeta_1e-14_0.5_0.5", qbeta(1e-14, 0.5, 0.5))
f("c4_qbeta_1e-8_0.5_0.5", qbeta(1e-8, 0.5, 0.5))
f("c4_qbeta_1e-12_0.5_0.5", qbeta(1e-12, 0.5, 0.5))
f("c4_qbeta_1e-16_0.5_0.5", qbeta(1e-16, 0.5, 0.5))
# the documented-floor cell (R's own floor + warning; true root below
# the double grid - I_x(a,1) = x^a puts the median at e^-6.9e7)
f("c4_qbeta_0.5_1e-8_1", qbeta(0.5, 1e-8, 1))
# the qf(1-1e-12,1,1) plateau (api-contract line cells)
for (p in c(1-1e-11, 1-1e-12, 1-1e-14, 1-1e-16)) {
  f(sprintf("c4_qf_%.17g_1_1", p), qf(p, 1, 1))
}
