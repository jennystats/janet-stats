# verification probe (janet side). Usage:
#   JANET_PATH=<modules-dir> janet probe.janet > out.txt
# Every cell %.17g, keyed, one per line: name<TAB>value
(import jennystats/typed :prefix "")
(import jennystats/kit :prefix "")
(import jennystats/stats/dist :prefix "")
(import jennystats/stats/utils :prefix "")

(def UP {:lower-tail false})

(defmacro cell [name expr]
  ~(try (printf "%s\t%.17g" ,name ,expr)
     ([e] (printf "%s\tTHROW: %s" ,name (string e)))))

# rounded-complement lane
(cell "qgamma_1e-16_5_1_up" (qgamma 1e-16 5 1 UP))
(cell "qgamma_1e-12_5_1_up" (qgamma 1e-12 5 1 UP))
(cell "qgamma_5e-17_5_1_up" (qgamma 5e-17 5 1 UP))
(cell "qchisq_1e-16_5_up" (qchisq 1e-16 5 UP))
(cell "qchisq_1e-16_10_up" (qchisq 1e-16 10 UP))
(cell "qf_1e-16_400001_5_up" (qf 1e-16 400001 5 UP))
# control: the complement pair (ours upper at p vs ours lower at 1-p)
(cell "qgamma_1-1e-12_5_1" (qgamma (- 1 1e-12) 5 1))
(cell "qchisq_1-1e-16_10" (qchisq (- 1 1e-16) 10))

# huge-df / gamma-p truncation
(cell "gamma-p_1e5_1e5" (gamma-p 1e5 1e5))
(cell "gamma-p_1e6_1e6" (gamma-p 1e6 1e6))
(cell "gamma-p_1e7_1e7" (gamma-p 1e7 1e7))
(cell "qchisq_0.5_1e5" (qchisq 0.5 1e5))
(cell "qchisq_0.5_4e5" (qchisq 0.5 4e5))
(cell "qchisq_0.5_1e6" (qchisq 0.5 1e6))
(cell "qchisq_0.5_8e6" (qchisq 0.5 8e6))
(cell "qchisq_0.5_1e7" (qchisq 0.5 1e7))
(cell "qchisq_0.5_1e10" (qchisq 0.5 1e10))
(cell "qchisq_0.025_1e7" (qchisq 0.025 1e7))
(cell "qf_0.025_1e9_1e9" (qf 0.025 1e9 1e9))
(cell "qf_0.95_400001_4e5" (qf 0.95 400001 4e5))
(cell "qf_0.5_4e5_400001" (qf 0.5 4e5 400001))
(cell "qf_0.025_4e5_4e5" (qf 0.025 4e5 4e5))

# p=0.5 large-shape
(cell "qgamma_0.5_2e5_0.5" (qgamma 0.5 2e5 0.5))
(cell "qgamma_0.5_5e5_0.5" (qgamma 0.5 5e5 0.5))
(cell "qgamma_0.5_1e6_0.5" (qgamma 0.5 1e6 0.5))
(cell "qchisq_0.5_2e5" (qchisq 0.5 2e5))
# control: upper branch (p=0.55 far from the boundary)
(cell "qgamma_0.55_2e5_0.5" (qgamma 0.55 2e5 0.5))

# qbeta lower-lane floor
(cell "qbeta_1e-8_0.5_0.5" (qbeta 1e-8 0.5 0.5))
(cell "qbeta_1e-12_0.5_0.5" (qbeta 1e-12 0.5 0.5))
(cell "qbeta_1e-16_0.5_0.5" (qbeta 1e-16 0.5 0.5))
(cell "qbeta_0.5_1e-8_1" (qbeta 0.5 1e-8 1))
(cell "qbeta_0.5_1e-8_1e-8" (qbeta 0.5 1e-8 1e-8))
(cell "qf_5e-324_4_7_up" (qf 5e-324 4 7 UP))

# discrete lane
(cell "qpois_1e-16_100" (qpois 1e-16 100))
(cell "qpois_1e-16_1000" (qpois 1e-16 1000))
(cell "qpois_1e-16_10000" (qpois 1e-16 10000))
(cell "qpois_2e-15_100" (qpois 2e-15 100))
(cell "qpois_1-1e-16_100_up" (qpois (- 1 1e-16) 100 UP))
(cell "qbinom_1e-16_1000_0.9_up" (qbinom 1e-16 1000 0.9 UP))
(cell "qbinom_1e-16_100_0.5_up" (qbinom 1e-16 100 0.5 UP))
(cell "qbinom_1e-16_50_0.4_up" (qbinom 1e-16 50 0.4 UP))
(cell "qnbinom_1e-16_100_0.3_up" (qnbinom 1e-16 100 0.3 UP))
(cell "qpois_1e-300_5_up" (qpois 1e-300 5 UP))
(cell "qnbinom_1e-300_10_0.3_up" (qnbinom 1e-300 10 0.3 UP))
(cell "pbinom_960_1000_0.9" (pbinom 960 1000 0.9))
(cell "qbinom_1-5e-15_48_0.5_up" (qbinom (- 1 5e-15) 48 0.5 UP))

# qlnorm endpoints flipped
(cell "qlnorm_0_0_1_up" (qlnorm 0 0 1 UP))
(cell "qlnorm_1_0_1_up" (qlnorm 1 0 1 UP))

# qtukey near-1
(cell "qtukey_1-1e-15_4_12" (qtukey (- 1 1e-15) 4 12))
(cell "qtukey_1-1e-16_4_12" (qtukey (- 1 1e-16) 4 12))
(cell "qtukey_1-1e-16_4_12_up" (qtukey (- 1 1e-16) 4 12 UP))
(cell "qtukey_1e-16_5_10_up" (qtukey 1e-16 5 10 UP))

# qt polish band
(cell "qt_0.975_9.5e6" (qt 0.975 9.5e6))
(cell "qt_0.975_9.99e6" (qt 0.975 9.99e6))
(cell "qt_0.975_1e7" (qt 0.975 1e7))
(cell "qt_0.975_1e8" (qt 0.975 1e8))

# reference: suite constant cell already covered (qchisq_1-1e-16_10)

# handle gate (error capture, not numeric)
(each [nm form] [["dpois_partial" ~(dpois {:n 2} 1)]
                 ["dexp_partial" ~(dexp {:a 1} 1)]]
  (def r (try (string (eval form)) ([e] (string "ERR: " e))))
  (printf "%s\t%s" nm r))
