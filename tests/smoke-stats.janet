# smoke-stats - conformance: stats/utils + stats/dist vs the
# RECORDED R fixtures (committed; no R needed at smoke time - pure
# janet, runs everywhere). Fixture invariant: controlled
# CSV (no quoting/NA/embedded commas) - the trivial reader is complete.
(import stats)
(import typed :prefix "")
(import kit :prefix "")

(defn assert-close [a b tol msg]
  (def diff (math/abs (- a b)))
  (def rel (if (zero? b) diff (/ diff (math/abs b))))
  (unless (or (<= rel tol) (<= diff 1e-13))
    (errorf "FAIL %s: %.17g != %.17g rel=%.2e" msg a b rel)))

(defn parse-fixture [path]
  (def out @[])
  (each line (string/split "\n" (string (slurp path)))
    (when (and (not (empty? line)) (not (string/find "fn," line)))
      (array/push out line)))
  out)

(def fixtures-dir
  (let [here (string (dyn *current-file*))]
    (string (dirname here) "/../tests/fixtures/stats")))# ---- utils fixtures: fn,a,b,x,expected -------------------------------------
(each line (parse-fixture (string fixtures-dir "/tier-b-utils.csv"))
  (def fields (string/split "," line))
  (def fn-kw (fields 0))
  (def av (scan-number (get fields 1)))
  (def bv (scan-number (get fields 2)))
  (def xv (scan-number (get fields 3)))
  (def ev (scan-number (get fields 4)))
  (def got (case fn-kw
             "lgamma" (stats/lgamma av)
             "erf" (stats/erf av)
             "erfc" (stats/erfc av)
             "gamma_p" (stats/gamma-p av xv)
             "beta_inc" (stats/beta-inc av bv xv)
             (errorf "unknown fn %s" fn-kw)))
  (assert-close got ev 1e-9 (string fn-kw " " av " " bv " " xv)))

# ---- dist fixtures: fn,args,expected ----------------------------------------
(each line (parse-fixture (string fixtures-dir "/tier-b-dist.csv"))
  (def fields (string/split "," line))
  (def fn-kw (fields 0))
  (def parts (map scan-number (string/split "," (string/join (slice fields 1 -2) ","))))
  (def ev (scan-number (last fields)))
  (def got (case fn-kw
             "dnorm" (stats/dnorm (parts 0) (parts 1) (parts 2) {:impl :pure})
             "dnorm_log" (stats/dnorm (parts 0) (parts 1) (parts 2) {:log true})
             "pnorm" (stats/pnorm (parts 0) (parts 1) (parts 2) {:impl :pure})
             "pnorm_log" (stats/pnorm (parts 0) (parts 1) (parts 2) {:log true})
             "qnorm" (stats/qnorm (parts 0) (parts 1) (parts 2) {:impl :pure})
             "dt" (stats/dt (parts 0) (parts 1) {:impl :pure})
             "dt_log" (stats/dt (parts 0) (parts 1) {:log true})
             "pt" (stats/pt (parts 0) (parts 1) {:impl :pure})
             "qt" (stats/qt (parts 0) (parts 1) {:impl :pure})
             "dchisq" (stats/dchisq (parts 0) (parts 1) {:impl :pure})
             "pchisq" (stats/pchisq (parts 0) (parts 1) {:impl :pure})
             "qchisq" (stats/qchisq (parts 0) (parts 1) {:impl :pure})
             "df" (stats/df-f (parts 0) (parts 1) (parts 2) {:impl :pure})
             "pf" (stats/pf (parts 0) (parts 1) (parts 2) {:impl :pure})
             "qf" (stats/qf (parts 0) (parts 1) (parts 2) {:impl :pure})
             "dbinom" (stats/dbinom (parts 0) (parts 1) (parts 2) {:impl :pure})
             "pbinom" (stats/pbinom (parts 0) (parts 1) (parts 2) {:impl :pure})
             "qbinom" (stats/qbinom (parts 0) (parts 1) (parts 2) {:impl :pure})
             "dpois" (stats/dpois (parts 0) (parts 1) {:impl :pure})
             "ppois" (stats/ppois (parts 0) (parts 1) {:impl :pure})
             "qpois" (stats/qpois (parts 0) (parts 1) {:impl :pure})
             "dexp" (stats/dexp (parts 0) (parts 1) {:impl :pure})
             "dexp_log" (stats/dexp (parts 0) (parts 1) {:log true})
             "pexp" (stats/pexp (parts 0) (parts 1) {:impl :pure})
             "qexp" (stats/qexp (parts 0) (parts 1) {:impl :pure})
             "dweibull" (stats/dweibull (parts 0) (parts 1) (parts 2) {:impl :pure})
             "dweibull_log" (stats/dweibull (parts 0) (parts 1) (parts 2) {:log true})
             "pweibull" (stats/pweibull (parts 0) (parts 1) (parts 2) {:impl :pure})
             "qweibull" (stats/qweibull (parts 0) (parts 1) (parts 2) {:impl :pure})
             "dlnorm" (stats/dlnorm (parts 0) (parts 1) (parts 2) {:impl :pure})
             "dlnorm_log" (stats/dlnorm (parts 0) (parts 1) (parts 2) {:log true})
             "plnorm" (stats/plnorm (parts 0) (parts 1) (parts 2) {:impl :pure})
             "qlnorm" (stats/qlnorm (parts 0) (parts 1) (parts 2) {:impl :pure})
             "dnbinom" (stats/dnbinom (parts 0) (parts 1) (parts 2) {:impl :pure})
             "dnbinom_log" (stats/dnbinom (parts 0) (parts 1) (parts 2) {:log true})
             "pnbinom" (stats/pnbinom (parts 0) (parts 1) (parts 2) {:impl :pure})
             "qnbinom" (stats/qnbinom (parts 0) (parts 1) (parts 2) {:impl :pure})
             "ptukey" (stats/ptukey (parts 0) (parts 1) (parts 2) {:impl :pure})
             "qtukey" (stats/qtukey (parts 0) (parts 1) (parts 2) {})
             "ptukey_nr" (stats/ptukey (parts 0) (parts 1) (parts 2)
                                       {:impl :pure :nranges (parts 3)})
             "qtukey_nr" (stats/qtukey (parts 0) (parts 1) (parts 2)
                                       {:nranges (parts 3)})
             (errorf "unknown fn %s" fn-kw)))
  (def tol (if (or (string/find "q" (string/slice fn-kw 0 1))
                   (string/find "_log" fn-kw))
             1e-5 1e-9))
  (assert-close got ev tol fn-kw))

# ---- API sanity: {:impl} plumbing, vectorization ---------------------------
(assert (= (stats/pnorm 0 0 1 {:impl :pure}) 0.5) "pnorm center")
# erf/erfc poles match R (oracle-finding-erf-poles.md, resolved Q1(a)
# 2026-10-04): the CF machinery returns nan at ±inf; R and the libm
# kernels give the IEEE limits
(assert (= 1.0 (stats/erf math/inf)) "erf(inf)=1")
(assert (= -1.0 (stats/erf (- math/inf))) "erf(-inf)=-1")
(assert (= 0.0 (stats/erfc math/inf)) "erfc(inf)=0")
(assert (= 2.0 (stats/erfc (- math/inf))) "erfc(-inf)=2")
(def vp (stats/pnorm [-1 0 1] 0 1 {:impl :pure}))
(assert (and (= (length vp) 3)
             (= (vp 1) 0.5)
             (<= (math/abs (- (vp 0) 0.158655253931457)) 1e-12)
             (<= (math/abs (- (vp 2) 0.841344746068543)) 1e-12))
  "vectorized pnorm")
# P1.8: every tier landed - the tranche error is RETIRED; :num answers
(assert (<= (math/abs (- (stats/dnorm 1.2 0 1 {:impl :num})
                         (stats/dnorm 1.2 0 1 {:impl :pure}))) 1e-15)
  ":num answers (tranche error retired)")
# slice 3 tranche 1: pnorm vector :num IS landed - parity with pure
(def vp-num (stats/pnorm [-1 0 1] 0 1 {:impl :num}))
(assert (and (= (length vp-num) 3)
             (<= (math/abs (- (vp-num 0) (vp 0))) 1e-12)
             (<= (math/abs (- (vp-num 2) (vp 2))) 1e-12)
             (= (vp-num 1) 0.5))
  "pnorm vector :num == pure (slice 3)")
# P1.12: the DEFAULT is :impl - the shipping path gets smoke coverage
# (on this machine it routes the kernel; on absent-libjn machines it
# resolves :pure - either way these hold)
(assert (<= (math/abs (- (stats/pnorm 1.5 0 1 {})
                         (stats/pnorm 1.5 0 1 {:impl :pure}))) 1e-15)
  "default scalar == pure formula")
(def vp-def (stats/pnorm [-1 0 1] 0 1 {}))
(assert (and (= (length vp-def) 3)
             (<= (math/abs (- (vp-def 0) (vp 0))) 1e-12)
             (<= (math/abs (- (vp-def 2) (vp 2))) 1e-12)
             (= (vp-def 1) 0.5))
  "default vector pnorm == pure values")

# ---- tier A: descriptive vs recorded R fixtures ------------------------------
(each line (parse-fixture (string fixtures-dir "/tier-a-descriptive.csv"))
  (def fields (string/split "," line))
  (def fn-kw (fields 0))
  (def argstr (string/join (slice fields 1 -2) ","))
  (def ev (scan-number (last fields)))
  (defn nums [s] (map scan-number (string/split " " s)))
  (def got
    (case fn-kw
      "mean" (stats/mean (nums argstr) {:impl :pure})
      "variance" (let [xs (nums argstr)] (stats/variance xs (stats/mean xs {:impl :pure}) {:impl :pure}))
      "sd" (let [xs (nums argstr)] (stats/sd xs (stats/mean xs {:impl :pure}) {:impl :pure}))
      "skewness" (let [xs (nums argstr)]
                   (stats/skewness xs (stats/mean xs {:impl :pure}) (stats/sd xs (stats/mean xs {:impl :pure}) {:impl :pure}) {:impl :pure}))
      "kurtosis" (let [xs (nums argstr)]
                   (stats/kurtosis xs (stats/mean xs {:impl :pure}) (stats/sd xs (stats/mean xs {:impl :pure}) {:impl :pure}) {:impl :pure}))
      "median" (stats/median (nums argstr) {:impl :pure})
      "iqr" (stats/iqr (nums argstr) {:impl :pure})
      "quantile" (let [[p vs] (string/split "|" argstr)]
                   (stats/quantile-t7 (stats/sorted-copy (nums vs)) (scan-number p) {:impl :pure}))
      "pearson" (let [[a b] (string/split "|" argstr)] (stats/pearson (nums a) (nums b) {:impl :pure}))
      "spearman" (let [[a b] (string/split "|" argstr)] (stats/spearman (nums a) (nums b) {:impl :pure}))
      "kendall" (let [[a b] (string/split "|" argstr)] (stats/kendall (nums a) (nums b) {:impl :pure}))
      "mad" (stats/mad (nums argstr) nil {:impl :pure})
      "trimmed_mean" (let [[t vs] (string/split "|" argstr)]
                       (stats/trimmed-mean (nums vs) (scan-number t) {:impl :pure}))
      # ---- cor-test ----
      "ct_pearson_est" (let [[a b] (string/split "|" argstr)] ((stats/cor-test (nums a) (nums b) {}) :estimate))
      "ct_pearson_t" (let [[a b] (string/split "|" argstr)] ((stats/cor-test (nums a) (nums b) {}) :statistic))
      "ct_pearson_p" (let [[a b] (string/split "|" argstr)] ((stats/cor-test (nums a) (nums b) {}) :p-value))
      "ct_pearson_lo" (let [[a b] (string/split "|" argstr)]
                        (get ((stats/cor-test (nums a) (nums b) {}) :conf-interval) 0))
      "ct_pearson_hi" (let [[a b] (string/split "|" argstr)]
                        (get ((stats/cor-test (nums a) (nums b) {}) :conf-interval) 1))
      "ct_pearson_less_p" (let [[a b] (string/split "|" argstr)]
                            ((stats/cor-test (nums a) (nums b) {:alternative :less}) :p-value))
      "ct_pearson_less_hi" (let [[a b] (string/split "|" argstr)]
                             (get ((stats/cor-test (nums a) (nums b) {:alternative :less})
                                  :conf-interval) 1))
      "ct_pearson_great_p" (let [[a b] (string/split "|" argstr)]
                             ((stats/cor-test (nums a) (nums b) {:alternative :greater}) :p-value))
      "ct_pearson_great_lo" (let [[a b] (string/split "|" argstr)]
                              (get ((stats/cor-test (nums a) (nums b) {:alternative :greater})
                                   :conf-interval) 0))
      "ct_pearson90_lo" (let [[a b] (string/split "|" argstr)]
                          (get ((stats/cor-test (nums a) (nums b) {:conf-level 0.9})
                               :conf-interval) 0))
      "ct_pearson90_hi" (let [[a b] (string/split "|" argstr)]
                          (get ((stats/cor-test (nums a) (nums b) {:conf-level 0.9})
                               :conf-interval) 1))
      "ct_pearson3_est" (let [[a b] (string/split "|" argstr)] ((stats/cor-test (nums a) (nums b) {}) :estimate))
      "ct_pearson3_p" (let [[a b] (string/split "|" argstr)] ((stats/cor-test (nums a) (nums b) {}) :p-value))
      "ct_kendall_est" (let [[a b] (string/split "|" argstr)]
                         ((stats/cor-test (nums a) (nums b) {:method :kendall}) :estimate))
      "ct_kendall_T" (let [[a b] (string/split "|" argstr)]
                       ((stats/cor-test (nums a) (nums b) {:method :kendall}) :statistic))
      "ct_kendall_p" (let [[a b] (string/split "|" argstr)]
                       ((stats/cor-test (nums a) (nums b) {:method :kendall}) :p-value))
      "ct_kendall_less_p" (let [[a b] (string/split "|" argstr)]
                            ((stats/cor-test (nums a) (nums b)
                                             {:method :kendall :alternative :less}) :p-value))
      "ct_kendall_great_p" (let [[a b] (string/split "|" argstr)]
                             ((stats/cor-test (nums a) (nums b)
                                              {:method :kendall :alternative :greater}) :p-value))
      "ct_kendall_ap_z" (let [[a b] (string/split "|" argstr)]
                          ((stats/cor-test (nums a) (nums b)
                                           {:method :kendall :exact false}) :statistic))
      "ct_kendall_ap_p" (let [[a b] (string/split "|" argstr)]
                          ((stats/cor-test (nums a) (nums b)
                                           {:method :kendall :exact false}) :p-value))
      "ct_kendall_ties_z" (let [[a b] (string/split "|" argstr)]
                            ((stats/cor-test (nums a) (nums b) {:method :kendall}) :statistic))
      "ct_kendall_ties_p" (let [[a b] (string/split "|" argstr)]
                            ((stats/cor-test (nums a) (nums b) {:method :kendall}) :p-value))
      "ct_kendall_cont_z" (let [[a b] (string/split "|" argstr)]
                            ((stats/cor-test (nums a) (nums b)
                                             {:method :kendall :exact false
                                              :continuity true}) :statistic))
      "ct_kendall49_p" (let [[a b] (string/split "|" argstr)]
                         ((stats/cor-test (nums a) (nums b) {:method :kendall}) :p-value))
      "ct_kendall50_z" (let [[a b] (string/split "|" argstr)]
                         ((stats/cor-test (nums a) (nums b) {:method :kendall}) :statistic))
      "ct_kendall50_p" (let [[a b] (string/split "|" argstr)]
                         ((stats/cor-test (nums a) (nums b) {:method :kendall}) :p-value))
      "ct_spearman_est" (let [[a b] (string/split "|" argstr)]
                          ((stats/cor-test (nums a) (nums b) {:method :spearman}) :estimate))
      "ct_spearman_S" (let [[a b] (string/split "|" argstr)]
                        ((stats/cor-test (nums a) (nums b) {:method :spearman}) :statistic))
      "ct_spearman_p" (let [[a b] (string/split "|" argstr)]
                       ((stats/cor-test (nums a) (nums b) {:method :spearman}) :p-value))
      "ct_spearman_less_p" (let [[a b] (string/split "|" argstr)]
                             ((stats/cor-test (nums a) (nums b)
                                              {:method :spearman :alternative :less}) :p-value))
      "ct_spearman_great_p" (let [[a b] (string/split "|" argstr)]
                              ((stats/cor-test (nums a) (nums b)
                                               {:method :spearman :alternative :greater}) :p-value))
      "ct_spearman_ap_p" (let [[a b] (string/split "|" argstr)]
                           ((stats/cor-test (nums a) (nums b)
                                            {:method :spearman :exact false}) :p-value))
      "ct_spearman_ties_p" (let [[a b] (string/split "|" argstr)]
                             ((stats/cor-test (nums a) (nums b) {:method :spearman}) :p-value))
      "ct_spearman_cont_p" (let [[a b] (string/split "|" argstr)]
                             ((stats/cor-test (nums a) (nums b)
                                              {:method :spearman :exact false
                                               :continuity true}) :p-value))
      (errorf "unknown fn %s" fn-kw)))
  (assert-close got ev 1e-9
                (string fn-kw " " (string/slice argstr 0 (min 30 (length argstr))))))

(print "stats smoke ok")

# ---- tier C: classical tests vs recorded R fixtures --------------------------
(each line (parse-fixture (string fixtures-dir "/tier-c-tests.csv"))
  (def fields (string/split "," line))
  (def fn-kw (fields 0))
  (def argstr (string/join (slice fields 1 -2) ","))
  (def ev (scan-number (last fields)))
  (defn nums [s] (map scan-number (string/split " " s)))
  (defn mh-tbls [groups]
    (map (fn [grp]
           (def ps (map scan-number (string/split " " grp)))
           @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]])
         groups))
  (def got
    (case fn-kw
      "t_one" (let [xs (nums argstr)] ((stats/t-test-one xs 5.0 {:impl :pure}) :statistic))
      "t_one_p" (let [xs (nums argstr)] ((stats/t-test-one xs 5.0 {:impl :pure}) :p-value))
      "t_one_df" (let [xs (nums argstr)] ((stats/t-test-one xs 5.0 {:impl :pure}) :parameter))
      "t_one_ci_lo" (let [xs (nums argstr)] (get ((stats/t-test-one xs 5.0 {:impl :pure}) :conf-interval) 0))
      "t_one_ci_hi" (let [xs (nums argstr)] (get ((stats/t-test-one xs 5.0 {:impl :pure}) :conf-interval) 1))
      "t_welch" (let [[a b] (string/split "|" argstr)] ((stats/t-test-two (nums a) (nums b) {:impl :pure}) :statistic))
      "t_welch_p" (let [[a b] (string/split "|" argstr)] ((stats/t-test-two (nums a) (nums b) {:impl :pure}) :p-value))
      "t_welch_df" (let [[a b] (string/split "|" argstr)] ((stats/t-test-two (nums a) (nums b) {:impl :pure}) :parameter))
      "t_pooled" (let [[a b] (string/split "|" argstr)] ((stats/t-test-two (nums a) (nums b) {:var-equal true}) :statistic))
      "t_pooled_p" (let [[a b] (string/split "|" argstr)] ((stats/t-test-two (nums a) (nums b) {:var-equal true}) :p-value))
      "t_pooled_df" (let [[a b] (string/split "|" argstr)] ((stats/t-test-two (nums a) (nums b) {:var-equal true}) :parameter))
      "t_paired" (let [[a b] (string/split "|" argstr)] ((stats/t-test-paired (nums a) (nums b) {:impl :pure}) :statistic))
      "t_paired_p" (let [[a b] (string/split "|" argstr)] ((stats/t-test-paired (nums a) (nums b) {:impl :pure}) :p-value))
      "chisq_gof" (let [xs (nums argstr)] ((stats/chisq-gof xs :uniform {:impl :pure}) :statistic))
      "chisq_gof_p" (let [xs (nums argstr)] ((stats/chisq-gof xs :uniform {:impl :pure}) :p-value))
      "chisq_gof_df" (let [xs (nums argstr)] ((stats/chisq-gof xs :uniform {:impl :pure}) :parameter))
      "chisq_gof_exp" (let [[a b] (string/split "|" argstr)]
                        ((stats/chisq-gof (nums a)
                                          (map |(* $ (/ (reduce + 0 (nums a)) (reduce + 0 (nums b))))
                                               (nums b))
                                          {:impl :pure})
                          :statistic))
      "chisq_gof_exp_p" (let [[a b] (string/split "|" argstr)]
                          ((stats/chisq-gof (nums a)
                                            (map |(* $ (/ (reduce + 0 (nums a)) (reduce + 0 (nums b))))
                                                 (nums b))
                                            {:impl :pure})
                            :p-value))
      "chisq_test" (let [[a b] (string/split "|" argstr)] ((stats/chisq-test @[(nums a) (nums b)]) :statistic))
      "chisq_test_p" (let [[a b] (string/split "|" argstr)] ((stats/chisq-test @[(nums a) (nums b)]) :p-value))
      "chisq_test_df" (let [[a b] (string/split "|" argstr)] ((stats/chisq-test @[(nums a) (nums b)]) :parameter))
      "fisher_p" (let [[a b c d] (string/split "|" argstr)]
                   ((stats/fisher-test-2x2 @[@[(scan-number a) (scan-number b)] @[(scan-number c) (scan-number d)]] {:impl :pure}) :p-value))
      "fisher_odds" (let [[a b c d] (string/split "|" argstr)]
                      ((stats/fisher-test-2x2 @[@[(scan-number a) (scan-number b)] @[(scan-number c) (scan-number d)]] {:impl :pure}) :estimate))
      "wilcoxon_sr_V" (let [xs (nums argstr)] ((stats/wilcoxon-sr xs 0 {:impl :pure}) :statistic))
      "wilcoxon_sr_p" (let [xs (nums argstr)] ((stats/wilcoxon-sr xs 0 {:impl :pure}) :p-value))
      "mw_W" (let [[a b] (string/split "|" argstr)] ((stats/mann-whitney (nums a) (nums b) {:impl :pure}) :statistic))
      "mw_p" (let [[a b] (string/split "|" argstr)] ((stats/mann-whitney (nums a) (nums b) {:impl :pure}) :p-value))
      "kw_H" (let [[a b c] (string/split "|" argstr)]
               ((stats/kruskal-wallis @[(nums a) (nums b) (nums c)]) :statistic))
      "kw_p" (let [[a b c] (string/split "|" argstr)]
               ((stats/kruskal-wallis @[(nums a) (nums b) (nums c)]) :p-value))
      # ---- p.adjust / binom-test / poisson-test / sign-test ----
      "padj_bonferroni" (let [[ps idx] (string/split "|" argstr)]
                          ((stats/p-adjust (nums ps) {:method :bonferroni}) (scan-number idx)))
      "padj_holm" (let [[ps idx] (string/split "|" argstr)]
                    ((stats/p-adjust (nums ps) {:method :holm}) (scan-number idx)))
      "padj_hochberg" (let [[ps idx] (string/split "|" argstr)]
                       ((stats/p-adjust (nums ps) {:method :hochberg}) (scan-number idx)))
      "padj_bh" (let [[ps idx] (string/split "|" argstr)]
                  ((stats/p-adjust (nums ps) {:method :bh}) (scan-number idx)))
      "padj_by" (let [[ps idx] (string/split "|" argstr)]
                 ((stats/p-adjust (nums ps) {:method :by}) (scan-number idx)))
      "padj_bh_n8" (let [[ps idx] (string/split "|" argstr)]
                    ((stats/p-adjust (nums ps) {:method :bh :n 8}) (scan-number idx)))
      "binom_p" (let [ps (map scan-number (string/split "|" argstr))]
                  ((stats/binom-test (ps 0) (ps 1) (ps 2)
                                     (if (= 4 (length ps)) {:conf-level (ps 3)} {})) :p-value))
      "binom_ci_lo" (let [ps (map scan-number (string/split "|" argstr))]
                      (get ((stats/binom-test (ps 0) (ps 1) (ps 2)
                                              (if (= 4 (length ps)) {:conf-level (ps 3)} {}))
                            :conf-interval) 0))
      "binom_ci_hi" (let [ps (map scan-number (string/split "|" argstr))]
                      (get ((stats/binom-test (ps 0) (ps 1) (ps 2)
                                              (if (= 4 (length ps)) {:conf-level (ps 3)} {}))
                            :conf-interval) 1))
      "binom_less_p" (let [ps (map scan-number (string/split "|" argstr))]
                       ((stats/binom-test (ps 0) (ps 1) (ps 2) {:alternative :less}) :p-value))
      "binom_great_p" (let [ps (map scan-number (string/split "|" argstr))]
                        ((stats/binom-test (ps 0) (ps 1) (ps 2) {:alternative :greater}) :p-value))
      "pois_p" (let [ps (map scan-number (string/split "|" argstr))]
                 ((stats/poisson-test (ps 0) (ps 1) (ps 2) {}) :p-value))
      "pois_ci_lo" (let [ps (map scan-number (string/split "|" argstr))]
                     (get ((stats/poisson-test (ps 0) (ps 1) (ps 2) {}) :conf-interval) 0))
      "pois_ci_hi" (let [ps (map scan-number (string/split "|" argstr))]
                    (get ((stats/poisson-test (ps 0) (ps 1) (ps 2) {}) :conf-interval) 1))
      "pois_est" (let [ps (map scan-number (string/split "|" argstr))]
                  ((stats/poisson-test (ps 0) (ps 1) (ps 2) {}) :estimate))
      "pois2_p" (let [ps (map scan-number (string/split "|" argstr))]
                 ((stats/poisson-test @[(ps 0) (ps 1)] @[(ps 2) (ps 3)] (ps 4) {}) :p-value))
      "pois2_ci_lo" (let [ps (map scan-number (string/split "|" argstr))]
                      (get ((stats/poisson-test @[(ps 0) (ps 1)] @[(ps 2) (ps 3)] (ps 4) {})
                            :conf-interval) 0))
      "pois2_ci_hi" (let [ps (map scan-number (string/split "|" argstr))]
                      (get ((stats/poisson-test @[(ps 0) (ps 1)] @[(ps 2) (ps 3)] (ps 4) {})
                            :conf-interval) 1))
      "pois2_est" (let [ps (map scan-number (string/split "|" argstr))]
                   ((stats/poisson-test @[(ps 0) (ps 1)] @[(ps 2) (ps 3)] (ps 4) {}) :estimate))
      "sign_p" (let [[ds mu] (string/split "|" argstr)]
                 ((stats/sign-test (nums ds) {:mu (scan-number mu)}) :p-value))
      "sign_ci_lo" (let [[ds mu] (string/split "|" argstr)]
                     (get ((stats/sign-test (nums ds) {:mu (scan-number mu)})
                           :conf-interval) 0))
      "sign_ci_hi" (let [[ds mu] (string/split "|" argstr)]
                    (get ((stats/sign-test (nums ds) {:mu (scan-number mu)})
                          :conf-interval) 1))
      "sign_nzero" (let [[ds mu] (string/split "|" argstr)]
                     ((stats/sign-test (nums ds) {:mu (scan-number mu)}) :n-zero))
      # ---- prop-test (one/two-proportion) ----
      "prop_stat" (let [ps (map scan-number (string/split "|" argstr))]
                    ((stats/prop-test (ps 0) (ps 1) (ps 2)
                                      {:correct (= 1 (ps 3))}) :statistic))
      "prop_p" (let [ps (map scan-number (string/split "|" argstr))]
                 ((stats/prop-test (ps 0) (ps 1) (ps 2)
                                   {:correct (= 1 (ps 3))}) :p-value))
      "prop_ci_lo" (let [ps (map scan-number (string/split "|" argstr))]
                    (get ((stats/prop-test (ps 0) (ps 1) (ps 2)
                                           {:correct (= 1 (ps 3))})
                          :conf-interval) 0))
      "prop_ci_hi" (let [ps (map scan-number (string/split "|" argstr))]
                    (get ((stats/prop-test (ps 0) (ps 1) (ps 2)
                                           {:correct (= 1 (ps 3))})
                          :conf-interval) 1))
      "prop_less_p" (let [ps (map scan-number (string/split "|" argstr))]
                     ((stats/prop-test (ps 0) (ps 1) (ps 2)
                                       {:alternative :less}) :p-value))
      "prop_great_p" (let [ps (map scan-number (string/split "|" argstr))]
                      ((stats/prop-test (ps 0) (ps 1) (ps 2)
                                        {:alternative :greater}) :p-value))
      "prop_ci_less_hi" (let [ps (map scan-number (string/split "|" argstr))]
                          (get ((stats/prop-test (ps 0) (ps 1) (ps 2)
                                                 {:alternative :less})
                                :conf-interval) 1))
      "prop_ci_great_lo" (let [ps (map scan-number (string/split "|" argstr))]
                           (get ((stats/prop-test (ps 0) (ps 1) (ps 2)
                                                  {:alternative :greater})
                                 :conf-interval) 0))
      "prop_ci90_lo" (let [ps (map scan-number (string/split "|" argstr))]
                      (get ((stats/prop-test (ps 0) (ps 1) (ps 2)
                                             {:conf-level (ps 3)})
                            :conf-interval) 0))
      "prop_ci90_hi" (let [ps (map scan-number (string/split "|" argstr))]
                      (get ((stats/prop-test (ps 0) (ps 1) (ps 2)
                                             {:conf-level (ps 3)})
                            :conf-interval) 1))
      "prop2_stat" (let [ps (string/split "|" argstr)]
                    ((stats/prop-test (nums (ps 0)) (nums (ps 1)) nil
                                      {:correct (= 1 (scan-number (ps 2)))})
                     :statistic))
      "prop2_p" (let [ps (string/split "|" argstr)]
                  ((stats/prop-test (nums (ps 0)) (nums (ps 1)) nil
                                    {:correct (= 1 (scan-number (ps 2)))})
                   :p-value))
      "prop2_ci_lo" (let [ps (string/split "|" argstr)]
                     (get ((stats/prop-test (nums (ps 0)) (nums (ps 1)) nil
                                             {:correct
                                              (= 1 (scan-number (ps 2)))})
                           :conf-interval) 0))
      "prop2_ci_hi" (let [ps (string/split "|" argstr)]
                     (get ((stats/prop-test (nums (ps 0)) (nums (ps 1)) nil
                                             {:correct
                                              (= 1 (scan-number (ps 2)))})
                           :conf-interval) 1))
      "prop2_less_p" (let [ps (string/split "|" argstr)]
                       ((stats/prop-test (nums (ps 0)) (nums (ps 1)) nil
                                         {:alternative :less}) :p-value))
      "prop2_great_p" (let [ps (string/split "|" argstr)]
                        ((stats/prop-test (nums (ps 0)) (nums (ps 1)) nil
                                          {:alternative :greater}) :p-value))
      "prop2_ci_less_hi" (let [ps (string/split "|" argstr)]
                          (get ((stats/prop-test (nums (ps 0)) (nums (ps 1)) nil
                                                  {:alternative :less})
                                :conf-interval) 1))
      "prop2_ci_great_lo" (let [ps (string/split "|" argstr)]
                            (get ((stats/prop-test (nums (ps 0)) (nums (ps 1)) nil
                                                    {:alternative :greater})
                                  :conf-interval) 0))
      "prop2g_stat" (let [ps (string/split "|" argstr)]
                      ((stats/prop-test (nums (ps 0)) (nums (ps 1)) (nums (ps 2))
                                        {:correct (= 1 (scan-number (ps 3)))})
                       :statistic))
      "prop2g_p" (let [ps (string/split "|" argstr)]
                   ((stats/prop-test (nums (ps 0)) (nums (ps 1)) (nums (ps 2))
                                     {:correct (= 1 (scan-number (ps 3)))})
                    :p-value))
      # ---- epi-2x2 (epiR) + cohens-d/hedges-g (effectsize) ----
      "epi_rr_est" (let [ps (map scan-number (string/split "|" argstr))]
                     (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]] {})
                       :rr) :estimate))
      "epi_rr_lo" (let [ps (map scan-number (string/split "|" argstr))]
                    (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]] {})
                           :rr) :conf-interval) 0))
      "epi_rr_hi" (let [ps (map scan-number (string/split "|" argstr))]
                    (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]] {})
                           :rr) :conf-interval) 1))
      "epi_or_est" (let [ps (map scan-number (string/split "|" argstr))]
                     (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]] {})
                       :or) :estimate))
      "epi_or_lo" (let [ps (map scan-number (string/split "|" argstr))]
                    (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]] {})
                           :or) :conf-interval) 0))
      "epi_or_hi" (let [ps (map scan-number (string/split "|" argstr))]
                    (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]] {})
                           :or) :conf-interval) 1))
      "epi_cc_or_est" (let [ps (map scan-number (string/split "|" argstr))]
                        (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]]
                                         {:method :case-control}) :or) :estimate))
      "epi_cc_or_lo" (let [ps (map scan-number (string/split "|" argstr))]
                       (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]]
                                        {:method :case-control}) :or)
                            :conf-interval) 0))
      "epi_cc_or_hi" (let [ps (map scan-number (string/split "|" argstr))]
                       (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]]
                                        {:method :case-control}) :or)
                            :conf-interval) 1))
      "epi_c90_or_lo" (let [ps (map scan-number (string/split "|" argstr))]
                        (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]]
                                         {:conf-level (ps 4)}) :or) :conf-interval) 0))
      "epi_c90_or_hi" (let [ps (map scan-number (string/split "|" argstr))]
                        (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]]
                                         {:conf-level (ps 4)}) :or) :conf-interval) 1))
      "epi_rd_est" (let [ps (map scan-number (string/split "|" argstr))]
                     (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]] {})
                       :rd) :estimate))
      "epi_rd_lo" (let [ps (map scan-number (string/split "|" argstr))]
                    (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]] {})
                           :rd) :conf-interval) 0))
      "epi_rd_hi" (let [ps (map scan-number (string/split "|" argstr))]
                    (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]] {})
                           :rd) :conf-interval) 1))
      "epi_c90_rd_lo" (let [ps (map scan-number (string/split "|" argstr))]
                        (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]]
                                         {:conf-level (ps 4)}) :rd) :conf-interval) 0))
      "epi_c90_rd_hi" (let [ps (map scan-number (string/split "|" argstr))]
                        (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]]
                                         {:conf-level (ps 4)}) :rd) :conf-interval) 1))
      "epi_rdsc_lo" (let [ps (map scan-number (string/split "|" argstr))]
                      (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]] {})
                            :rd-score) :conf-interval) 0))
      "epi_rdsc_hi" (let [ps (map scan-number (string/split "|" argstr))]
                      (get (((stats/epi-2x2 @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]] {})
                            :rd-score) :conf-interval) 1))
      "cd_pool_est" (let [[x y] (string/split "|" argstr)]
                      ((stats/cohens-d (nums x) (nums y) {}) :estimate))
      "cd_pool_lo" (let [[x y] (string/split "|" argstr)]
                     (get ((stats/cohens-d (nums x) (nums y) {}) :conf-interval) 0))
      "cd_pool_hi" (let [[x y] (string/split "|" argstr)]
                     (get ((stats/cohens-d (nums x) (nums y) {}) :conf-interval) 1))
      "cd_unp_est" (let [[x y] (string/split "|" argstr)]
                    ((stats/cohens-d (nums x) (nums y) {:pooled-sd false}) :estimate))
      "cd_unp_lo" (let [[x y] (string/split "|" argstr)]
                    (get ((stats/cohens-d (nums x) (nums y) {:pooled-sd false})
                          :conf-interval) 0))
      "cd_unp_hi" (let [[x y] (string/split "|" argstr)]
                    (get ((stats/cohens-d (nums x) (nums y) {:pooled-sd false})
                          :conf-interval) 1))
      "cd_one_est" (let [x argstr] ((stats/cohens-d (nums x) nil {}) :estimate))
      "cd_one_lo" (let [x argstr]
                    (get ((stats/cohens-d (nums x) nil {}) :conf-interval) 0))
      "cd_one_hi" (let [x argstr]
                    (get ((stats/cohens-d (nums x) nil {}) :conf-interval) 1))
      "cd_paired_est" (let [[x y] (string/split "|" argstr)]
                        ((stats/cohens-d (nums x) (nums y) {:paired true}) :estimate))
      "cd_paired_lo" (let [[x y] (string/split "|" argstr)]
                       (get ((stats/cohens-d (nums x) (nums y) {:paired true})
                             :conf-interval) 0))
      "cd_paired_hi" (let [[x y] (string/split "|" argstr)]
                       (get ((stats/cohens-d (nums x) (nums y) {:paired true})
                             :conf-interval) 1))
      "cd_mu_est" (let [ps (string/split "|" argstr)]
                   ((stats/cohens-d (nums (ps 0)) (nums (ps 1))
                                    {:mu (scan-number (ps 2))}) :estimate))
      "cd_mu_lo" (let [ps (string/split "|" argstr)]
                   (get ((stats/cohens-d (nums (ps 0)) (nums (ps 1))
                                         {:mu (scan-number (ps 2))}) :conf-interval) 0))
      "cd_mu_hi" (let [ps (string/split "|" argstr)]
                   (get ((stats/cohens-d (nums (ps 0)) (nums (ps 1))
                                         {:mu (scan-number (ps 2))}) :conf-interval) 1))
      "cd_c90_est" (let [ps (string/split "|" argstr)]
                     ((stats/cohens-d (nums (ps 0)) (nums (ps 1))
                                      {:conf-level (scan-number (ps 2))}) :estimate))
      "cd_c90_lo" (let [ps (string/split "|" argstr)]
                    (get ((stats/cohens-d (nums (ps 0)) (nums (ps 1))
                                          {:conf-level (scan-number (ps 2))})
                         :conf-interval) 0))
      "cd_c90_hi" (let [ps (string/split "|" argstr)]
                    (get ((stats/cohens-d (nums (ps 0)) (nums (ps 1))
                                          {:conf-level (scan-number (ps 2))})
                         :conf-interval) 1))
      "cd_less_hi" (let [[x y] (string/split "|" argstr)]
                    (get ((stats/cohens-d (nums x) (nums y)
                                          {:alternative :less}) :conf-interval) 1))
      "cd_great_lo" (let [[x y] (string/split "|" argstr)]
                      (get ((stats/cohens-d (nums x) (nums y)
                                            {:alternative :greater}) :conf-interval) 0))
      "hg_pool_est" (let [[x y] (string/split "|" argstr)]
                      ((stats/hedges-g (nums x) (nums y) {}) :estimate))
      "hg_pool_lo" (let [[x y] (string/split "|" argstr)]
                     (get ((stats/hedges-g (nums x) (nums y) {}) :conf-interval) 0))
      "hg_pool_hi" (let [[x y] (string/split "|" argstr)]
                     (get ((stats/hedges-g (nums x) (nums y) {}) :conf-interval) 1))
      "hg_unp_est" (let [[x y] (string/split "|" argstr)]
                     ((stats/hedges-g (nums x) (nums y) {:pooled-sd false}) :estimate))
      "hg_unp_lo" (let [[x y] (string/split "|" argstr)]
                    (get ((stats/hedges-g (nums x) (nums y) {:pooled-sd false})
                          :conf-interval) 0))
      "hg_unp_hi" (let [[x y] (string/split "|" argstr)]
                    (get ((stats/hedges-g (nums x) (nums y) {:pooled-sd false})
                          :conf-interval) 1))
      "hg_one_est" (let [x argstr] ((stats/hedges-g (nums x) nil {}) :estimate))
      "hg_one_lo" (let [x argstr]
                    (get ((stats/hedges-g (nums x) nil {}) :conf-interval) 0))
      "hg_one_hi" (let [x argstr]
                    (get ((stats/hedges-g (nums x) nil {}) :conf-interval) 1))
      "hg_small_est" (let [[x y] (string/split "|" argstr)]
                       ((stats/hedges-g (nums x) (nums y) {}) :estimate))
      "hg_small_lo" (let [[x y] (string/split "|" argstr)]
                      (get ((stats/hedges-g (nums x) (nums y) {}) :conf-interval) 0))
      "hg_small_hi" (let [[x y] (string/split "|" argstr)]
                      (get ((stats/hedges-g (nums x) (nums y) {}) :conf-interval) 1))
      # ---- mcnemar-test + mantelhaen-test ----
      "mc_stat" (let [parts (string/split "|" argstr)
                      ps (map scan-number (string/split " " (parts 0)))]
                  ((stats/mcnemar-test @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]]
                                       {:correct (= 1 (scan-number (parts 1)))})
                   :statistic))
      "mc_p" (let [parts (string/split "|" argstr)
                   ps (map scan-number (string/split " " (parts 0)))]
               ((stats/mcnemar-test @[@[(ps 0) (ps 1)] @[(ps 2) (ps 3)]]
                                    {:correct (= 1 (scan-number (parts 1)))})
                :p-value))
      "mh_stat" (let [tbls (string/split "|" argstr)]
                  ((stats/mantelhaen-test (mh-tbls tbls) {}) :statistic))
      "mh_p" (let [tbls (string/split "|" argstr)]
               ((stats/mantelhaen-test (mh-tbls tbls) {}) :p-value))
      "mh_nocc_stat" (let [tbls (string/split "|" argstr)]
                       ((stats/mantelhaen-test (mh-tbls tbls) {:correct false}) :statistic))
      "mh_nocc_p" (let [tbls (string/split "|" argstr)]
                    ((stats/mantelhaen-test (mh-tbls tbls) {:correct false}) :p-value))
      "mh_great_p" (let [tbls (string/split "|" argstr)]
                     ((stats/mantelhaen-test (mh-tbls tbls) {:alternative :greater}) :p-value))
      "mh_great_lo" (let [tbls (string/split "|" argstr)]
                      (get ((stats/mantelhaen-test (mh-tbls tbls)
                                                   {:alternative :greater})
                           :conf-interval) 0))
      "mh_less_p" (let [tbls (string/split "|" argstr)]
                    ((stats/mantelhaen-test (mh-tbls tbls) {:alternative :less}) :p-value))
      "mh_less_hi" (let [tbls (string/split "|" argstr)]
                     (get ((stats/mantelhaen-test (mh-tbls tbls)
                                                  {:alternative :less})
                          :conf-interval) 1))
      "mh_or" (let [tbls (string/split "|" argstr)]
                ((stats/mantelhaen-test (mh-tbls tbls) {}) :estimate))
      "mh_or_lo" (let [tbls (string/split "|" argstr)]
                   (get ((stats/mantelhaen-test (mh-tbls tbls) {}) :conf-interval) 0))
      "mh_or_hi" (let [tbls (string/split "|" argstr)]
                   (get ((stats/mantelhaen-test (mh-tbls tbls) {}) :conf-interval) 1))
      "mh_c90_lo" (let [parts (string/split "|" argstr)]
                    (get ((stats/mantelhaen-test
                           (mh-tbls (slice parts 0 (- (length parts) 1)))
                           {:conf-level (scan-number (last parts))})
                         :conf-interval) 0))
      "mh_c90_hi" (let [parts (string/split "|" argstr)]
                    (get ((stats/mantelhaen-test
                           (mh-tbls (slice parts 0 (- (length parts) 1)))
                           {:conf-level (scan-number (last parts))})
                         :conf-interval) 1))
      "mh_rr" (let [tbls (string/split "|" argstr)]
                (((stats/mantelhaen-test (mh-tbls tbls) {}) :rr) :estimate))
      "mh_rr_lo" (let [tbls (string/split "|" argstr)]
                   (get (((stats/mantelhaen-test (mh-tbls tbls) {}) :rr) :conf-interval) 0))
      "mh_rr_hi" (let [tbls (string/split "|" argstr)]
                   (get (((stats/mantelhaen-test (mh-tbls tbls) {}) :rr) :conf-interval) 1))
      "mh_homog_rr" (let [tbls (string/split "|" argstr)]
                      ((((stats/mantelhaen-test (mh-tbls tbls) {}) :homogeneity) :rr)
                       :statistic))
      "mh_homog_rr_p" (let [tbls (string/split "|" argstr)]
                        ((((stats/mantelhaen-test (mh-tbls tbls) {}) :homogeneity) :rr)
                         :p-value))
      "mh_homog_or" (let [tbls (string/split "|" argstr)]
                      ((((stats/mantelhaen-test (mh-tbls tbls) {}) :homogeneity) :or)
                       :statistic))
      "mh_homog_or_p" (let [tbls (string/split "|" argstr)]
                        ((((stats/mantelhaen-test (mh-tbls tbls) {}) :homogeneity) :or)
                         :p-value))
      (errorf "unknown fn %s" fn-kw)))
  # fisher_odds: R's $estimate is its own uniroot approximation of the
  # conditional MLE (ours solves E=a exactly; verified R's is looser).
  # cd_*/hg_* CI bounds: effectsize solves the ncp via optim (abstol
  # 1e-9) - ITS bounds carry ~1e-8 solver noise; ours is the exact
  # monotone inversion. ABS tolerance (bounds can sit near 0 where
  # the rel ratio explodes)
  (def ci-noise? (and (or (string/has-prefix? "cd_" fn-kw)
                          (string/has-prefix? "hg_" fn-kw))
                      (or (string/has-suffix? "_lo" fn-kw)
                          (string/has-suffix? "_hi" fn-kw))))
  (if ci-noise?
    (unless (<= (math/abs (- got ev)) 1e-6)
      (errorf "FAIL %s: %.17g != %.17g (abs-noise lane)" fn-kw got ev))
    (assert-close got ev
                  (if (= fn-kw "fisher_odds") 1e-3 1e-9)
                  (string fn-kw " " (string/slice argstr 0 (min 30 (length argstr)))))))

# ---- tier D: regression/survival vs recorded R fixtures ----------------------
# Data files hold the design at %.17g; the results CSV maps fn -> expected.
# Tolerance 1e-8 rel: cox z/p bounded by R coxph's own stopping criterion
# (worst verified row cox_efron_p 1.8e-9); all others <= 4.3e-11.
(defn load-num-csv [path]
  (def out @[])
  (each line (string/split "\n" (string (slurp path)))
    (when (and (not (empty? line)) (not (string/find "\"" line)))
      (array/push out (map scan-number (string/split "," line)))))
  out)

(def dsets @{})
(defn dset [name]
  (unless (get dsets name)
    (put dsets name (load-num-csv (string fixtures-dir "/tier-d-" name ".csv"))))
  (get dsets name))

(def dfits @{})
(defn d-groups [name]
  # aov datasets: [g y] rows, g 1-based -> k sample arrays
  (def rows (dset name))
  (def k (max ;(map |(get $ 0) rows)))
  (def gs (seq [_ :range [0 k]] @[]))
  (each r rows (array/push (gs (- (get r 0) 1)) (get r 1)))
  gs)
(defn run-dfit [algo name]
  (def key (string algo "/" name))
  (unless (get dfits key)
    (put dfits key
         (case algo
           "lm" (let [r (dset name)]
                  (stats/lm (seq [row :in r] @[1.0 (get row 0) (get row 1)])
                            (seq [row :in r] (get row 2)) {:impl :pure}))
           "glm_bin" (let [r (dset name)]
                       (stats/glm (seq [row :in r] @[1.0 (get row 0) (get row 1)])
                                  (seq [row :in r] (get row 2)) :binomial {:impl :pure}))
           "glm_pois" (let [r (dset name)]
                        (stats/glm (seq [row :in r] @[1.0 (get row 0) (get row 1)])
                                   (seq [row :in r] (get row 2)) :poisson {:impl :pure}))
           "cox_efron" (let [r (dset name)]
                         (stats/coxph (seq [row :in r] @[(get row 0) (get row 1)])
                                      (seq [row :in r] (get row 2))
                                      (seq [row :in r] (get row 3)) {:ties :efron}))
           "cox_breslow" (let [r (dset name)]
                           (stats/coxph (seq [row :in r] @[(get row 0) (get row 1)])
                                        (seq [row :in r] (get row 2))
                                        (seq [row :in r] (get row 3)) {:ties :breslow}))
           "clogit" (let [r (dset name)]
                      (stats/clogit (seq [row :in r] @[(get row 1) (get row 2)])
                                    (seq [row :in r] (get row 3))
                                    (seq [row :in r] (get row 0)) {:impl :pure}))
           "km" (let [r (dset name)]
                  (stats/survfit-km (seq [row :in r] (get row 0))
                                    (seq [row :in r] (get row 1)) {}))
           "lr" (let [r (dset name)]
                  (stats/survdiff-logrank (seq [row :in r] (get row 0))
                                          (seq [row :in r] (get row 1))
                                          (seq [row :in r] (get row 2)) {}))
           "aov" (stats/aov (d-groups name) {:impl :pure})
           "oneway" (stats/oneway-test (d-groups name) {:impl :pure})
           "oneway_e" (stats/oneway-test (d-groups name)
                                         {:impl :pure :var-equal true})
           "tukey" (stats/tukey-hsd (run-dfit "aov" name) {:impl :pure})
           "tukey_c9" (stats/tukey-hsd (run-dfit "aov" name)
                                       {:impl :pure :conf-level 0.9})
           "tukey_o" (stats/tukey-hsd (run-dfit "aov" name)
                                      {:impl :pure :ordered true})
           (errorf "unknown algo %s" algo))))
  (get dfits key))

(each line (parse-fixture (string fixtures-dir "/tier-d-tests.csv"))
  (def fields (string/split "," line))
  (def fn-kw (fields 0))
  (def name (fields 1))
  (def ev (scan-number (get fields 2)))
  (var algo nil)
  (var metric nil)
  (each a ["lm" "glm_bin" "glm_pois" "cox_efron" "cox_breslow" "clogit"
           "km" "lr" "aov" "oneway_e" "oneway" "tukey_c9" "tukey_o" "tukey"]
    (when (and (not algo) (string/has-prefix? (string a "_") fn-kw))
      (set algo a)
      (set metric (string/slice fn-kw (+ 1 (length a))))))
  (def fit (run-dfit algo name))
  (defn idx-val [len] (scan-number (string/slice metric len)))
  (def got
    (cond
      (= algo "aov")
      (get-in fit
              (case metric
                "ssg" [:table :group :ss]
                "ssr" [:table :residual :ss]
                "msg" [:table :group :ms]
                "msr" [:table :residual :ms]
                "f" [:table :group :f]
                "p" [:table :group :p]
                "mse" [:mse]
                (errorf "unknown aov metric %s in %s" metric fn-kw)))
      (or (= algo "oneway") (= algo "oneway_e"))
      (case metric
        "f" (fit :statistic)
        "df2" (get (fit :parameter) 1)
        "p" (fit :p-value)
        (errorf "unknown oneway metric %s in %s" metric fn-kw))
      (or (= algo "tukey") (= algo "tukey_c9") (= algo "tukey_o"))
      (let [parts (string/split "_" metric)
            col (find |(or (= $ "diff") (= $ "lwr") (= $ "upr") (= $ "p")) parts)
            i (scan-number (get parts (- (length parts) 2)))
            j (scan-number (get parts (- (length parts) 1)))
            pair (find |(and (= (get $ :i) i) (= (get $ :j) j)) (fit :pairs))]
        (unless pair (errorf "tukey pair %d-%d not found in %s" i j fn-kw))
        (case col
          "diff" (pair :diff)
          "lwr" (pair :lwr)
          "upr" (pair :upr)
          "p" (pair :p)
          (errorf "unknown tukey column %s in %s" col fn-kw)))
      (= algo "km")
      (cond
        (= metric "n") (fit :n)
        (string/has-prefix? "time_" metric) (get (fit :time) (idx-val 5))
        (string/has-prefix? "surv_" metric) (get (fit :surv) (idx-val 5))
        (string/has-prefix? "nrisk_" metric) (get (fit :n-risk) (idx-val 6))
        (string/has-prefix? "nevent_" metric) (get (fit :n-event) (idx-val 7))
        (string/has-prefix? "ncensor_" metric) (get (fit :n-censor) (idx-val 8))
        (errorf "unknown km metric %s in %s" metric fn-kw))
      (= algo "lr")
      (cond
        (= metric "chisq") (fit :statistic)
        (= metric "p") (fit :p-value)
        (= metric "var") (fit :variance)
        (string/has-prefix? "obs_" metric) (get (fit :observed) (idx-val 4))
        (string/has-prefix? "exp_" metric) (get (fit :expected) (idx-val 4))
        (string/has-prefix? "n_" metric) (get (fit :n) (idx-val 2))
        (errorf "unknown lr metric %s in %s" metric fn-kw))
      (cond
        (string/has-prefix? "coef_" metric) (get (fit :coefficients) (idx-val 5))
        (string/has-prefix? "se_" metric) (get (fit :se) (idx-val 3))
        (string/has-prefix? "t_" metric) (get (fit :t) (idx-val 2))
        (string/has-prefix? "z_" metric) (get (fit :z) (idx-val 2))
        (string/has-prefix? "p_" metric) (get (fit :p) (idx-val 2))
        (= metric "sigma") (fit :sigma)
        (= metric "r2") (fit :r2)
        (= metric "adj_r2") (fit :adj-r2)
        (= metric "f") (fit :f)
        (= metric "deviance") (fit :deviance)
        (= metric "loglik") (fit :loglik)
        (errorf "unknown metric %s in %s" metric fn-kw))))
  # the tukey CI rows carry qtukey's secant eps (the tier-b q-fn class)
  (def tol (if (and (string/has-prefix? "tukey" fn-kw)
                    (or (string/find "lwr" fn-kw) (string/find "upr" fn-kw)))
            1e-5 1e-8))
  (assert-close got ev tol (string fn-kw " " name)))

# ---- tier E: resampling vs recorded R fixtures -------------------------------
# Same rows as the conformance suite's tier-E coverage but
# inline: smokes stay zero-dep per the charter TEST LIBRARY DECISION.
(defn load-e-nums [path]
  (def out @[])
  (each line (string/split "\n" (string (slurp path)))
    (when (and (not (empty? line)) (not (string/find "\"" line)))
      (array/push out (map scan-number (string/split "," line)))))
  out)

(defn e-exp [fn-kw]
  (var ev nil)
  (each line (string/split "\n" (string (slurp (string fixtures-dir "/tier-e-tests.csv"))))
    (when (and (not (empty? line)) (not (string/find "fn," line)))
      (def fields (string/split "," line))
      (when (= (fields 0) fn-kw) (set ev (scan-number (get fields 2))))))
  (unless ev (errorf "no tier-e fixture row for %s" fn-kw))
  ev)

(defn e-col [path]
  # single-column file -> flat number array (load-e-nums keeps row shape)
  (map |(get $ 0) (load-e-nums path)))

(def e-mean-i (fn [data idx] (var s 0) (each i idx (+= s (get data i))) (/ s (* 1.0 (length idx)))))
(def e-median-i
  (fn [data idx]
    (def vs (sorted (map |(get data $) idx)))
    (def m (length vs))
    (if (odd? m) (get vs (div m 2))
      (/ (+ (get vs (dec (div m 2))) (get vs (div m 2))) 2))))

(def e-data-mean (e-col (string fixtures-dir "/tier-e-boot-mean-data.csv")))
(def e-reps-mean (e-col (string fixtures-dir "/tier-e-boot-mean-reps.csv")))
(def e-boot-mean {:t0 (e-exp "boot_mean_t0") :t e-reps-mean :B (length e-reps-mean)})
(def e-l-mean ((stats/jackknife e-data-mean e-mean-i) :L))
(assert-close ((stats/boot-ci e-boot-mean {:impl :pure}) 0) (e-exp "boot_mean_perc_lo") 1e-9 "boot mean perc lo")
(assert-close ((stats/boot-ci e-boot-mean {:impl :pure}) 1) (e-exp "boot_mean_perc_hi") 1e-9 "boot mean perc hi")
(assert-close ((stats/boot-ci e-boot-mean {:type :bca :L e-l-mean}) 0) (e-exp "boot_mean_bca_lo") 1e-9 "boot mean bca lo")
(assert-close ((stats/boot-ci e-boot-mean {:type :bca :L e-l-mean}) 1) (e-exp "boot_mean_bca_hi") 1e-9 "boot mean bca hi")

(def e-data-med (e-col (string fixtures-dir "/tier-e-boot-median-data.csv")))
(def e-reps-med (e-col (string fixtures-dir "/tier-e-boot-median-reps.csv")))
(def e-boot-med {:t0 (e-exp "boot_median_t0") :t e-reps-med :B (length e-reps-med)})
(def e-l-med ((stats/jackknife e-data-med e-median-i) :L))
(assert-close ((stats/boot-ci e-boot-med {:impl :pure}) 0) (e-exp "boot_median_perc_lo") 1e-9 "boot median perc lo")
(assert-close ((stats/boot-ci e-boot-med {:impl :pure}) 1) (e-exp "boot_median_perc_hi") 1e-9 "boot median perc hi")
(assert-close ((stats/boot-ci e-boot-med {:type :bca :L e-l-med}) 0) (e-exp "boot_median_bca_lo") 1e-9 "boot median bca lo")
(assert-close ((stats/boot-ci e-boot-med {:type :bca :L e-l-med}) 1) (e-exp "boot_median_bca_hi") 1e-9 "boot median bca hi")

(defn e-split [name]
  (def rs (load-e-nums (string fixtures-dir "/tier-e-" name ".csv")))
  [(map |(get $ 1) (filter |(= (get $ 0) 0) rs))
   (map |(get $ 1) (filter |(= (get $ 0) 1) rs))])
(def [e-p1a e-p1b] (e-split "perm1"))
(def [e-p2a e-p2b] (e-split "perm2"))
(def e-t1 (stats/perm-test-two e-p1a e-p1b {:impl :pure}))
(assert (= (e-t1 :B) 126) "perm1 C(9,4)")
(assert-close (e-t1 :p-value) (e-exp "perm1_two") 1e-9 "perm1 two-sided p")
(assert-close ((stats/perm-test-two e-p1a e-p1b {:alternative :less}) :p-value)
              (e-exp "perm1_less") 1e-9 "perm1 less p")
(def e-t2 (stats/perm-test-two e-p2a e-p2b {:impl :pure}))
(assert (= (e-t2 :B) 462) "perm2 C(11,5)")
(assert-close (e-t2 :p-value) (e-exp "perm2_two") 1e-9 "perm2 two-sided p")
(assert-close ((stats/perm-test-two e-p2a e-p2b {:alternative :greater}) :p-value)
              (e-exp "perm2_greater") 1e-9 "perm2 greater p")

# RNG paths are never compared across languages: determinism + MC band
(def e-o1 (stats/bootstrap e-data-mean e-mean-i {:B 400 :seed "replay"}))
(def e-o2 (stats/bootstrap e-data-mean e-mean-i {:B 400 :seed "replay"}))
(assert (deep= (e-o1 :t) (e-o2 :t)) "same seed identical replicates")
(def e-m1 (stats/perm-test-two e-p1a e-p1b {:exact false :B 2000 :seed "replay"}))
(def e-m2 (stats/perm-test-two e-p1a e-p1b {:exact false :B 2000 :seed "replay"}))
(assert (= (e-m1 :p-value) (e-m2 :p-value)) "same seed identical MC p")
(def e-band (* 4 (math/sqrt (/ (* (e-t1 :p-value) (- 1 (e-t1 :p-value))) 2000))))
(assert (<= (math/abs (- (e-m1 :p-value) (e-t1 :p-value))) e-band) "MC p in 4-sigma band")

# ---- tier F: bayesian v1 vs recorded R fixtures -------------------------------
# Same rows as the conformance suite's tier-F coverage, inline.
(defn f-exp [fn-kw]
  (var ev nil)
  (each line (string/split "\n" (string (slurp (string fixtures-dir "/tier-f-tests.csv"))))
    (when (and (not ev) (not (empty? line)) (not (string/find "fn," line)))
      (def fields (string/split "," line))
      (when (= (fields 0) fn-kw) (set ev (scan-number (get fields 2))))))
  (unless ev (errorf "no tier-f fixture row for %s" fn-kw))
  ev)

(defn f-args [fn-kw]
  (var out nil)
  (each line (string/split "\n" (string (slurp (string fixtures-dir "/tier-f-tests.csv"))))
    (when (and (not out) (not (empty? line)) (not (string/find "fn," line)))
      (def fields (string/split "," line))
      (when (= (fields 0) fn-kw) (set out (get fields 1)))))
  (unless out (errorf "no tier-f fixture row for %s" fn-kw))
  out)

(defn f-parts [fn-kw] (map scan-number (string/split "|" (f-args fn-kw))))
(defn f-close [got fn-kw tol]
  (def ev (f-exp fn-kw))
  (def diff (math/abs (- got ev)))
  (def rel (if (zero? ev) diff (/ diff (math/abs ev))))
  (assert (or (<= rel tol) (<= diff 1e-12))
          (string "tier-f " fn-kw ": " got " vs " ev)))

# dist plumbing (beta + gamma, R-exact)
(each [fn-kw call] [["dbeta" (fn [x a b] (stats/dbeta x a b {:impl :pure}))]
                    ["pbeta" (fn [x a b] (stats/pbeta x a b {:impl :pure}))]
                    ["qbeta" (fn [p a b] (stats/qbeta p a b {:impl :pure}))]
                    ["dgamma" (fn [x s r] (stats/dgamma x s r {:impl :pure}))]
                    ["pgamma" (fn [x s r] (stats/pgamma x s r {:impl :pure}))]
                    ["qgamma" (fn [p s r] (stats/qgamma p s r {:impl :pure}))]]
  (def [a1 b1 c1] (f-parts fn-kw))
  (f-close (call a1 b1 c1) fn-kw 1e-9))

# conjugate posteriors (analytic)
(def [bk bn ba bb] (f-parts "post_bb_mean"))
(def bbpost (stats/post-beta-binomial bk bn ba bb {:impl :pure}))
(f-close (bbpost :mean) "post_bb_mean" 1e-9)
(f-close ((bbpost :ci) 0) "post_bb_ci_lo" 1e-9)
(f-close ((bbpost :ci) 1) "post_bb_ci_hi" 1e-9)
(f-close (bbpost :log-marginal) "post_bb_logmarg" 1e-9)
(def nnsegs (string/split "|" (f-args "post_nn_mean")))
(def nxs (map scan-number (string/split " " (get nnsegs 0))))
(def [nmu nt ns] (map scan-number (slice nnsegs 1)))
(def nnpost (stats/post-normal-normal nxs nmu nt ns {:impl :pure}))
(f-close (nnpost :mean) "post_nn_mean" 1e-9)
(f-close (nnpost :sd) "post_nn_sd" 1e-9)
(f-close (nnpost :log-marginal) "post_nn_logmarg" 1e-9)
(def gpsegs (string/split "|" (f-args "post_gp_shape")))
(def gks (map scan-number (string/split " " (get gpsegs 0))))
(def gex (map scan-number (string/split " " (get gpsegs 1))))
(def [gpa gpb] (map scan-number (slice gpsegs 2)))
(def gppost (stats/post-gamma-poisson gks gex gpa gpb {:impl :pure}))
(f-close (gppost :shape) "post_gp_shape" 1e-9)
(f-close (gppost :rate) "post_gp_rate" 1e-9)
(f-close ((gppost :ci) 0) "post_gp_ci_lo" 1e-9)
(f-close ((gppost :ci) 1) "post_gp_ci_hi" 1e-9)

# bayes-lm flat-prior limit = lm (tier-d lm1 data; R-verified there)
(def f-rows (load-e-nums (string fixtures-dir "/tier-d-lm1.csv")))
(def f-X (seq [r :in f-rows] @[1.0 (get r 0) (get r 1)]))
(def f-y (seq [r :in f-rows] (get r 2)))
(def f-blm (stats/bayes-lm f-X f-y {:impl :pure}))
(def f-lm (stats/lm f-X f-y {:impl :pure}))
(for j 0 3
  (assert (<= (math/abs (- (get (f-blm :beta-mean) j) (get (f-lm :coefficients) j))) 1e-12)
          (string "bayes-lm coef " j " = lm"))
  (assert (<= (math/abs (- (get (f-blm :se) j) (get (f-lm :se) j))) 1e-12)
          (string "bayes-lm se " j " = lm"))
  (assert (<= (math/abs (- (get (f-blm :p) j) (get (f-lm :p) j))) 1e-12)
          (string "bayes-lm p " j " = lm")))
(assert (= (f-blm :df) (f-lm :df)) "bayes-lm df = lm df")

# Nelder-Mead vs R optim optimum
(def f-f1 (fn [x] (+ (* (- (x 0) 2.5) (- (x 0) 2.5))
                     (* 3 (- (x 1) -1.25) (- (x 1) -1.25)) 0.7)))
(def f-nm (stats/nelder-mead f-f1 @[0 0] {:tol 1e-13 :max-it 20000}))
(f-close (get (f-nm :x) 0) "nm1_x1" 1e-4)
(f-close (get (f-nm :x) 1) "nm1_x2" 1e-4)

# EB beta-binomial (MLE vs R optim; loglik tight, params flat-surface 1e-4)
(def f-ebsegs (string/split "|" (f-args "eb_loglik")))
(def f-eks (map scan-number (string/split " " (get f-ebsegs 0))))
(def f-ens (map scan-number (string/split " " (get f-ebsegs 1))))
(def f-eb (stats/eb-beta-binomial f-eks f-ens {:impl :pure}))
(f-close (f-eb :log-marginal) "eb_loglik" 1e-9)
(f-close (f-eb :a) "eb_a" 1e-4)
(f-close (f-eb :b) "eb_b" 1e-4)
(def f-eg (get (f-eb :groups) 3))
(f-close (f-eg :mean) "eb_g4_mean" 1e-3)
(f-close ((f-eg :ci) 0) "eb_g4_ci_lo" 1e-3)
(f-close ((f-eg :ci) 1) "eb_g4_ci_hi" 1e-3)

# simulation + convergence: known truth, no RNG-path compares
(def f-lp (fn [x] (* -0.5 x x)))
(def f-mh (stats/mh f-lp 0 {:steps 8000 :burn 2000 :seed "sm"}))
(def f-ms (stats/draws-summary (f-mh :draws) {:impl :pure}))
(assert (< (math/abs (f-ms :mean)) 0.07) "MH N(0,1) mean band")
(assert (and (> (f-ms :sd) 0.92) (< (f-ms :sd) 1.08)) "MH N(0,1) sd band")
(def f-mh2 (stats/mh f-lp 0 {:steps 8000 :burn 2000 :seed "sm"}))
(assert (deep= (f-mh :draws) (f-mh2 :draws)) "MH same-seed determinism")

(def f-gb (stats/gibbs-normal @[5.1 4.9 5.2 5.0 4.8 5.05]
                              {:steps 8000 :burn 2000 :seed "sm"
                               :m0 3 :k0 1 :a0 4 :b0 2}))
# grid-integration posterior: E[mu]=4.7215 sd(mu)=0.3141 E[s2]=0.6915
(assert (< (math/abs (- (f-gb :mu-mean) 4.7215)) 0.04) "gibbs E[mu] band")
(assert (< (math/abs (- (f-gb :mu-sd) 0.3141)) 0.04) "gibbs sd(mu) band")
(assert (< (math/abs (- (f-gb :s2-mean) 0.6915)) 0.06) "gibbs E[s2] band")

(def f-iid (seq [c :range [0 4]]
              (let [r (math/rng (string "smiid" c))]
                (seq [_ :range [0 400]] (stats/qnorm (math/rng-uniform r) 0 1 {:impl :pure})))))
(assert (< (math/abs (- (stats/rhat f-iid {:impl :pure}) 1)) 0.01) "iid rhat ~ 1")
(assert (> (stats/ess f-iid {:impl :pure}) 1200) "iid ess near m*n")
(def f-stuck @[(seq [_ :range [0 300]] 1.0) (seq [_ :range [0 300]] 3.0)])
(assert (> (stats/rhat f-stuck {:impl :pure}) 1.5) "diverged chains rhat >> 1")

(def f-hr (math/rng "smhpd"))
(def f-hd (stats/rnorm f-hr 8000 0 1))
(def f-h (stats/hpd f-hd 0.95))
(assert (< (math/abs (+ (f-h 0) 1.96)) 0.15) "hpd lo ~ -1.96")
(assert (< (math/abs (- (f-h 1) 1.96)) 0.15) "hpd hi ~ 1.96")

# ---- tier G: multivariate vs recorded R fixtures (FINAL tier) -----------------
(defn g-rows []
  (filter |(and (not (empty? $)) (not (string/find "fn," $)))
          (string/split "\n" (string (slurp (string fixtures-dir "/tier-g-tests.csv"))))))
(defn g-exp [fn-kw]
  (var ev nil)
  (each line (g-rows)
    (def fields (string/split "," line))
    (when (= (fields 0) fn-kw)
      (set ev (scan-number (last fields)))))
  (unless ev (errorf "no tier-g fixture row for %s" fn-kw))
  ev)
(defn g-mat [name]
  (def out @[])
  (each line (string/split "\n" (string (slurp (string fixtures-dir "/tier-g-" name ".csv"))))
    (unless (empty? line)
      (array/push out (map scan-number (string/split "," line)))))
  out)
(defn g-close [got fn-kw tol]
  (def ev (g-exp fn-kw))
  (def diff (math/abs (- got ev)))
  (def rel (if (zero? ev) diff (/ diff (math/abs ev))))
  (assert (or (<= rel tol) (<= diff 1e-12)) (string "tier-g " fn-kw ": " got " vs " ev)))

# PCA vs prcomp (center-only + scaled)
(def g-pca1 (stats/pca (g-mat "pca1") {:impl :pure}))
(for j 0 3
  (g-close (get (g-pca1 :sdev) j) (string "pca1_sdev_" j) 1e-9)
  (g-close (get (g-pca1 :var-explained) j) (string "pca1_varexp_" j) 1e-9))
(for i 0 3 (g-close (math/abs (get (get (g-pca1 :rotation) i) 0))
                    (string "pca1_rot_0_" i) 1e-9))
(def g-pca2 (stats/pca (g-mat "pca2") {:scale true}))
(for j 0 4 (g-close (get (g-pca2 :sdev) j) (string "pca2_sdev_" j) 1e-9))
(g-close (get (g-pca2 :center) 1) "pca2_center_1" 1e-9)

# k-means vs R Lloyd, same explicit start
(def g-km1x (g-mat "km1"))
(def g-km1 (stats/k-means g-km1x 3 {:centers @[@[0 0] @[5 5] @[9 1]] :iters 50}))
(g-close (g-km1 :tot-withinss) "km1_tot" 1e-9)
(for c 0 3 (g-close (get (g-km1 :withinss) c) (string "km1_within_" c) 1e-9))
(assert (deep= (g-km1 :labels) @[0 0 0 0 1 1 1 1 2 2 2 2]) "km1 labels vs R")

# seeded invariants: determinism + Lloyd-fixed point
(def g-a (stats/k-means g-km1x 3 {:seed "sm" :iters 50}))
(def g-b (stats/k-means g-km1x 3 {:seed "sm" :iters 50}))
(assert (deep= (g-a :labels) (g-b :labels)) "k-means same-seed determinism")
(def g-c2 (stats/k-means g-km1x 3 {:centers (g-a :centers) :iters 50}))
(assert (deep= (g-c2 :labels) (g-a :labels)) "converged result is a Lloyd fixed point")

# ---- ECDF + KS - recorded fixtures + inline spot checks -----
# the ecdf-ks.qnt contract inline (spec pair D = 2/3; R-parity p = 0.4)
(def cdf-fn (stats/ecdf @[1 1 2 3 3 3] {}))
(assert (= 0 (cdf-fn 0.5)) "ecdf below data")
(assert (<= (math/abs (- (cdf-fn 1) (/ 2 6))) 1e-15) "ecdf tie block shares one value")
(assert (= 1 (cdf-fn 3)) "ecdf F(max) = 1")
(def ks-sp (stats/ks-test @[1 1 2] @[1 3 3] {}))
(assert (<= (math/abs (- (ks-sp :statistic) (/ 2 3))) 1e-15) "spec pair D = 2/3")
(assert (<= (math/abs (- (ks-sp :p-value) 0.4)) 1e-12) "spec pair p = 0.4 (R parity)")
(assert (ks-sp :exact) "small samples -> the exact lane")
(assert (ks-sp :ties) "tie pools flagged")
# default :impl + both-impl parity on the walk (S1.1)
(def ks-par
  (let [xs (seq [i :range [0 30]] (math/sin (* i 0.3)))
        ys (seq [i :range [0 25]] (+ 0.5 (math/cos (* i 0.21))))]
    [(stats/ks-test xs ys {:impl :pure}) (stats/ks-test xs ys {:impl :num})]))
(assert (deep= ((ks-par 0) :statistic) ((ks-par 1) :statistic))
        "ks two-sample D: pure == num walk")
(assert (= ((ks-par 0) :p-value) ((ks-par 1) :p-value)) "ks p: pure == num")
# one-sample vs the named CDF, default resolution
(def ks-one (stats/ks-test (seq [i :range [0 40]] (math/sin (* i 0.2))) :pnorm {}))
(assert (and (>= (ks-one :p-value) 0) (<= (ks-one :p-value) 1)) "one-sample p in [0,1]")
# the full fixture row sweep (pure impl - the harness covers both)
(each line (parse-fixture (string fixtures-dir "/tier-h6-tests.csv"))
  (def fields (string/split "," line))
  (def fn-kw (fields 0))
  (def argstr (string/join (slice fields 1 -2) ","))
  (def ev (scan-number (last fields)))
  (defn nums [s] (map scan-number (string/split " " s)))
  (def got
    (case fn-kw
      "ecdf_val" (let [[data q] (string/split "|" argstr)]
                   ((stats/ecdf (nums data) {:impl :pure}) (scan-number q)))
      "ks2_d" (let [[a b] (string/split "|" argstr)]
                ((stats/ks-test (nums a) (nums b) {:impl :pure}) :statistic))
      "ks2_p" (let [[a b] (string/split "|" argstr)]
                ((stats/ks-test (nums a) (nums b) {:impl :pure}) :p-value))
      "ks2_exact" (let [[a b] (string/split "|" argstr)]
                    (if ((stats/ks-test (nums a) (nums b) {:impl :pure}) :exact) 1 0))
      "ks2_g_d" (let [[a b] (string/split "|" argstr)]
                  ((stats/ks-test (nums a) (nums b)
                                  {:impl :pure :alternative :greater}) :statistic))
      "ks2_g_p" (let [[a b] (string/split "|" argstr)]
                  ((stats/ks-test (nums a) (nums b)
                                  {:impl :pure :alternative :greater}) :p-value))
      "ks2_l_d" (let [[a b] (string/split "|" argstr)]
                  ((stats/ks-test (nums a) (nums b)
                                  {:impl :pure :alternative :less}) :statistic))
      "ks2_l_p" (let [[a b] (string/split "|" argstr)]
                  ((stats/ks-test (nums a) (nums b)
                                  {:impl :pure :alternative :less}) :p-value))
      "ks1_d" (let [[z dist p1 p2] (string/split "|" argstr)]
                ((stats/ks-test (nums z) (keyword dist)
                                {:impl :pure
                                 :mean (scan-number p1) :sd (scan-number p2)
                                 :df (scan-number p1)
                                 :shape (scan-number p1) :rate (scan-number p2)})
                 :statistic))
      "ks1_p" (let [[z dist p1 p2] (string/split "|" argstr)]
                ((stats/ks-test (nums z) (keyword dist)
                                {:impl :pure
                                 :mean (scan-number p1) :sd (scan-number p2)
                                 :df (scan-number p1)
                                 :shape (scan-number p1) :rate (scan-number p2)})
                 :p-value))
      "ks1_exact" (let [[z dist p1 p2] (string/split "|" argstr)]
                    (if ((stats/ks-test (nums z) (keyword dist)
                                        {:impl :pure
                                         :mean (scan-number p1) :sd (scan-number p2)
                                         :df (scan-number p1)
                                         :shape (scan-number p1) :rate (scan-number p2)})
                         :exact)
                      1 0))
      "ks1_g_d" (let [[z dist p1 p2] (string/split "|" argstr)]
                  ((stats/ks-test (nums z) (keyword dist)
                                  {:impl :pure :alternative :greater
                                   :mean (scan-number p1) :sd (scan-number p2)
                                   :df (scan-number p1)
                                   :shape (scan-number p1) :rate (scan-number p2)})
                   :statistic))
      "ks1_g_p" (let [[z dist p1 p2] (string/split "|" argstr)]
                  ((stats/ks-test (nums z) (keyword dist)
                                  {:impl :pure :alternative :greater
                                   :mean (scan-number p1) :sd (scan-number p2)
                                   :df (scan-number p1)
                                   :shape (scan-number p1) :rate (scan-number p2)})
                   :p-value))
      "ks1_l_d" (let [[z dist p1 p2] (string/split "|" argstr)]
                  ((stats/ks-test (nums z) (keyword dist)
                                  {:impl :pure :alternative :less
                                   :mean (scan-number p1) :sd (scan-number p2)
                                   :df (scan-number p1)
                                   :shape (scan-number p1) :rate (scan-number p2)})
                   :statistic))
      "ks1_l_p" (let [[z dist p1 p2] (string/split "|" argstr)]
                  ((stats/ks-test (nums z) (keyword dist)
                                  {:impl :pure :alternative :less
                                   :mean (scan-number p1) :sd (scan-number p2)
                                   :df (scan-number p1)
                                   :shape (scan-number p1) :rate (scan-number p2)})
                   :p-value))
      (errorf "unknown fn %s" fn-kw)))
  (def diff6 (math/abs (- got ev)))
  (def rel6 (if (zero? ev) diff6 (/ diff6 (math/abs ev))))
  (assert (or (<= rel6 1e-9) (<= diff6 1e-12))
          (string/format "ecdf-ks row %s: %.17g != %.17g" fn-kw got ev)))

# ---- Kaplan-Meier + log-rank - recorded spot checks ----
(def skm (stats/survfit-km @[1 2 2 3 4 5] @[1 0 1 0 1 0] {}))
(assert (= 6 (skm :n)) "km cohort")
(assert (deep= @[6 5 3 2 1] (skm :n-risk)) "km n-risk (same-time censorings at risk)")
(assert (deep= @[1 1 0 1 0] (skm :n-event)) "km events")
(assert (deep= @[0 1 1 0 1] (skm :n-censor)) "km censorings")
(assert (<= (math/abs (- (get (skm :surv) 1) (/ 2 3))) 1e-15) "km S(2)=2/3")
(assert (<= (math/abs (- (get (skm :surv) 4) (/ 1 3))) 1e-15) "km S stays 1/3 at the tail censor")
(def skm-n (stats/survfit-km @[1 2 2 3 4 5] @[1 0 1 0 1 0] {:impl :num}))
(assert (deep= (skm :surv) (skm-n :surv)) "km vscan :prod == pure loop")
# status 1/2 convention (the coxph precedent): identical curve
(def skm12 (stats/survfit-km @[1 2 2 3 4 5] @[2 1 2 1 2 1] {}))
(assert (deep= (skm :surv) (skm12 :surv)) "status 1/2 == 0/1")
(def slr (stats/survdiff-logrank @[1 2 2 3 4 5 2 3 3 6] @[1 0 1 0 1 0 1 1 0 1]
                                 @[1 1 1 1 1 2 2 2 2 2] {}))
(assert (<= (math/abs (- (slr :statistic) 0.79178082191780785)) 1e-9) "log-rank chisq vs R")
(assert (<= (math/abs (- (slr :p-value) 0.37356218672091324)) 1e-9) "log-rank p vs R")
(assert (= (+ ((slr :observed) 0) ((slr :observed) 1))
            (+ ((slr :expected) 0) ((slr :expected) 1)))
        "log-rank O total == E total (conservation)")
(def slr-n (stats/survdiff-logrank @[1 2 2 3 4 5 2 3 3 6] @[1 0 1 0 1 0 1 1 0 1]
                                   @[1 1 1 1 1 2 2 2 2 2] {:impl :num}))
(assert (= (slr :statistic) (slr-n :statistic)) "log-rank num == pure")

(print "stats smoke ok")
