# project.janet - jpm package manifest for janet-stats.
(declare-project
  :name "janet-stats"
  :description "Statistics library for Janet: distribution functions, hypothesis tests, regression, survival, resampling and Bayesian estimation - pure Janet, R-convention semantics"
  :license "AGPL-3.0"
  :version "0.1.1")

(declare-source
  :source ["jennystats/stats"])
