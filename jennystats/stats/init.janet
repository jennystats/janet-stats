# stats/init.janet - flat re-export of the janet-stats package.
# The F density exports as df-f (the R name df would read as
# degrees-of-freedom elsewhere): (stats/df-f x d1 d2 opts).
(import ./utils :prefix "" :export true)
(import ./dist :prefix "" :export true)
(import ./descriptive :prefix "" :export true)
(import ./ranksum :prefix "" :export true)
(import ./tests :prefix "" :export true)
(import ./regsurv :prefix "" :export true)
(import ./resample :prefix "" :export true)
(import ./bayes :prefix "" :export true)
(import ./multivar :prefix "" :export true)
(import ./linalg :prefix "" :export true)
