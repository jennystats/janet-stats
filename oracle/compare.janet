# comparator: post-vs-R and pre-vs-R rel errors, keyed join.
(import typed :prefix "")
(import kit :prefix "")

(defn/typed load-cells {:args [:string] :ret :table} [path]
  (def t @{})
  (each line (string/split "\n" (string (slurp path)))
    (when (not (empty? line))
      (def parts (string/split "\t" line))
      (put t (get parts 0) (string (get parts 1)))))
  t)

(def r-t (load-cells "r.txt"))
(def post-t (load-cells "post.txt"))
(def pre-t (load-cells "pre.txt"))

(defn/typed parse1 {:args [:string] :ret :number} [s]
  (if (string/find "THROW" s) math/nan
    (try (let [v (parse s)] (if (number? v) v math/nan)) ([_] math/nan))))

(defn/typed relerr {:args [:number :number] :ret :number} [v r0]
  (def d (math/abs (- v r0)))
  (cond (= r0 0) (if (= v 0) 0 math/inf)
    (= v math/nan) math/nan
    (/ d (math/abs r0))))

(each k (sort (keys r-t))
  (def rv (get r-t k))
  (def pv (get post-t k "MISSING"))
  (def qv (get pre-t k "MISSING"))
  (def r0 (parse1 rv))
  (def v (parse1 pv))
  (def w (if (= qv "MISSING") math/nan (parse1 qv)))
  (printf "%-30s post=%-25s pre=%-25s R=%-25s rel-post=%s rel-pre=%s"
          k pv qv rv
          (if (or (= v math/nan) (= r0 math/nan)) "n/a" (string/format "%.3g" (relerr v r0)))
          (if (or (= w math/nan) (= r0 math/nan)) "n/a" (string/format "%.3g" (relerr w r0)))))
