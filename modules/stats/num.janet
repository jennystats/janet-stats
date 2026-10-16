# stats/num.janet - the :impl RESOLUTION CONTRACT, single owner.
# Formal spec: the janet-num repository's specs/impl-resolution.qnt
# (quint-verified; this module must satisfy all six invariants).
# DEFAULT = :impl - :num when libjn is
# present, :pure when absent (the ONLY dynamic path); explicit :num is
# STRICT - absent libjn throws actionable, never a silent :pure;
# explicit :pure never loads the library; the probe runs ONCE per
# process. Since 1.0 the default is the dynamic :impl resolution.
# STATE TABLE: module vars are
# snapshotted by import - a table is the live per-process state the
# tests manipulate (num/state).
(import typed :prefix "")

(def state
  @{:presence :unknown    # :unknown | :present | :absent
    :default-impl :impl})   # the 1.0 default - dynamic

(defn/typed probe {:args [] :ret :keyword} []
  "THE one loader invocation per process (quint: probeAtMostOnce);
  janet does not cache failed requires - :presence is the cache.
  The gate lives here too (a direct second (num/probe) call must not
  re-run the loader and re-warn)."
  (def cached (state :presence))
  (when (not= cached :unknown) cached)
  # module envs inherit ROOT dyns, not the
  # calling script's setdyn syspath self-heal - a plain
  # `janet script.janet` could not see ~/.local/lib/janet here and
  # silently degraded to :pure on a libjn machine. Heal THIS env's
  # syspath (janet-num.janet 1.0.1+ heals its own env for its typed
  # import the same way).
  # when presence is UNVERIFIABLE (HOME
  # unset/empty and the loader failed), the silent :pure resolution
  # gets a LOUD one-per-process warning - the user knows what ran.
  # the heal predicate is ALIGNED with janet-num's own
  # heal (a deliberate copy, not a shared module - the vendored
  # refresh stays byte-simple; the formal contract is
  # the janet-num repository's specs/syspath-heal.qnt): no heal when HOME is nil/empty,
  # no clobber of a non-default inherited syspath, first EXISTING dir
  # of the short list, single dir.
  (def home (get (os/environ) "HOME"))
  (def inherited (string (or (dyn *syspath*) "")))
  (when (and home (not (empty? home))
              (or (empty? inherited)
                  (= inherited "/usr/local/lib/janet")))
    (def candidates [(string home "/.local/lib/janet")
                     "/usr/local/lib/janet"
                     "/usr/lib/janet"])
    (def first-dir
      (find (fn [p]
              # a regular file at a candidate path must not win the
              # scan - same predicate as janet-num's heal
              (def s (try (os/stat p) ([_] nil)))
              (and s (= (s :mode) :directory)))
            candidates))
    (when first-dir (setdyn *syspath* first-dir)))
  (var load-err nil)
  (def r (if (try (do (require "janet-num") true)
               ([e] (set load-err (string e)) false))
          :present :absent))
  # the require failure is CAPTURED,
  # never swallowed - the only silent catch in the stats tree is
  # gone. When janet-num was FOUND but FAILED (its A1 error text
  # distinguishes "cannot load"/"cannot find typed" from the
  # not-found class), its own actionable error is printed verbatim -
  # a wrong remedy is worse than none. The on-disk check requires
  # :file mode (a DIRECTORY in $JANET_NUM_LIB no longer satisfies).
  (when (= r :absent)
    (def pinned (get (os/environ) "JANET_NUM_LIB"))
    (defn stat-file [p]
      (def s (try (os/stat p) ([_] nil)))
      (and p s (= (s :mode) :file)))
    (def on-disk
      (find stat-file
            [(when pinned pinned)
             (string (or home "") "/.local/lib/janet-num/libjn.so")]))
    (cond
      # found-but-FAILED: any load error that is NOT the not-found
      # class (janet-num's own remedy text when it has one, the raw
      # dlopen error otherwise - surfacing verbatim is still honest,
      # the false
      # "not findable" was the defect)
      (and load-err (not (string/find "could not find module" load-err)))
      (eprintf "WARNING stats: the janet-num module was found but FAILED to load:\n%s\nThis run is pure janet (the R-verified reference lane)." load-err)
      on-disk
      (eprintf "WARNING stats: libjn is on disk but the janet-num module was not\nfindable - check $JANET_PATH/$JANET_NUM_LIB include the module dir,\nor run libs/install-janet-num.sh. This run is pure janet\n(the R-verified reference lane).")))
  (put state :presence r)
  r)

(defn/typed present? {:args [] :ret :boolean} []
  "True when the janet-num library is loadable (probes once per
  process, then serves the cached answer)."
  (when (= (state :presence) :unknown) (probe))
  (= (state :presence) :present))

(defn/typed resolve-impl {:args [:any] :ret :keyword} [opts]
  "Resolve the opts table's :impl to :pure | :num. Explicit :num on
  an absent library throws actionable (strict - never a silent
  :pure fallback); :pure never probes (quint: pureNeverProbes);
  :impl never errors (implNeverErrors)."
  (case (get opts :impl (state :default-impl))
    :pure :pure
    :num (if (present?) :num
           (error "stats: {:impl :num} requires libjn - run libs/install-janet-num.sh (needs gcc) or set $JANET_NUM_LIB"))
    :impl (if (present?) :num :pure)
    (error "stats: :impl must be :pure, :num, or :impl")))

(defn/typed kernels {:args [] :ret :any} []
  "The kernel accessor (name -> unwrapped value); defn/typed module
  bindings are ref-tables {:value ...} - unwrap on access. Only
  valid when present (resolve-impl gates every call site)."
  (def m (require "janet-num"))
  (fn [nm]
    (def w (get m nm))
    (if (table? w) (get w :value) w)))

(defn/typed stage {:args [:any :any] :ret :any} [x k]
  "Transparent input dispatch: a janet array stages into a packed
  handle (the one boundary copy); a handle passes through - pack
  once, chain kernel calls, unpack at the pipeline's end. k = the
  kernels accessor."
  (if (indexed? x) ((k 'pack-f64) x) x))

(defn/typed count {:args [:any] :ret :number} [x]
  "The resident-input length: janet array or packed handle (:n)."
  (if (indexed? x) (length x) (x :n)))

(defn check-impl-landed
  "Guard for fns WITH a landed num path: resolves (validates,
  strict-throws, probes once); with x (a resident input), :pure plus
  a packed handle throws actionable (handles are num-resident - a
  pure loop must never mis-iterate one). Scalar calls run the
  identical pure formulas (honest scalar policy). Plain defn:
  typed cannot spell &opt."
  [opts &opt x]
  (def impl (resolve-impl opts))
  (when (and x (= :pure impl) (struct? x))
    (error "stats: :pure path consumes scalars/arrays; a packed handle needs {:impl :num} or {:impl :impl}"))
  impl)
