# Spec - the controlled-CSV fixture dialect (stats fixtures, v1)

Non-Quint specification (interface contract). Consumers: the test
harnesses (the conformance suite + the smoke). The two shipped parser
bugs (middle-field truncation; `string/slice` on short strings) were
both grammar-ambiguity failures - this spec is the fix for that bug
class. Source of truth: the committed fixtures under
`tests/fixtures/stats` (R-generated, base R only, fixed seeds).

## 1. Dialect classes

Two file classes only:

- **RESULTS files** (`tier-<t>-tests.csv`): records of
  `fn,args,expected` (tier-D+ variant: `fn,dataset,expected`).
- **DATA files** (`tier-<t>-<dataset>.csv`): bare numeric matrices,
  row-major, one line per row, `%.17g` doubles.

## 2. Grammar (results files)

```
file     := header LF record*
header   := "fn" "," ("args" | "dataset") "," "expected" LF
record   := fn "," args "," expected LF
fn       := ident            ; MUST NOT contain ',' or '"' or LF
args     := text             ; MAY contain ','; MUST NOT contain '"' or LF
expected := double          ; %.17g, parseable by scan-number
```

**THE JOIN RULE (normative)**: a record is EXACTLY 3 logical fields,
but `args` may contain commas - the naive `fields[1]` read TRUNCATES
at the first comma (the twice-caught bug). Every parser MUST
reconstruct:

- `fn       = fields[0]`
- `args     = join(fields[1..n-2], ",")`
- `expected = fields[n-1]`  (the LAST field, always)

No quoting, no escaping, no NA, no embedded newlines, no comment
lines - THAT is what makes the dialect controlled (the R
`read.csv`-guessing drift class does not apply to it).

## 3. Grammar (data files)

```
file  := (row LF)+
row   := double ("," double)*     ; %.17g, row-major matrix row
```

- NO header is required; when a generator emits one (R `write.csv`
  quotes headers, e.g. tier-D), the header line CONTAINS `"` -
  consumers skip any line containing a double-quote (the shipped
  quote-skip rule). Headerless generators (tier-G `write.table`,
  `col.names=FALSE`) emit bare rows.
- Column count is constant within a file; the consumer treats it as
  the matrix width.

## 4. Number format

- Producers emit `sprintf("%.17g", x)` - round-trip exact for doubles.
- Scientific notation (`1.9e-06`, `-2.5E-3`) is valid; consumers parse
  with `scan-number`.
- Consumers MUST NOT re-format expected values (comparison happens in
  parsed double space, tolerance-based).

## 5. Parser contract (every consumer: harnesses, smoke, future tiers)

1. **Skip rules**: lines that are empty; lines containing `fn,` (the
   header); lines containing `"` (data-file headers).
2. **Join rule** (§2) for results records.
3. **Display-slice guard**: truncating an args string for messages
   MUST clamp - `(string/slice s 0 (min K (length s)))`; a bare
   `(string/slice s 0 K)` CRASHES on short strings (the shipped bug).
4. **Trailing-line tolerance**: `string/split` yields a trailing `""`
   - filter empties.
5. **Dataset reference convention**: the `dataset` field names the
   sibling file `tier-<tier>-<dataset>.csv` in the same directory;
   load via the data-file grammar.

## 6. Producer contract (R side)

- Base R ONLY (`write.csv`, `writeLines`, `sprintf`) - zero CRAN,
  even as tooling (the fixture-lane resolution).
- Fixed seeds; NO trailing `cat()` to stdout (the generator's stdout
  IS the results file via redirect - a stray print corrupts it; the
  shipped incident).
- `%.17g` everywhere; LF endings; ASCII.

## 7. Dialect invariants (the conformance bar)

I1. No field contains `"`, LF, or an NA token.
I2. `expected` is always the last field of its record.
I3. `fn` never contains `,`.
I4. Every line is total: header, blank, or a complete record/row.
I5. Round-trip: parse(format(x)) == x for every double emitted.

A parser or producer change that breaks any invariant is a dialect
change - it requires this spec's revision, not a silent edit.
