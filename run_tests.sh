#!/bin/sh
# run_tests.sh - janet-stats smoke + conformance suite. Repo modules
# come FIRST so the tests always exercise this tree (typed resolves
# in-repo), not an installed copy. The $HOME/.local/lib/janet seed
# serves the OPTIONAL janet-num companion when the env has none.
DIR="$(cd "$(dirname "$0")" && pwd)"
[ -n "$JANET_PATH" ] || export JANET_PATH="$HOME/.local/lib/janet"
export JANET_PATH="$DIR:$DIR/tests/lib:$JANET_PATH"
fail=0
if janet "$DIR/tests/smoke-stats.janet"; then
  echo "PASS smoke-stats"
else
  echo "FAIL smoke-stats"
  fail=1
fi
if janet "$DIR/tests/suite-stats.janet"; then
  echo "PASS suite-stats"
else
  echo "FAIL suite-stats"
  fail=1
fi
if [ "$fail" -eq 0 ]; then
  echo "RESULT: pass"
else
  echo "RESULT: fail"
  exit 1
fi
