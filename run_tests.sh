#!/bin/sh
# run_tests.sh - janet-stats smoke + conformance suite. The repo root
# comes FIRST so the tests exercise this tree's stats modules; typed and
# kit resolve from the janet-kit dependency on $JANET_PATH. The
# $HOME/.local/lib/janet seed serves that dependency when the env has none.
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
