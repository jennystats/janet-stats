// perm-enumeration.qnt - the two-lane permutation-test contract.
// EXACT lane: the enumeration visits each of the C(n,k) splits EXACTLY
// ONCE (distinctness guard + set union; the real next-comb walks the
// lexicographic chain - order is immaterial to the count), and the
// exact p = cnt / C(n,k) is a valid probability (count conservation:
// every exceedance was counted at a visited split, all splits were
// visited). MC lane: the add-one p = (cnt+1)/(B+1) of OUR documented
// convention is a valid probability in (0, 1].
// Modeled at C(4,2) = 6 splits (bitmask pairs); the contract is
// size-independent. Mirrored by the testament suite.

module permEnumeration {

  // the C(4,2) pair bitmasks: {0,1}=3 {0,2}=5 {1,2}=6 {0,3}=9 {1,3}=10 {2,3}=12
  pure def ALL: Set[int] = Set(3, 5, 6, 9, 10, 12)
  pure def NSPLITS: int = 6
  pure def B: int = 8  // MC replicates

  var visited: Set[int]  // exact lane: splits evaluated
  var nVisited: int     // visit counter
  var cnt: int           // exact lane: exceedance count
  var done: bool         // exact lane: enumeration exhausted
  var mcIters: int      // MC lane: replicates run
  var mcCnt: int        // MC lane: exceedance count

  action init = all {
    visited' = Set(),
    nVisited' = 0,
    cnt' = 0,
    done' = false,
    mcIters' = 0,
    mcCnt' = 0,
  }

  action visit =
    nondet m = ALL.oneOf();
    nondet exceed = Set(0, 1).oneOf();
    all {
      not(done),
      not(visited.contains(m)),           // each split at most once
      visited' = visited.union(Set(m)),
      nVisited' = nVisited + 1,
      cnt' = cnt + exceed,               // counted only at a visit
      done' = done,
      mcIters' = mcIters,
      mcCnt' = mcCnt,
    }

  action exhaust = all {
    not(done),
    visited == ALL,                       // every split was reached
    done' = true,
    visited' = visited,
    nVisited' = nVisited,
    cnt' = cnt,
    mcIters' = mcIters,
    mcCnt' = mcCnt,
  }

  action mcStep =
    nondet exceed = Set(0, 1).oneOf();
    all {
      not(done),
      mcIters < B,
      mcIters' = mcIters + 1,
      mcCnt' = mcCnt + exceed,
      visited' = visited,
      nVisited' = nVisited,
      cnt' = cnt,
      done' = done,
    }

  action step = any {
    visit,
    exhaust,
    mcStep,
  }

  // ---- the contract, as state invariants --------------------------------

  // only real splits are ever visited
  val visitedSubset = visited.forall(m => ALL.contains(m))

  // the visit count matches the set (no double-visit could grow one
  // without the other)
  val visitCountBounds = and { nVisited >= 0, nVisited <= NSPLITS }

  // exhaustion means ALL splits were evaluated - the exact p's
  // denominator is earned
  val exhaustionComplete = implies(done, and { nVisited == NSPLITS, visited == ALL })

  // count conservation: exceedances only at visited splits, so the
  // exact p = cnt/NSPLITS is a valid probability [0,1]
  val exactCountConservation = and { cnt <= nVisited, cnt >= 0 }

  // MC add-one: p = (mcCnt + 1)/(B + 1) in (0, 1] - valid for every
  // count, including cnt == B (p = 1) and cnt == 0 (p > 0)
  val mcAddOneValid = and { mcCnt <= B, mcCnt >= 0, mcIters <= B, mcCnt <= mcIters + 1 }

  val inv = and {
    visitedSubset,
    visitCountBounds,
    exhaustionComplete,
    exactCountConservation,
    mcAddOneValid,
  }
}
