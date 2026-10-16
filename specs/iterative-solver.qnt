// iterative-solver.qnt - the termination/error contract shared by the
// error-terminating iterative estimators: glm IRLS (deviance),
// coxph + clogit Newton-Raphson (partial likelihood). Nelder-Mead is
// NOT modeled: it returns a :converged FLAG instead of erroring (a
// different outcome family - see api-contract.md's taxonomy).
//
// The contract (mirrored by the testament suite):
//   1. Results are returned ONLY when converged - never silently
//      non-converged.
//   2. Failure has EXACTLY two causes, both actionable: a non-finite
//      objective (boundary degeneracy) and the iteration bound.
//   3. Convergence is declared only when the criterion held, and only
//      from the family's minimum declaring iteration (IRLS needs the
//      previous objective: 2+; NR: 1+).
//   4. Iterations never exceed the bound; converging exactly AT the
//      bound is success (the code checks convergence before the loop
//      condition re-evaluates).
// Quint verifies the MODEL; the suite executes the invariants against
// the janet code (the impl-resolution.qnt pattern).

module iterativeSolver {

  pure def MAXIT: int = 6

  var phase: str        // "running" | "converged" | "failed"
  var reason: str        // "" | "nonfinite" | "maxit"
  var iters: int        // completed steps
  var minSteps: int     // 1 (NR) | 2 (IRLS: needs the previous objective)
  var lastFinite: bool  // was the last objective finite?
  var lastMet: bool     // did the convergence criterion hold last step?
  var finished: bool    // did the caller receive results?

  action init =
    nondet ms = Set(1, 2).oneOf();
    all {
      phase' = "running",
      reason' = "",
      iters' = 0,
      minSteps' = ms,
      lastFinite' = true,
      lastMet' = false,
      finished' = false,
    }

  // the estimator's step behaviour is NONDETERMINISTIC - the contract
  // must hold for every behaviour an estimator can show

  // objective went non-finite: immediate failure, before any bound
  action stepNonFinite = all {
    phase == "running",
    iters' = iters + 1,
    phase' = "failed",
    reason' = "nonfinite",
    lastFinite' = false,
    lastMet' = false,
    minSteps' = minSteps,
    finished' = false,
  }

  // criterion held and the family may declare: converged (at any
  // iteration up to AND INCLUDING the bound)
  action stepConverged = all {
    phase == "running",
    iters' = iters + 1,
    iters + 1 >= minSteps,
    lastFinite' = true,
    lastMet' = true,
    phase' = "converged",
    reason' = "",
    minSteps' = minSteps,
    finished' = false,
  }

  // criterion held but the family cannot declare yet (IRLS iteration
  // 1), budget remains: continue
  action stepContinueMet = all {
    phase == "running",
    iters' = iters + 1,
    iters + 1 < minSteps,
    iters + 1 < MAXIT,
    lastFinite' = true,
    lastMet' = true,
    phase' = "running",
    reason' = "",
    minSteps' = minSteps,
    finished' = false,
  }

  // criterion missed, budget remains: continue
  action stepContinueMissed = all {
    phase == "running",
    iters' = iters + 1,
    lastFinite' = true,
    lastMet' = false,
    iters + 1 < MAXIT,
    phase' = "running",
    reason' = "",
    minSteps' = minSteps,
    finished' = false,
  }

  // budget exhausted without a declarable criterion: max-iterations
  // failure (criterion met but undeclarable at the bound is included -
  // unreachable at MAXIT >= 2 but part of the family semantics)
  action stepBoundMissed =
    nondet met = Set(true, false).oneOf();
    all {
      phase == "running",
      iters' = iters + 1,
      iters + 1 == MAXIT,
      implies(met, iters + 1 < minSteps),
      lastFinite' = true,
      lastMet' = met,
      phase' = "failed",
      reason' = "maxit",
      minSteps' = minSteps,
      finished' = false,
    }

  // the caller receives results only from a converged run
  action finish = all {
    phase == "converged",
    finished' = true,
    phase' = phase,
    reason' = reason,
    iters' = iters,
    minSteps' = minSteps,
    lastFinite' = lastFinite,
    lastMet' = lastMet,
  }

  action step = any {
    stepNonFinite,
    stepConverged,
    stepContinueMet,
    stepContinueMissed,
    stepBoundMissed,
    finish,
  }

  // ---- the contract, as state invariants --------------------------------

  val boundHeld = iters <= MAXIT

  val phaseOK = Set("running", "converged", "failed").contains(phase)

  // 1: results only when converged
  val resultsOnlyWhenConverged = implies(finished, phase == "converged")

  // 2: failures only from the two defined causes, each with its
  //    signature: non-finite objective, or the bound with no
  //    declarable criterion at it
  val failureCausesExact = implies(
    phase == "failed",
    or {
      and { reason == "nonfinite", not(lastFinite) },
      and {
        reason == "maxit",
        iters == MAXIT,
        or { not(lastMet), iters < minSteps },
      },
    })

  // 3: convergence is earned - criterion held, at the family's
  //    minimum declaring iteration
  val convergedIsEarned = implies(
    phase == "converged",
    and { lastMet, lastFinite, iters >= minSteps, iters >= 1 })

  val reasonDiscipline = implies(phase != "failed", reason == "")

  val inv = and {
    boundHeld,
    phaseOK,
    resultsOnlyWhenConverged,
    failureCausesExact,
    convergedIsEarned,
    reasonDiscipline,
  }
}
