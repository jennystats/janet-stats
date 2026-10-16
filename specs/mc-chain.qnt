// mc-chain.qnt - the MH/Gibbs bookkeeping contract: draws are
// appended ONLY post-burn, exactly one per iteration (accept or
// reject alike); the Robbins-Monro log-sd adaptation applies ONLY
// during burn-in (post-burn adaptation biases the chain - the
// classic silent bug); accept accounting splits at the boundary and
// conserves the total. Mirrored by the testament suite (the
// impl-resolution.qnt pattern).

module mcChain {

  pure def STEPS: int = 10

  var i: int            // iterations completed
  var burn: int         // burn-in length (fixed per chain)
  var adapted: int      // Robbins-Monro adaptation runs (burn-in ONLY)
  var draws: int        // draws appended (post-burn only)
  var acceptTotal: int // all accepted proposals
  var acceptBurn: int  // accepted during burn-in
  var acceptPost: int  // accepted post-burn

  action init =
    nondet b = 0.to(STEPS).oneOf();
    all {
      i' = 0,
      burn' = b,
      adapted' = 0,
      draws' = 0,
      acceptTotal' = 0,
      acceptBurn' = 0,
      acceptPost' = 0,
    }

  // burn-in iteration, proposal ACCEPTED (adaptation runs on EVERY
  // burn-in iteration - the default adaptive=true chain)
  action acceptBurnStep = all {
    i < burn,
    i' = i + 1,
    burn' = burn,
    adapted' = adapted + 1,
    draws' = draws,
    acceptTotal' = acceptTotal + 1,
    acceptBurn' = acceptBurn + 1,
    acceptPost' = acceptPost,
  }

  action rejectBurnStep = all {
    i < burn,
    i' = i + 1,
    burn' = burn,
    adapted' = adapted + 1,
    draws' = draws,
    acceptTotal' = acceptTotal,
    acceptBurn' = acceptBurn,
    acceptPost' = acceptPost,
  }

  // post-burn iteration: a draw is appended on EVERY iteration -
  // acceptance only touches the accounting
  action acceptPostStep = all {
    or { i >= burn, burn == 0 },
    i < STEPS,
    i' = i + 1,
    burn' = burn,
    adapted' = adapted,
    draws' = draws + 1,
    acceptTotal' = acceptTotal + 1,
    acceptBurn' = acceptBurn,
    acceptPost' = acceptPost + 1,
  }

  action rejectPostStep = all {
    or { i >= burn, burn == 0 },
    i < STEPS,
    i' = i + 1,
    burn' = burn,
    adapted' = adapted,
    draws' = draws + 1,
    acceptTotal' = acceptTotal,
    acceptBurn' = acceptBurn,
    acceptPost' = acceptPost,
  }

  action step = any {
    acceptBurnStep,
    rejectBurnStep,
    acceptPostStep,
    rejectPostStep,
  }

  // ---- the contract, as state invariants --------------------------------

  val iBound = i <= STEPS

  // draws conservation: exactly one per post-burn iteration, none
  // before the boundary
  val drawsConservation = or {
    and { i <= burn, draws == 0 },
    and { i > burn, draws == i - burn },
  }

  // acceptance accounting conserves the total and the split boundary
  val acceptSplit = acceptBurn + acceptPost == acceptTotal
  val acceptBurnBound = and { acceptBurn <= burn, acceptBurn <= i }

  // THE adaptation boundary: Robbins-Monro ran exactly during the
  // burn-in iterations and NEVER after (post-burn adaptation biases
  // the chain - the classic silent bug this pins)
  val adaptationOnlyInBurn = or {
    and { i <= burn, adapted == i },
    and { i > burn, adapted == burn },
  }

  val inv = and {
    iBound,
    drawsConservation,
    acceptSplit,
    acceptBurnBound,
    adaptationOnlyInBurn,
  }
}
