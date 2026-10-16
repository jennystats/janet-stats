// km-logrank.qnt - FORWARD CONTRACT for Kaplan-Meier + log-rank
// (oracles survival::survfit/survdiff).
// Written in the SPECIFICATION PHASE, before implementation (the
// corrected sequencing); the implementation must satisfy these
// invariants, and the pre-written suite mirrors
// turn green at landing.
//
// THE CONVENTION THIS SPEC PINS (R survfit, decode-from-source at
// slice time): at a tied time point, DEATHS are processed FIRST -
// the KM factor (Y - d)/Y uses the risk set Y INCLUDING those
// censored at the same time; censorings leave the risk set only
// AFTER the factor. Getting this order wrong is THE classic KM bug.
//
// Log-rank (14b): the same event times carry the two-group
// decomposition Y = Y1 + Y2, and observed deaths attribute exactly:
// O1 + O2 == d at every event time (the O-E statistic is built on
// it).

module kmLogrank {

  pure def N: int = 6   // cohort size (model-checking scale)

  var atrisk: int       // total at risk (BEFORE the current time)
  var y1: int           // group-1 at risk (y2 = atrisk - y1)
  var sNum: int         // KM S numerator
  var sDen: int         // KM S denominator
  var deaths: int       // total observed deaths
  var censored: int     // total observed censorings
  var o1: int           // group-1 observed deaths
  var o2: int           // group-2 observed deaths
  var done: bool

  action init =
    nondet y = 0.to(N).oneOf();
    all {
      atrisk' = N,
      y1' = y,
      sNum' = 1,
      sDen' = 1,
      deaths' = 0,
      censored' = 0,
      o1' = 0,
      o2' = 0,
      done' = false,
    }

  // one distinct time point: d1/d2 deaths (groups), c1/c2 censorings.
  // The KM factor is applied FIRST (risk set = atrisk, INCLUDING the
  // same-time censorings), THEN censorings leave the set.
  action eventTime =
    nondet d1 = 0.to(N).oneOf();
    nondet d2 = 0.to(N).oneOf();
    nondet c1 = 0.to(N).oneOf();
    nondet c2 = 0.to(N).oneOf();
    all {
      not(done),
      d1 + d2 + c1 + c2 >= 1,        // a time point processes someone
      d1 <= y1,
      d2 <= atrisk - y1,
      c1 <= y1 - d1,
      c2 <= atrisk - y1 - d2,
      // KM factor: deaths before censorings - the risk set for the
      // factor is the PRE-censoring atrisk
      sNum' = sNum * (atrisk - d1 - d2),
      sDen' = sDen * atrisk,
      // log-rank attribution conserves the deaths at this time
      o1' = o1 + d1,
      o2' = o2 + d2,
      deaths' = deaths + d1 + d2,
      censored' = censored + c1 + c2,
      // only now do censorings (and deaths) leave the risk sets
      atrisk' = atrisk - d1 - d2 - c1 - c2,
      y1' = y1 - d1 - c1,
      done' = false,
    }

  action finish = all {
    not(done),
    atrisk == 0,
    done' = true,
    atrisk' = atrisk,
    y1' = y1,
    sNum' = sNum,
    sDen' = sDen,
    deaths' = deaths,
    censored' = censored,
    o1' = o1,
    o2' = o2,
  }

  // the exhausted state stutters (quint requires no deadlocks)
  action doneIdle = all {
    done,
    done' = true,
    atrisk' = atrisk,
    y1' = y1,
    sNum' = sNum,
    sDen' = sDen,
    deaths' = deaths,
    censored' = censored,
    o1' = o1,
    o2' = o2,
  }

  action step = any {
    eventTime,
    finish,
    doneIdle,
  }

  // ---- the contract, as state invariants --------------------------------

  // S is a survival function: bounded by [0, 1], non-increasing
  // (each factor (Y-d)/Y <= 1), and the denominator stays positive
  val sBounded = and { sNum <= sDen, sNum >= 0, sDen >= 1 }

  // risk sets stay non-negative and the two-group decomposition holds
  val riskSetsOK = and { atrisk >= 0, y1 >= 0, y1 <= atrisk }

  // COHORT CONSERVATION: everyone is accounted for exactly once -
  // still at risk, dead, or censored
  val cohortConservation = deaths + censored + atrisk == N

  // deaths at risk times only, attributed exactly across the groups
  val deathAttribution = and {
    o1 + o2 == deaths,
    o1 <= deaths,
    o2 <= deaths,
    deaths <= N,
    censored <= N,
  }

  // S(0) = 1 at the start: no factor applied yet means sNum == sDen
  // until the first event (captured by: before any death, S == 1)
  val sStartsAtOne = implies(deaths == 0, sNum == sDen)

  val inv = and {
    sBounded,
    riskSetsOK,
    cohortConservation,
    deathAttribution,
    sStartsAtOne,
  }
}
