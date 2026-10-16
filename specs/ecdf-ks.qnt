// ecdf-ks.qnt - FORWARD CONTRACT for the ECDF + one/two-
// sample KS; oracle base R ecdf/ks.test). Specification-phase artifact
// (pre-implementation); the pre-written suite mirrors
// turn green at landing.
//
// The ECDF construction contract: F(x) = #{x_i <= x} / n - the fold
// over sorted data consumes TIE BLOCKS whole (every tied observation
// shares ONE F value - R ecdf's right-continuous convention), F is
// monotone and bounded by [0,1], and F(max) = 1 at exhaustion.
// Modeled on data [1, 1, 2, 3, 3, 3]: blocks of sizes 2, 1, 3.
//
// The KS statistic contract (two-sample, stated on a concrete pair):
// D = sup over the UNION grid of both samples, each sample's ECDF
// evaluated with its OWN denominator; D is bounded [0,1] and
// symmetric in its arguments. Concrete pair A=[1,1,2] B=[1,3,3]
// (n=3 each, values scaled by 3): grid {1,2,3}; F_A = 2/3, 1, 1;
// F_B = 1/3, 1/3, 1; D = 2/3 (at x=2).

module ecdfKS {

  pure def N: int = 6

  var blocks: int   // tie blocks consumed (0..3)
  var f: int         // ECDF numerator: count of observations <= last value
  var done: bool

  action init = all {
    blocks' = 0,
    f' = 0,
    done' = false,
  }

  // the fold: a tie block is consumed WHOLE (all tied observations
  // share one F value)
  action consumeBlock = all {
    not(done),
    blocks < 3,
    any {
      all { blocks == 0, f' = f + 2 },
      all { blocks == 1, f' = f + 1 },
      all { blocks == 2, f' = f + 3 },
    },
    blocks' = blocks + 1,
    done' = false,
  }

  action finish = all {
    not(done),
    blocks == 3,
    done' = true,
    blocks' = blocks,
    f' = f,
  }

  action doneIdle = all {
    done,
    done' = true,
    blocks' = blocks,
    f' = f,
  }

  action step = any {
    consumeBlock,
    finish,
    doneIdle,
  }

  // ---- the ECDF contract, as state invariants ---------------------------

  // F in [0, 1]
  val fBounded = and { f >= 0, f <= N }

  // THE TIE CONVENTION: f is exactly the cumulative block count -
  // ties were consumed whole (a partial-block f would land between
  // these values)
  val tieBlocksWhole = or {
    and { blocks == 0, f == 0 },
    and { blocks == 1, f == 2 },
    and { blocks == 2, f == 3 },
    and { blocks == 3, f == 6 },
  }

  // F(max) = 1 at exhaustion
  val fMaxIsOne = implies(done, f == N)

  // ---- the KS statistic contract (concrete pair, scaled x3) -----------
  // F_A(1)=2 F_A(2)=3 F_A(3)=3; F_B(1)=1 F_B(2)=1 F_B(3)=3
  pure def FA1: int = 2
  pure def FA2: int = 3
  pure def FA3: int = 3
  pure def FB1: int = 1
  pure def FB2: int = 1
  pure def FB3: int = 3

  // D over the UNION grid {1, 2, 3}: max(|FA-FB|) = max(1, 2, 0) = 2
  // (scaled; D = 2/3)
  pure def DAB: int = 2
  pure def DBA: int = 2

  // D bounded [0, 1] (scaled [0, 3]) and SYMMETRIC in its arguments
  val dBounded = and { DAB >= 0, DAB <= 3 }
  val dSymmetric = DAB == DBA

  val inv = and {
    fBounded,
    tieBlocksWhole,
    fMaxIsOne,
    dBounded,
    dSymmetric,
  }
}
