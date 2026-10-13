"""Simulation checks for the Smith-Waterman processing element."""

from amaranth.sim import Settle, Simulator

from sw_pe import SmithWatermanPE


def evaluate(dut, *, diag, up, left, match):
    """Return the combinational result for one PE input vector."""
    sim = Simulator(dut)
    observed = []

    def process():
        yield dut.diag.eq(diag)
        yield dut.up.eq(up)
        yield dut.left.eq(left)
        yield dut.match.eq(match)
        yield Settle()
        observed.append((yield dut.score))

    sim.add_process(process)
    sim.run()
    return observed[0]


def test_match_uses_diagonal_score():
    assert evaluate(SmithWatermanPE(), diag=8, up=5, left=4, match=1) == 10


def test_gap_from_up_can_win():
    assert evaluate(SmithWatermanPE(), diag=2, up=12, left=3, match=0) == 11


def test_local_alignment_clamps_at_zero():
    assert evaluate(SmithWatermanPE(), diag=0, up=0, left=0, match=0) == 0


def test_score_saturates_at_width_limit():
    pe = SmithWatermanPE(score_width=4)
    assert evaluate(pe, diag=15, up=0, left=0, match=1) == 15
