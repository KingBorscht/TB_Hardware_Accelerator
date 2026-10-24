import pytest
from amaranth.sim import Simulator

from sw_pe import (
    TB_DIAG,
    TB_E_EXTEND_BIT,
    TB_F_EXTEND_BIT,
    TB_LEFT,
    TB_STOP,
    TB_UP,
    SmithWatermanPE,
    min_score_width,
)
from sw_pe_testbench import (
    check_score_width,
    run_exhaustive,
    run_matrix_check,
    run_testbench,
)

ZERO = dict(h_diag=0, h_up=0, f_up=0, h_left=0, e_left=0, match=0)


def evaluate(dut, **inputs):
    vector = ZERO | inputs
    sim = Simulator(dut)
    observed = {}

    async def process(ctx):
        for name, value in vector.items():
            ctx.set(getattr(dut, name), value)
        for name in ("h", "e", "f", "tb"):
            observed[name] = ctx.get(getattr(dut, name))

    sim.add_testbench(process)
    sim.run()
    return observed


def test_match_adds_match_score_from_diagonal():
    out = evaluate(SmithWatermanPE(), h_diag=8, match=1)
    assert out["h"] == 10 and out["tb"] & 0b11 == TB_DIAG


def test_mismatch_subtracts_four():
    out = evaluate(SmithWatermanPE(), h_diag=10, match=0)
    assert out["h"] == 6 and out["tb"] & 0b11 == TB_DIAG


def test_local_alignment_clamps_at_zero_and_stops():
    out = evaluate(SmithWatermanPE(), h_diag=1, match=0)
    assert out["h"] == 0 and out["tb"] & 0b11 == TB_STOP


def test_first_gap_base_costs_open_plus_extend():
    out = evaluate(SmithWatermanPE(), h_left=20)
    assert out["e"] == 14 and out["h"] == 14
    assert out["tb"] & 0b11 == TB_LEFT
    assert not (out["tb"] >> TB_E_EXTEND_BIT) & 1


def test_gap_extension_costs_extend_only():
    out = evaluate(SmithWatermanPE(), f_up=14)
    assert out["f"] == 12 and out["h"] == 12
    assert out["tb"] & 0b11 == TB_UP
    assert (out["tb"] >> TB_F_EXTEND_BIT) & 1


def test_extend_tie_chooses_open():
    out = evaluate(SmithWatermanPE(), h_left=20, e_left=16)
    assert out["e"] == 14 and not (out["tb"] >> TB_E_EXTEND_BIT) & 1


def test_diagonal_wins_three_way_tie():
    out = evaluate(SmithWatermanPE(), h_diag=10, h_up=18, h_left=18, match=1)
    assert out["h"] == 12 and out["tb"] & 0b11 == TB_DIAG


def test_score_saturates_at_width_limit():
    assert evaluate(SmithWatermanPE(score_width=4), h_diag=15, match=1)["h"] == 15


def test_rejects_invalid_parameters():
    with pytest.raises(ValueError):
        SmithWatermanPE(score_width=0)
    with pytest.raises(ValueError):
        SmithWatermanPE(gap_extend=-1)
    with pytest.raises(ValueError):
        SmithWatermanPE(score_width=3, gap_open=4, gap_extend=4)


def test_locked_score_width_is_13_bits():
    assert check_score_width() == 13
    assert min_score_width(3072) == 13


def test_directed_and_random_vectors_bit_exact():
    assert run_testbench() > 0


def test_exhaustive_small_width_bit_exact():
    assert run_exhaustive(score_width=3) == 8 ** 5 * 2


def test_matrix_fill_matches_independent_aligner_and_traceback():
    assert run_matrix_check(pair_count=15) > 0