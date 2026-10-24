from itertools import product
from random import Random

from amaranth.sim import Simulator

from sw_pe import (
    DEFAULT_GAP_EXTEND,
    DEFAULT_GAP_OPEN,
    DEFAULT_L_MAX,
    DEFAULT_MATCH_SCORE,
    DEFAULT_MISMATCH_PENALTY,
    DEFAULT_SCORE_WIDTH,
    TB_DIAG,
    TB_E_EXTEND_BIT,
    TB_F_EXTEND_BIT,
    TB_LEFT,
    TB_STOP,
    TB_UP,
    SmithWatermanPE,
    min_score_width,
)

INPUTS = ("h_diag", "h_up", "f_up", "h_left", "e_left", "match")
OUTPUTS = ("h", "e", "f", "tb")


def reference_cell(*, h_diag, h_up, f_up, h_left, e_left, match,
                   match_score=DEFAULT_MATCH_SCORE,
                   mismatch_penalty=DEFAULT_MISMATCH_PENALTY,
                   gap_open=DEFAULT_GAP_OPEN,
                   gap_extend=DEFAULT_GAP_EXTEND,
                   score_width=DEFAULT_SCORE_WIDTH):
    limit = (1 << score_width) - 1
    first = gap_open + gap_extend

    if match:
        diag = min(limit, h_diag + match_score)
    else:
        diag = max(0, h_diag - mismatch_penalty)

    e_open, e_ext = max(0, h_left - first), max(0, e_left - gap_extend)
    f_open, f_ext = max(0, h_up - first), max(0, f_up - gap_extend)
    e_extended, f_extended = e_ext > e_open, f_ext > f_open
    e = e_ext if e_extended else e_open
    f = f_ext if f_extended else f_open

    h = max(diag, e, f)
    if h == 0:
        source = TB_STOP
    elif diag == h:
        source = TB_DIAG
    elif e == h:
        source = TB_LEFT
    else:
        source = TB_UP

    tb = source | (int(e_extended) << TB_E_EXTEND_BIT) | (int(f_extended) << TB_F_EXTEND_BIT)
    return {"h": h, "e": e, "f": f, "tb": tb}


def gotoh_best_score(query, ref, *, match_score=DEFAULT_MATCH_SCORE,
                     mismatch_penalty=DEFAULT_MISMATCH_PENALTY,
                     gap_open=DEFAULT_GAP_OPEN, gap_extend=DEFAULT_GAP_EXTEND):
    neg = -(10 ** 9)
    rows, cols = len(query) + 1, len(ref) + 1
    h = [[0] * cols for _ in range(rows)]
    e = [[neg] * cols for _ in range(rows)]
    f = [[neg] * cols for _ in range(rows)]
    best = 0
    for i in range(1, rows):
        for j in range(1, cols):
            e[i][j] = max(h[i][j - 1] - gap_open - gap_extend, e[i][j - 1] - gap_extend)
            f[i][j] = max(h[i - 1][j] - gap_open - gap_extend, f[i - 1][j] - gap_extend)
            sub = match_score if query[i - 1] == ref[j - 1] else -mismatch_penalty
            h[i][j] = max(0, h[i - 1][j - 1] + sub, e[i][j], f[i][j])
            best = max(best, h[i][j])
    return best


def _pe_kwargs(pe):
    return dict(match_score=pe.match_score, mismatch_penalty=pe.mismatch_penalty,
                gap_open=pe.gap_open, gap_extend=pe.gap_extend,
                score_width=pe.score_width)


def _directed_vectors(pe):
    limit = (1 << pe.score_width) - 1
    first = pe.gap_open + pe.gap_extend
    ext = pe.gap_extend
    vectors = []
    for match in (0, 1):
        for value in (0, 1, pe.mismatch_penalty - 1, pe.mismatch_penalty,
                      limit - pe.match_score, limit - 1, limit):
            vectors.append(dict(ZERO, h_diag=value, match=match))
        for value in (first - 1, first, first + 1):
            vectors.append(dict(ZERO, h_left=value, match=match))
            vectors.append(dict(ZERO, h_up=value, match=match))
        for base in (first, first + 3, 40):
            vectors.append(dict(ZERO, h_up=base, f_up=base - first + ext,
                                h_left=base, e_left=base - first + ext, match=match))
        vectors.append(dict(ZERO, h_diag=limit, h_up=limit, f_up=limit,
                            h_left=limit, e_left=limit, match=match))
        vectors.append(dict(ZERO, match=match))
    vectors.append(dict(ZERO, h_diag=10, h_up=10 + first, h_left=10 + first, match=1))
    vectors.append(dict(ZERO, h_up=20 + first, h_left=20 + first, match=0))
    return vectors


ZERO = dict(h_diag=0, h_up=0, f_up=0, h_left=0, e_left=0, match=0)


def _random_vectors(pe, count, rng):
    full = 1 << pe.score_width
    small = min(2 * (pe.gap_open + pe.gap_extend) + 2, full)
    vectors = []
    for index in range(count):
        span = full if index % 2 == 0 else small
        vector = {name: rng.randrange(span) for name in INPUTS[:-1]}
        vector["match"] = rng.randrange(2)
        vectors.append(vector)
    return vectors


def _drive_and_check(pe, vectors, label):
    sim = Simulator(pe)

    async def process(ctx):
        for index, vector in enumerate(vectors):
            for name in INPUTS:
                ctx.set(getattr(pe, name), vector[name])
            actual = {name: ctx.get(getattr(pe, name)) for name in OUTPUTS}
            expected = reference_cell(**vector, **_pe_kwargs(pe))
            assert actual == expected, (
                f"{label} vector {index}: expected {expected}, got {actual}; "
                f"inputs={vector}"
            )

    sim.add_testbench(process)
    sim.run()
    return len(vectors)


def run_testbench(vector_count=2000, seed=0x5A17):
    pe = SmithWatermanPE()
    rng = Random(seed)
    vectors = _directed_vectors(pe) + _random_vectors(pe, vector_count, rng)
    return _drive_and_check(pe, vectors, "scope-v2")


def run_exhaustive(score_width=3):
    pe = SmithWatermanPE(score_width=score_width)
    values = range(1 << score_width)
    vectors = [
        dict(h_diag=a, h_up=b, f_up=c, h_left=d, e_left=e, match=m)
        for a, b, c, d, e, m in product(values, values, values, values, values, (0, 1))
    ]
    return _drive_and_check(pe, vectors, f"exhaustive-w{score_width}")


def _mutate(rng, ref, sub_rate=0.08, indel_rate=0.06):
    out = []
    for base in ref:
        roll = rng.random()
        if roll < indel_rate / 2:
            continue
        if roll < indel_rate:
            out.append(base)
            out.extend(rng.choice("ACGT") for _ in range(rng.randrange(1, 4)))
        elif roll < indel_rate + sub_rate:
            out.append(rng.choice([b for b in "ACGT" if b != base]))
        else:
            out.append(base)
    return "".join(out)


def _walk_traceback(query, ref, tb, start, pe):
    first = pe.gap_open + pe.gap_extend
    i, j = start
    state = "H"
    score = 0
    while True:
        if state == "H":
            if i == 0 or j == 0:
                break
            source = tb[i][j] & 0b11
            if source == TB_STOP:
                break
            if source == TB_DIAG:
                same = query[i - 1] == ref[j - 1]
                score += pe.match_score if same else -pe.mismatch_penalty
                i, j = i - 1, j - 1
            elif source == TB_LEFT:
                state = "E"
            else:
                state = "F"
        elif state == "E":
            extended = (tb[i][j] >> TB_E_EXTEND_BIT) & 1
            score -= pe.gap_extend if extended else first
            j -= 1
            state = "E" if extended else "H"
        else:
            extended = (tb[i][j] >> TB_F_EXTEND_BIT) & 1
            score -= pe.gap_extend if extended else first
            i -= 1
            state = "F" if extended else "H"
        assert i >= 0 and j >= 0, "traceback walked off the matrix"
    return score


def run_matrix_check(pair_count=40, ref_length=36, seed=0xA11C):
    pe = SmithWatermanPE()
    rng = Random(seed)
    pairs = []
    for _ in range(pair_count):
        ref = "".join(rng.choice("ACGT") for _ in range(ref_length))
        pairs.append((_mutate(rng, ref), ref))
    pairs.append(("A" * 12, "C" * 12))

    sim = Simulator(pe)

    async def process(ctx):
        for pair_index, (query, ref) in enumerate(pairs):
            rows, cols = len(query) + 1, len(ref) + 1
            h = [[0] * cols for _ in range(rows)]
            e = [[0] * cols for _ in range(rows)]
            f = [[0] * cols for _ in range(rows)]
            tb = [[0] * cols for _ in range(rows)]
            best, best_cell = 0, (0, 0)
            for i in range(1, rows):
                for j in range(1, cols):
                    ctx.set(pe.h_diag, h[i - 1][j - 1])
                    ctx.set(pe.h_up, h[i - 1][j])
                    ctx.set(pe.f_up, f[i - 1][j])
                    ctx.set(pe.h_left, h[i][j - 1])
                    ctx.set(pe.e_left, e[i][j - 1])
                    ctx.set(pe.match, int(query[i - 1] == ref[j - 1]))
                    h[i][j] = ctx.get(pe.h)
                    e[i][j] = ctx.get(pe.e)
                    f[i][j] = ctx.get(pe.f)
                    tb[i][j] = ctx.get(pe.tb)
                    if h[i][j] > best:
                        best, best_cell = h[i][j], (i, j)

            expected = gotoh_best_score(
                query, ref, match_score=pe.match_score,
                mismatch_penalty=pe.mismatch_penalty,
                gap_open=pe.gap_open, gap_extend=pe.gap_extend)
            assert best == expected, (
                f"pair {pair_index}: PE matrix best {best} != Gotoh {expected}\n"
                f"query={query}\nref  ={ref}")

            walked = _walk_traceback(query, ref, tb, best_cell, pe)
            assert walked == best, (
                f"pair {pair_index}: traceback re-scores to {walked}, "
                f"matrix best is {best}\nquery={query}\nref  ={ref}")

    sim.add_testbench(process)
    sim.run()
    return len(pairs)


def check_score_width():
    needed = min_score_width(DEFAULT_L_MAX, DEFAULT_MATCH_SCORE)
    assert needed == DEFAULT_SCORE_WIDTH, (
        f"L_max {DEFAULT_L_MAX} needs {needed} bits, "
        f"PE default is {DEFAULT_SCORE_WIDTH}")
    assert DEFAULT_MATCH_SCORE * DEFAULT_L_MAX < (1 << DEFAULT_SCORE_WIDTH)
    return needed


if __name__ == "__main__":
    print(f"PASS: score width {check_score_width()} bits covers L_max {DEFAULT_L_MAX}")
    print(f"PASS: {run_testbench()} directed/random PE vectors verified")
    print(f"PASS: {run_exhaustive()} exhaustive 3-bit PE vectors verified")
    print(f"PASS: {run_matrix_check()} DP matrices match Gotoh, traceback re-scores correctly")