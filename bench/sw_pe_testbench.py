"""Standalone 200-vector simulation testbench for ``SmithWatermanPE``.

Run after activating the environment with:

    python sw_pe_testbench.py
"""

from random import Random

from amaranth.sim import Simulator

from sw_pe import SmithWatermanPE


def reference_score(*, diag, up, left, match, match_score,
                    mismatch_penalty, gap_penalty, score_width):
    """Software implementation of the PE recurrence."""
    score_limit = (1 << score_width) - 1
    diagonal = min(score_limit, diag + match_score) if match else max(
        0, diag - mismatch_penalty
    )
    vertical = max(0, up - gap_penalty)
    horizontal = max(0, left - gap_penalty)
    return max(diagonal, vertical, horizontal)


def run_testbench(vector_count=200, seed=0x5A17):
    """Run ``vector_count`` deterministic input vectors and return its count.

    A fixed seed keeps failures reproducible while the vectors cover values over
    the full input range, both match states, underflow, and saturation cases.
    """
    pe = SmithWatermanPE()
    rng = Random(seed)
    vectors = [
        {
            "diag": rng.randrange(1 << pe.score_width),
            "up": rng.randrange(1 << pe.score_width),
            "left": rng.randrange(1 << pe.score_width),
            "match": rng.randrange(2),
        }
        for _ in range(vector_count)
    ]

    sim = Simulator(pe)

    async def process(ctx):
        for index, vector in enumerate(vectors):
            ctx.set(pe.diag, vector["diag"])
            ctx.set(pe.up, vector["up"])
            ctx.set(pe.left, vector["left"])
            ctx.set(pe.match, vector["match"])

            actual = ctx.get(pe.score)
            expected = reference_score(
                **vector,
                match_score=pe.match_score,
                mismatch_penalty=pe.mismatch_penalty,
                gap_penalty=pe.gap_penalty,
                score_width=pe.score_width,
            )
            assert actual == expected, (
                f"vector {index}: expected {expected}, got {actual}; "
                f"inputs={vector}"
            )

    sim.add_testbench(process)
    sim.run()
    return vector_count


if __name__ == "__main__":
    count = run_testbench()
    print(f"PASS: {count} Smith-Waterman PE vectors verified")
