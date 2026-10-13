"""Smith-Waterman processing element written in Amaranth.

The surrounding band controller supplies only valid neighbour scores for each
anti-diagonal. This element performs the recurrence for one cell.
"""

from amaranth import Const, Elaboratable, Module, Mux, Signal


class SmithWatermanPE(Elaboratable):
    """One local-alignment Smith-Waterman processing element.

    Inputs are the three preceding-cell scores and whether the current pair of
    sequence symbols matches. ``score`` is the local-alignment score for this
    cell:

        max(0, diag + substitution, up - gap, left - gap)
    """

    def __init__(self, *, score_width=8, match_score=2,
                 mismatch_penalty=1, gap_penalty=1):
        if score_width < 1:
            raise ValueError("score_width must be at least 1")
        if match_score < 0 or mismatch_penalty < 0 or gap_penalty < 0:
            raise ValueError("scores and penalties must be non-negative")

        self.score_width = score_width
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty

        self.diag = Signal(score_width)
        self.up = Signal(score_width)
        self.left = Signal(score_width)
        self.match = Signal()
        self.score = Signal(score_width)

    def elaborate(self, platform):
        m = Module()
        limit = (1 << self.score_width) - 1

        diagonal_sum = Signal(self.score_width + 1)
        diagonal = Signal(self.score_width)
        vertical = Signal(self.score_width)
        horizontal = Signal(self.score_width)
        best = Signal(self.score_width)

        m.d.comb += diagonal_sum.eq(self.diag + self.match_score)
        m.d.comb += diagonal.eq(
            Mux(
                self.match,
                Mux(diagonal_sum > limit, Const(limit, self.score_width), diagonal_sum),
                Mux(self.diag >= self.mismatch_penalty,
                    self.diag - self.mismatch_penalty, 0),
            )
        )
        m.d.comb += vertical.eq(
            Mux(self.up >= self.gap_penalty, self.up - self.gap_penalty, 0)
        )
        m.d.comb += horizontal.eq(
            Mux(self.left >= self.gap_penalty, self.left - self.gap_penalty, 0)
        )
        m.d.comb += best.eq(Mux(diagonal >= vertical, diagonal, vertical))
        m.d.comb += self.score.eq(Mux(best >= horizontal, best, horizontal))

        return m
