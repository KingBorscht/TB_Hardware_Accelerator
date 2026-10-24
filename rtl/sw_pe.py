from amaranth import Cat, Const, Elaboratable, Module, Mux, Signal

DEFAULT_MATCH_SCORE = 2
DEFAULT_MISMATCH_PENALTY = 4
DEFAULT_GAP_OPEN = 4
DEFAULT_GAP_EXTEND = 2
DEFAULT_SCORE_WIDTH = 13
DEFAULT_L_MAX = 3072

TB_WIDTH = 4
TB_STOP = 0
TB_DIAG = 1
TB_LEFT = 2
TB_UP = 3
TB_E_EXTEND_BIT = 2
TB_F_EXTEND_BIT = 3


def min_score_width(max_read_length, match_score=DEFAULT_MATCH_SCORE):
    return (match_score * max_read_length).bit_length()


class SmithWatermanPE(Elaboratable):
    def __init__(self, *, score_width=DEFAULT_SCORE_WIDTH,
                 match_score=DEFAULT_MATCH_SCORE,
                 mismatch_penalty=DEFAULT_MISMATCH_PENALTY,
                 gap_open=DEFAULT_GAP_OPEN,
                 gap_extend=DEFAULT_GAP_EXTEND):
        if score_width < 1:
            raise ValueError("score_width must be at least 1")
        if min(match_score, mismatch_penalty, gap_open, gap_extend) < 0:
            raise ValueError("scores and penalties must be non-negative")
        if gap_open + gap_extend >= (1 << score_width):
            raise ValueError("gap_open + gap_extend must fit in score_width bits")

        self.score_width = score_width
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_open = gap_open
        self.gap_extend = gap_extend

        self.h_diag = Signal(score_width)
        self.h_up = Signal(score_width)
        self.f_up = Signal(score_width)
        self.h_left = Signal(score_width)
        self.e_left = Signal(score_width)
        self.match = Signal()

        self.h = Signal(score_width)
        self.e = Signal(score_width)
        self.f = Signal(score_width)
        self.tb = Signal(TB_WIDTH)

    def _sat_sub(self, m, value, amount, name):
        out = Signal(self.score_width, name=name)
        if amount == 0:
            m.d.comb += out.eq(value)
        else:
            m.d.comb += out.eq(Mux(value >= amount, value - amount, 0))
        return out

    def elaborate(self, platform):
        m = Module()
        width = self.score_width
        limit = (1 << width) - 1
        # a gap of length k costs gap_open + k * gap_extend
        gap_first = self.gap_open + self.gap_extend

        diag_sum = Signal(width + 1)
        diag_cand = Signal(width)
        m.d.comb += diag_sum.eq(self.h_diag + self.match_score)
        m.d.comb += diag_cand.eq(
            Mux(
                self.match,
                Mux(diag_sum > limit, Const(limit, width), diag_sum),
                Mux(self.h_diag >= self.mismatch_penalty,
                    self.h_diag - self.mismatch_penalty, 0),
            )
        )

        e_open = self._sat_sub(m, self.h_left, gap_first, "e_open")
        e_ext = self._sat_sub(m, self.e_left, self.gap_extend, "e_ext")
        e_extended = Signal()
        m.d.comb += e_extended.eq(e_ext > e_open)
        m.d.comb += self.e.eq(Mux(e_extended, e_ext, e_open))

        f_open = self._sat_sub(m, self.h_up, gap_first, "f_open")
        f_ext = self._sat_sub(m, self.f_up, self.gap_extend, "f_ext")
        f_extended = Signal()
        m.d.comb += f_extended.eq(f_ext > f_open)
        m.d.comb += self.f.eq(Mux(f_extended, f_ext, f_open))

        e_wins = Signal()
        diag_wins = Signal()
        all_zero = Signal()
        best_ef = Signal(width)
        m.d.comb += e_wins.eq(self.e >= self.f)
        m.d.comb += best_ef.eq(Mux(e_wins, self.e, self.f))
        m.d.comb += diag_wins.eq(diag_cand >= best_ef)
        m.d.comb += all_zero.eq((diag_cand == 0) & (best_ef == 0))
        m.d.comb += self.h.eq(Mux(diag_wins, diag_cand, best_ef))

        # tie-break order: diag, then left, then up; h == 0 stops.
        # Decided from the comparisons above so it doesn't wait on h.
        source = Signal(2)
        with m.If(all_zero):
            m.d.comb += source.eq(TB_STOP)
        with m.Elif(diag_wins):
            m.d.comb += source.eq(TB_DIAG)
        with m.Elif(e_wins):
            m.d.comb += source.eq(TB_LEFT)
        with m.Else():
            m.d.comb += source.eq(TB_UP)

        m.d.comb += self.tb.eq(Cat(source, e_extended, f_extended))

        return m