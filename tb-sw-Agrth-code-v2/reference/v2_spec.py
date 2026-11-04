"""SCOPE v2 specification model: banded affine-gap Smith-Waterman with traceback.

Scoring (SCOPE v2 decision 9, minimap2 map-ont preset): match +2, mismatch -4,
gap open 4, extend 2. A gap of length L costs OPEN + L * EXTEND, so a 1 bp gap
costs 6 and every further base 2 (minimap2's convention).

Recurrences (Gotoh), for read base i and reference base j:
  E[i][j] = max(H[i][j-1] - (OPEN + EXTEND), E[i][j-1] - EXTEND)   gap in the read  (LEFT)
  F[i][j] = max(H[i-1][j] - (OPEN + EXTEND), F[i-1][j] - EXTEND)   gap in the ref   (UP)
  H[i][j] = max(0, H[i-1][j-1] + s(i, j), E[i][j], F[i][j])

The band is one row at a time: row i computes the cells with |j - i - c_i| <= w,
where c_i is the band centre (a diagonal j - i) for that row. Cells outside the
band read as H = 0 and E = F = minus infinity, as in the v1 window.

Traceback pointer, 4 bits per cell (SCOPE v2 decision 12):
  bits 1..0  where H came from: STOP, DIAG, LEFT (from E), UP (from F)
  bit 2      E extended (1) or opened (0)
  bit 3      F extended (1) or opened (0)

Tie order for H: DIAG > LEFT > UP with LEFT = 2, UP = 3, the repo's PE order
(agreed with the club lead 2026-11-03; not yet checked against the PE code).
The v1 reference (sw_reference.py) keeps DIAG > UP > LEFT with UP = 2.

PROVISIONAL, must match the RTL PE before bit-exact diffing:
  - an E or F tie between opening and extending counts as opening.
  - the end cell is the first maximum in row-major order.

This is the slow, readable model that produces test vectors. v2_band_study.py has
a fast numpy version of the same scores.
"""

MATCH, MISMATCH, OPEN, EXTEND = 2, -4, 4, 2
STOP, DIAG, LEFT, UP = 0, 1, 2, 3
TIE_ORDER = (DIAG, LEFT, UP)
NEG = -(1 << 30)
SCORE_BITS = 13                      # SCOPE v2 decision 10


def full_band(m, n):
    """Centres and half-width that cover the whole matrix (ordinary SW)."""
    return [0] * (m + 1), m + n


def sw_affine(read, ref, centres=None, half=None, tie_order=TIE_ORDER):
    """Return (score, path, pointers, cells).

    centres[i] is the band centre for row i (index 0 unused); half is w.
    path is a list of (i, j, move) with 0-based read/ref positions, as in v1.
    pointers maps (i, j) to the 4-bit pointer of every computed cell.
    """
    m, n = len(read), len(ref)
    if centres is None:
        centres, half = full_band(m, n)
    H_prev, F_prev = {}, {}
    ptr = {}
    best, best_i, best_j, cells = 0, 0, 0, 0
    for i in range(1, m + 1):
        lo = max(1, i + centres[i] - half)
        hi = min(n, i + centres[i] + half)
        H_row, E_row, F_row = {}, {}, {}
        a = read[i - 1]
        for j in range(lo, hi + 1):
            e_open = H_row.get(j - 1, 0) - (OPEN + EXTEND)
            e_ext = E_row.get(j - 1, NEG) - EXTEND
            E, e_bit = (e_ext, 1) if e_ext > e_open else (e_open, 0)
            f_open = H_prev.get(j, 0) - (OPEN + EXTEND)
            f_ext = F_prev.get(j, NEG) - EXTEND
            F, f_bit = (f_ext, 1) if f_ext > f_open else (f_open, 0)
            D = H_prev.get(j - 1, 0) + (MATCH if a == ref[j - 1] else MISMATCH)
            cand = {DIAG: D, LEFT: E, UP: F}
            H, src = 0, STOP
            for move in tie_order:
                if cand[move] > H:
                    H, src = cand[move], move
            assert H < (1 << SCORE_BITS), "score does not fit in 13 bits"
            H_row[j], E_row[j], F_row[j] = H, E, F
            ptr[(i, j)] = src | (e_bit << 2) | (f_bit << 3)
            cells += 1
            if H > best:
                best, best_i, best_j = H, i, j
        H_prev, F_prev = H_row, F_row
    return best, traceback(ptr, best_i, best_j), ptr, cells


def traceback(ptr, i, j):
    """Follow the 4-bit pointers from (i, j) back to a STOP."""
    path, state = [], "H"
    while True:
        p = ptr.get((i, j), STOP)          # row 0, column 0 and cells outside the band stop

        if state == "H":
            src = p & 3
            if src == STOP:
                break
            if src == DIAG:
                path.append((i - 1, j - 1, DIAG))
                i, j = i - 1, j - 1
                continue
            state = "E" if src == LEFT else "F"
        if state == "E":                               # gap in the read: consume ref base j
            path.append((i - 1, j - 1, LEFT))
            state = "E" if (p >> 2) & 1 else "H"
            j -= 1
        else:                                          # gap in the ref: consume read base i
            path.append((i - 1, j - 1, UP))
            state = "F" if (p >> 3) & 1 else "H"
            i -= 1
    path.reverse()
    return path


def path_score(read, ref, path):
    """Score of a path, recomputed from the moves (checks the traceback)."""
    s, prev = 0, None
    for i, j, move in path:
        if move == DIAG:
            s += MATCH if read[i] == ref[j] else MISMATCH
        else:
            s -= EXTEND if move == prev else OPEN + EXTEND
        prev = move
    return s


if __name__ == "__main__":
    import random
    rng = random.Random(1)
    for t in range(300):
        m, n = rng.randint(1, 40), rng.randint(1, 40)
        read = "".join(rng.choice("ACGT") for _ in range(m))
        ref = "".join(rng.choice("ACGT") for _ in range(n))
        s, path, _, _ = sw_affine(read, ref)
        assert path_score(read, ref, path) == s, (read, ref)
    print("traceback score equals cell score on 300 random pairs")
