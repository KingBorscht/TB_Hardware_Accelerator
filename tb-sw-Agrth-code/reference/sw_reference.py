"""Smith-Waterman reference model (plain Python).

Local alignment, linear gap penalty, with traceback. This is the correctness
baseline: slow, simple, and deterministic. Tie-breaking is fixed and documented
because the hardware has to reproduce it bit-for-bit.

Conventions
    read index i = 0..m-1 (rows), reference index j = 0..n-1 (columns)
    diagonal d = j - i
    H has an extra zero row/column, so H[i+1][j+1] is the cell for (i, j)

Tie-breaking
    cell:      diagonal > up > left   (stop wins only when the score is 0)
    end cell:  first maximum in row-major order (smallest i, then smallest j)
"""

MATCH = 2
MISMATCH = -1
GAP = -2

STOP, DIAG, UP, LEFT = 0, 1, 2, 3


def sw_align(read, ref, match=MATCH, mismatch=MISMATCH, gap=GAP):
    """Full Smith-Waterman. Returns (score, path).

    path is the list of (i, j, move) cells from the start of the alignment to
    its end, where move is DIAG, UP (read base against a gap) or LEFT
    (reference base against a gap).
    """
    m, n = len(read), len(ref)
    H = [[0] * (n + 1) for _ in range(m + 1)]
    P = [[STOP] * (n + 1) for _ in range(m + 1)]
    best, best_i, best_j = 0, 0, 0

    for i in range(1, m + 1):
        a = read[i - 1]
        Hi, Hp, Pi = H[i], H[i - 1], P[i]
        for j in range(1, n + 1):
            diag = Hp[j - 1] + (match if a == ref[j - 1] else mismatch)
            up = Hp[j] + gap
            left = Hi[j - 1] + gap
            score, move = 0, STOP
            if diag > score:
                score, move = diag, DIAG
            if up > score:
                score, move = up, UP
            if left > score:
                score, move = left, LEFT
            Hi[j] = score
            Pi[j] = move
            if score > best:
                best, best_i, best_j = score, i, j

    path = []
    i, j = best_i, best_j
    while P[i][j] != STOP:
        move = P[i][j]
        path.append((i - 1, j - 1, move))
        if move == DIAG:
            i, j = i - 1, j - 1
        elif move == UP:
            i -= 1
        else:
            j -= 1
    path.reverse()
    return best, path


def path_diagonals(path):
    """Diagonals (j - i) visited by a traceback path."""
    return [j - i for i, j, _ in path]


def cigar(path):
    """Compact run-length string of the path: M (diag), I (up), D (left)."""
    if not path:
        return ""
    sym = {DIAG: "M", UP: "I", LEFT: "D"}
    out, run, prev = [], 0, None
    for _, _, move in path:
        s = sym[move]
        if s == prev:
            run += 1
        else:
            if prev is not None:
                out.append(f"{run}{prev}")
            prev, run = s, 1
    out.append(f"{run}{prev}")
    return "".join(out)


if __name__ == "__main__":
    ref = "ACGTTGCAGGCTAACGTCGATCGGATC"
    for name, read in [
        ("exact", "TGCAGGCTAACGTCG"),
        ("SNP", "TGCAGGCTTACGTCG"),
        ("deletion", "TGCAGGCACGTCGATC"),
        ("insertion", "TGCAGGCTAAGGACGTCGATC"),
    ]:
        score, path = sw_align(read, ref)
        print(f"{name:10s} score={score:3d} start=(read {path[0][0]}, ref {path[0][1]}) {cigar(path)}")
