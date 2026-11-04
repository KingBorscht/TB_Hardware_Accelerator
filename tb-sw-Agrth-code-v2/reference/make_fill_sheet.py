"""Make a "fill in the window" worksheet with a random reference and read.

Same layout and wording as docs/fill_table.tex. The read is the reference with
one letter deleted. Cells outside the window (no attention) are printed black;
only the white cells are filled in.

Run:  python reference/make_fill_sheet.py --seed 1      (from the project root)
Out:  docs/fill_table_window.tex   (page 1 worksheet, page 2 answers)
"""

import argparse
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MATCH, MISMATCH, GAP = 2, -1, -2
W_HALF = 1                      # window half-width: 3 diagonals
ORDINAL = ["first", "second", "third", "fourth"]


def table(read, ref, d0=None, w=None):
    """SW table; with d0 and w only cells inside the window are computed."""
    m, n = len(read), len(ref)
    H = [[0] * (n + 1) for _ in range(m + 1)]
    P = [[0] * (n + 1) for _ in range(m + 1)]
    best = (0, 0, 0)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if d0 is not None and abs(j - i - d0) > w:
                continue
            cands = [(H[i - 1][j - 1] + (MATCH if read[i - 1] == ref[j - 1] else MISMATCH), 1),
                     (H[i - 1][j] + GAP, 2),
                     (H[i][j - 1] + GAP, 3)]
            s, p = 0, 0
            for val, ptr in cands:
                if val > s:
                    s, p = val, ptr
            H[i][j], P[i][j] = s, p
            if s > best[0]:
                best = (s, i, j)
    return H, P, best


def path_of(P, i, j):
    out = []
    while i > 0 and j > 0 and P[i][j]:
        out.append((i, j, P[i][j]))
        if P[i][j] == 1:
            i, j = i - 1, j - 1
        elif P[i][j] == 2:
            i -= 1
        else:
            j -= 1
    return out[::-1]


def vote(read, ref):
    """Matches minus mismatches on every diagonal."""
    m, n = len(read), len(ref)
    z = {}
    for k in range(-(m - 1), n):
        pairs = [(i, i + k) for i in range(m) if 0 <= i + k < n]
        v = sum(read[i] == ref[j] for i, j in pairs)
        z[k] = 2 * v - len(pairs)
    return z


def pick(rng, n=9):
    """A reference and read for which the worksheet works cleanly."""
    while True:
        ref = "".join(rng.choice("ACGT") for _ in range(n))
        p = rng.randrange(3, n - 3)
        if ref[p] in (ref[p - 1], ref[p + 1]):
            continue                                    # keep the deleted letter unambiguous
        if any(ref[i] == ref[i + 1] == ref[i + 2] for i in range(n - 2)):
            continue
        read = ref[:p] + ref[p + 1:]
        z = vote(read, ref)
        top = max(z.values())
        if sum(1 for v in z.values() if v == top) != 1:
            continue                                    # the vote must have one clear winner
        d0 = max(z, key=z.get)
        Hf, Pf, bf = table(read, ref)
        Hw, Pw, bw = table(read, ref, d0, W_HALF)
        if bf[0] != bw[0] or bf[0] != MATCH * len(read) + GAP:
            continue
        if sum(1 for row in Hf for x in row if x == bf[0]) != 1:
            continue
        pw = path_of(Pw, bw[1], bw[2])
        steps = [s for _, _, s in pw]
        if steps.count(3) != 1 or 2 in steps or steps[0] == 3 or steps[-1] == 3:
            continue
        if path_of(Pf, bf[1], bf[2]) != pw:
            continue                                    # the full table must agree with the window
        return ref, read, p, d0, Hw, bw, pw


def tex_table(read, ref, cell):
    """cell(i, j) returns the LaTeX for one cell."""
    n = len(ref)
    cols = "|C|C||" + "C|" * n
    head = r"\rowh &  & " + " & ".join(r"\textbf{%s}" % c for c in ref) + r" \\ \hline"
    zero = r"\rowh & \textit{0} & " + " & ".join([r"\textit{0}"] * n) + r" \\ \hline\hline"
    lines = [r"\begin{tabular}{%s}" % cols, r"\hline", head, zero]
    for i in range(1, len(read) + 1):
        cells = " & ".join(cell(i, j) for j in range(1, n + 1))
        lines.append(r"\rowh \textbf{%s} & \textit{0} & %s \\ \hline" % (read[i - 1], cells))
    lines.append(r"\end{tabular}")
    return "\n".join(lines)


def explain(read, ref, H, j, left_black):
    """One worked cell of row 1, in the wording of fill_table.tex."""
    same = read[0] == ref[j - 1]
    up_left = ("upper-left is 0, same letters, so $0 + 2 = 2$" if same
               else "upper-left gives $0 - 1 = -1$")
    if left_black:
        left = "the cell to the left is black and counts as 0, so left gives $-2$"
    else:
        left = "left gives $%d - 2 = %d$" % (H[1][j - 1], H[1][j - 1] - 2)
    return "%s against %s: %s. Above gives $-2$, %s. Write %d." % (
        read[0], ref[j - 1], up_left, left, H[1][j])


TEMPLATE = r"""\documentclass[12pt]{article}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[margin=0.75in]{geometry}
\usepackage[table]{xcolor}
\usepackage{array}
\usepackage{enumitem}
\pagestyle{empty}

\newcolumntype{C}{>{\centering\arraybackslash}p{1.15cm}}
\newcommand{\rowh}{\rule{0pt}{0.7cm}}
\newcommand{\off}{\cellcolor{black!88}}

\begin{document}

\begin{center}
{\LARGE\bfseries Fill in the table, with attention}\\[4pt]
{\large Same rule. This time most of the table is switched off.}
\end{center}

\medskip
\noindent
The reference is \textbf{@REF@}. The read is \textbf{@READ@}. One letter is
missing from the read. We first slid the read along the reference to see where
the letters line up best, and kept only the cells near that place.
\textbf{Black cells are switched off: leave them empty, and count them as 0.}

\subsection*{The rule for one cell}

Look at three neighbours and work out three numbers:

\begin{center}
\renewcommand{\arraystretch}{1.4}
\begin{tabular}{@{}lll@{}}
from the \textbf{upper-left} cell: & same letters: add 2 & different letters: subtract 1 \\
from the cell \textbf{above}: & subtract 2 & \\
from the cell to the \textbf{left}: & subtract 2 & \\
\end{tabular}
\end{center}

\noindent
Write the \textbf{biggest} of the three. If the biggest is below 0, write
\textbf{0}.

\medskip
\noindent
The border, in italics, is all zeros. Two cells are done for you:
\begin{itemize}[nosep]
\item @EX1@
\item @EX2@
\end{itemize}

\medskip
\begin{center}
\renewcommand{\arraystretch}{1}
@BLANK@
\end{center}

\subsection*{When the white cells are full}

\begin{enumerate}[nosep]
\item What is the biggest number in the table, and where is it?
\item Start there. Ask which neighbour that number came from, and step to it.
      Repeat until you reach a 0. Draw the path.
\item One step of your path goes \textbf{left} instead of diagonally. That is
      the missing letter. Which one is it?
\item How many cells did you fill in? How many would the whole table have?
\end{enumerate}

\newpage

\begin{center}
{\LARGE\bfseries Answer}
\end{center}

\bigskip
\begin{center}
\renewcommand{\arraystretch}{1}
@WINDOW@
\end{center}

\bigskip
\begin{enumerate}
\item The biggest number is \textbf{@SCORE@}, in the bottom-right cell.
\item The path is the boxed cells: $@CHAIN@$.
\item @LEFTSTEP@ So the read skips the @WHICH@ of the reference:

\begin{center}
\ttfamily
\begin{tabular}{@{}l@{\ \ }l@{}}
reference: & @ALNREF@ \\
read:      & @ALNREAD@
\end{tabular}
\end{center}

\item @CELLS@ cells. The whole table has @ALL@.
\end{enumerate}

\end{document}
"""


BEGIN_MARK = "% --- attention sheet: generated by reference/make_fill_sheet.py ---"
END_MARK = "% --- end of attention sheet ---"


def merge_into(target, text):
    """Put the two generated pages at the end of an existing worksheet file."""
    src = Path(target).read_text(encoding="utf-8")
    body = text.split(r"\begin{document}", 1)[1].rsplit(r"\end{document}", 1)[0].strip()
    block = "\n".join([BEGIN_MARK, r"\newpage", r"\renewcommand{\rowh}{\rule{0pt}{0.7cm}}", "", body, END_MARK])
    if BEGIN_MARK in src:                                   # replace an earlier generated block
        head, rest = src.split(BEGIN_MARK, 1)
        src = head + block + rest.split(END_MARK, 1)[1]
    else:
        head, tail = src.rsplit(r"\end{document}", 1)
        src = head.rstrip() + "\n\n" + block + "\n\n" + r"\end{document}" + tail
    for line in (r"\usepackage[table]{xcolor}", r"\newcommand{\off}{\cellcolor{black!88}}"):
        if line not in src:                                 # what the black cells need
            src = src.replace(r"\begin{document}", line + "\n\n" + r"\begin{document}", 1)
    Path(target).write_text(src, encoding="utf-8", newline="\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=str(ROOT / "docs" / "fill_table_window.tex"))
    ap.add_argument("--into", help="append the two pages to this existing .tex file instead of writing --out")
    args = ap.parse_args()
    rng = random.Random(args.seed)
    ref, read, p, d0, H, best, path = pick(rng)
    m, n = len(read), len(ref)
    inside = lambda i, j: abs(j - i - d0) <= W_HALF
    on_path = {(i, j) for i, j, _ in path}

    row1 = [j for j in range(1, n + 1) if inside(1, j)]
    j1, j2 = row1[0], row1[1]
    done = {(1, j1), (1, j2)}

    def blank(i, j):
        if not inside(i, j):
            return r"\off"
        return str(H[i][j]) if (i, j) in done else ""

    def filled(i, j):
        if not inside(i, j):
            return r"\off"
        return r"\fbox{%d}" % H[i][j] if (i, j) in on_path else str(H[i][j])

    # the one step that goes left, described as in fill_table.tex
    k = next(idx for idx, (_, _, s) in enumerate(path) if s == 3)
    (ia, ja, _), (ib, jb, _), (ic, jc, _) = path[k + 1], path[k], path[k - 1]
    left_step = ("The step from the boxed %d in row %s back to the boxed %d in row %s goes diagonally, "
                 "and the step from that %d back to the boxed %d in row %s goes \\textbf{left}."
                 % (H[ia][ja], read[ia - 1], H[ib][jb], read[ib - 1], H[ib][jb], H[ic][jc], read[ic - 1]))
    nth = ref[:p + 1].count(ref[p])
    which = "%s \\textbf{%s}" % (ORDINAL[nth - 1], ref[p])
    chain = r" \leftarrow ".join(str(H[i][j]) for i, j, _ in reversed(path))

    cells = sum(1 for i in range(1, m + 1) for j in range(1, n + 1) if inside(i, j))
    fill = {
        "@REF@": ref, "@READ@": read,
        "@EX1@": explain(read, ref, H, j1, True), "@EX2@": explain(read, ref, H, j2, False),
        "@BLANK@": tex_table(read, ref, blank), "@WINDOW@": tex_table(read, ref, filled),
        "@SCORE@": str(best[0]), "@CHAIN@": chain, "@LEFTSTEP@": left_step, "@WHICH@": which,
        "@ALNREF@": " ".join(ref), "@ALNREAD@": " ".join(read[:p] + "_" + read[p:]).replace("_", r"\_"),
        "@CELLS@": str(cells), "@ALL@": str(m * n),
    }
    text = TEMPLATE
    for key, val in fill.items():
        text = text.replace(key, val)
    print("reference", ref, "| read", read, "| deleted", which.replace("\\textbf{", "").replace("}", ""))
    print("vote picks diagonal", d0, "| score", best[0], "| cells", cells, "of", m * n)
    if args.into:
        merge_into(args.into, text)
        print("added to", args.into)
    else:
        Path(args.out).write_text(text, encoding="utf-8", newline="\n")
        print("wrote", args.out)


if __name__ == "__main__":
    main()
