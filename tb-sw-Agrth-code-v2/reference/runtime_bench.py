"""Runtime of ordinary SW against attention variants, with results saved on every run.

Every version is registered in VARIANTS. A run times each registered version on
the same simulated reads, checks it against ordinary SW, and saves the results.
To add a new version, write a function that returns a Result and register it.

Registered versions (all plain Python in the same style):
  ordinary         full table with traceback
  attention        full vote (every diagonal, every base) + window SW
  attention_bits   bit-parallel vote: the same counts, computed with whole-integer
                   AND and popcount, as XNOR/popcount hardware would
  attention_prune  branch-and-bound vote: stops scoring a diagonal as soon as it
                   can no longer beat the best one found so far

All three votes choose exactly the same diagonal (the run checks this), so the
window SW and its result are the same; only the vote's cost differs.

Python timings are a software baseline. They say nothing direct about the FPGA.

Run:  python reference/runtime_bench.py                       (all versions)
      python reference/runtime_bench.py --variants ordinary attention_prune -n 50
Out:  results/runtime_<date>_<time>.csv    one row per read and version
      results/runtime_<date>_<time>.json   settings, machine, summary per version
"""

import argparse
import csv
import json
import platform
import random
import statistics
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from band_study import INDEL_LENGTHS, LOCI, ROOT, read_fasta, simulate_read
from sw_reference import GAP, MATCH, MISMATCH, STOP, DIAG, UP, LEFT, sw_align


@dataclass
class Result:
    score: int
    path: list
    vote_ns: int = 0          # time spent choosing the diagonal (0 for ordinary SW)
    align_ns: int = 0         # time spent in SW
    cells: int = 0            # SW cells computed
    comparisons: int = 0      # letter comparisons made by the vote
    d0: object = None         # chosen diagonal (None for ordinary SW)


# ------------------------------------------------------------------ votes
# Each vote returns (diagonal, comparisons). They all maximise
# Z(k) = matches - lam * mismatches and break ties toward the smallest k.

def vote_full(read, ref, lam):
    m, n = len(read), len(ref)
    best_k, best_z, comps = 0, None, 0
    for k in range(-(m - 1), n):
        i0, i1 = max(0, -k), min(m, n - k)
        v = 0
        for i in range(i0, i1):
            if read[i] == ref[i + k]:
                v += 1
        comps += i1 - i0
        z = (1 + lam) * v - lam * (i1 - i0)
        if best_z is None or z > best_z:
            best_k, best_z = k, z
    return best_k, comps


def letter_masks(seq):
    """One integer per letter; bit i is set where seq[i] is that letter."""
    masks = {c: 0 for c in "ACGT"}
    for i, c in enumerate(seq):
        masks[c] |= 1 << i
    return masks


def vote_bits(read, ref, lam):
    """Same counts as vote_full, from AND + popcount on whole integers."""
    m, n = len(read), len(ref)
    R, B = letter_masks(read), letter_masks(ref)
    best_k, best_z = 0, None
    for k in range(-(m - 1), n):
        if k >= 0:
            v = sum((R[c] & (B[c] >> k)).bit_count() for c in "ACGT")
        else:
            v = sum(((R[c] >> -k) & B[c]).bit_count() for c in "ACGT")
        N = min(m, n - k) - max(0, -k)
        z = (1 + lam) * v - lam * N
        if best_z is None or z > best_z:
            best_k, best_z = k, z
    return best_k, m * n          # the hardware equivalent still compares every pair


def vote_prune(read, ref, lam):
    """Exact branch and bound.

    A match adds 1 to Z and a mismatch subtracts lam, so a diagonal with r
    letters left to compare can gain at most r more. Diagonals are tried longest
    overlap first; a diagonal is abandoned as soon as its best possible final
    score cannot beat the leader (ties go to the smaller k, as in vote_full).
    """
    m, n = len(read), len(ref)
    order = sorted(range(-(m - 1), n), key=lambda k: (-(min(m, n - k) - max(0, -k)), k))
    best_k, best_z, comps = None, None, 0
    for k in order:
        i0, i1 = max(0, -k), min(m, n - k)
        N = i1 - i0
        if best_z is not None and (N, -k) <= (best_z, -best_k):
            continue                                   # even all matches would not win
        z, done = 0, True
        for i in range(i0, i1):
            comps += 1
            z += 1 if read[i] == ref[i + k] else -lam
            if best_z is not None and (z + (i1 - i - 1), -k) <= (best_z, -best_k):
                done = False                           # cannot catch up any more
                break
        if done and (best_z is None or (z, -k) > (best_z, -best_k)):
            best_k, best_z = k, z
    return best_k, comps


# ------------------------------------------------------------------ window SW

def window_sw(read, ref, d0, w):
    """SW with traceback on cells with |j - i - d0| <= w; cells outside read as 0.

    Same recurrence and tie-breaking as sw_reference.sw_align.
    """
    m, n = len(read), len(ref)
    H = [[0] * (n + 1) for _ in range(m + 1)]
    P = [[STOP] * (n + 1) for _ in range(m + 1)]
    best, best_i, best_j, cells = 0, 0, 0, 0
    for i in range(1, m + 1):
        a = read[i - 1]
        Hi, Hp, Pi = H[i], H[i - 1], P[i]
        for j in range(max(1, i + d0 - w), min(n, i + d0 + w) + 1):
            diag = Hp[j - 1] + (MATCH if a == ref[j - 1] else MISMATCH)
            up = Hp[j] + GAP
            left = Hi[j - 1] + GAP
            score, move = 0, STOP
            if diag > score:
                score, move = diag, DIAG
            if up > score:
                score, move = up, UP
            if left > score:
                score, move = left, LEFT
            Hi[j] = score
            Pi[j] = move
            cells += 1
            if score > best:
                best, best_i, best_j = score, i, j
    path, i, j = [], best_i, best_j
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
    return best, path, cells


# ------------------------------------------------------------------ registry

def clock(fn, *args):
    t0 = time.perf_counter_ns()
    out = fn(*args)
    return time.perf_counter_ns() - t0, out


def ordinary(read, ref, opts):
    t, (score, path) = clock(sw_align, read, ref)
    return Result(score, path, align_ns=t, cells=len(read) * len(ref))


def attention_with(vote_fn):
    def run(read, ref, opts):
        tv, (d0, comps) = clock(vote_fn, read, ref, opts.lam)
        ta, (score, path, cells) = clock(window_sw, read, ref, d0, opts.half)
        return Result(score, path, vote_ns=tv, align_ns=ta, cells=cells, comparisons=comps, d0=d0)
    return run


VARIANTS = {
    "ordinary": ordinary,
    "attention": attention_with(vote_full),
    "attention_bits": attention_with(vote_bits),
    "attention_prune": attention_with(vote_prune),
}


def best_of(fn, repeat, *args):
    """Run `repeat` times and keep the fastest timings."""
    best = None
    for _ in range(repeat):
        r = fn(*args)
        if best is None or r.vote_ns + r.align_ns < best.vote_ns + best.align_ns:
            best = r
    return best


# ------------------------------------------------------------------ run and save

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", nargs="+", default=list(VARIANTS), choices=list(VARIANTS))
    ap.add_argument("-n", type=int, default=30, help="reads per locus per condition")
    ap.add_argument("--half", type=int, default=8, help="window half-width w (W = 2w + 1)")
    ap.add_argument("--lam", type=int, default=1, help="vote penalty: 1 = matches - mismatches, 0 = matches")
    ap.add_argument("--repeat", type=int, default=1, help="run each version this many times, keep the fastest")
    ap.add_argument("--seed", type=int, default=0)
    opts = ap.parse_args()
    names = ["ordinary"] + [v for v in opts.variants if v != "ordinary"]   # ordinary is the yardstick

    rng = random.Random(opts.seed)
    regions = read_fasta(ROOT / "data" / "regions.fasta")
    conditions = [("none", 0)] + [(k, L) for k in ("del", "ins") for L in INDEL_LENGTHS]
    for _ in range(3):
        sw_align("ACGTACGTAC" * 3, "ACGTACGTAC" * 4)                       # warm-up

    rows, read_id = [], 0
    for locus, (centre, half) in LOCI.items():
        region = regions[locus]
        ref = region[centre - half:centre + half]
        for kind, length in conditions:
            for _ in range(opts.n):
                read, _ = simulate_read(rng, region, centre, half, kind, length)
                results = {name: best_of(VARIANTS[name], opts.repeat, read, ref, opts) for name in names}
                base = results["ordinary"]
                d0_full = results["attention"].d0 if "attention" in results else None
                for name, r in results.items():
                    rows.append({
                        "read": read_id, "locus": locus, "kind": kind, "indel_len": length, "variant": name,
                        "total_us": (r.vote_ns + r.align_ns) / 1000, "vote_us": r.vote_ns / 1000,
                        "align_us": r.align_ns / 1000, "score": r.score,
                        "same_score": r.score == base.score, "same_path": r.path == base.path,
                        "cells": r.cells, "vote_comparisons": r.comparisons, "d0": r.d0,
                        "same_d0_as_full_vote": (r.d0 == d0_full) if (r.d0 is not None and d0_full is not None) else "",
                    })
                read_id += 1

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    csv_path, json_path = out_dir / f"runtime_{stamp}.csv", out_dir / f"runtime_{stamp}.json"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    def summarise(sub):
        base_t = statistics.fmean(r["total_us"] for r in rows if r["variant"] == "ordinary"
                                  and r["indel_len"] in {x["indel_len"] for x in sub})
        t = statistics.fmean(r["total_us"] for r in sub)
        out = {
            "reads": len(sub), "total_us": round(t, 1),
            "vote_us": round(statistics.fmean(r["vote_us"] for r in sub), 1),
            "align_us": round(statistics.fmean(r["align_us"] for r in sub), 1),
            "speedup_vs_ordinary": round(base_t / t, 2),
            "same_score_share": round(sum(r["same_score"] for r in sub) / len(sub), 4),
            "cells": round(statistics.fmean(r["cells"] for r in sub)),
            "vote_comparisons": round(statistics.fmean(r["vote_comparisons"] for r in sub)),
        }
        flags = [r["same_d0_as_full_vote"] for r in sub if r["same_d0_as_full_vote"] != ""]
        if flags:
            out["same_d0_as_full_vote_share"] = round(sum(flags) / len(flags), 4)
        return out

    summary = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "settings": {"variants": names, "reads_per_condition": opts.n, "half_width": opts.half,
                     "window_cells": 2 * opts.half + 1, "vote_lambda": opts.lam, "repeat": opts.repeat,
                     "seed": opts.seed, "scoring": {"match": MATCH, "mismatch": MISMATCH, "gap": GAP}},
        "machine": {"python": platform.python_version(), "system": platform.platform(),
                    "processor": platform.processor()},
        "variants": {}, "csv": csv_path.name,
    }
    for name in names:
        sub = [r for r in rows if r["variant"] == name]
        summary["variants"][name] = {
            "overall": summarise(sub),
            "by_indel_length": {str(L): summarise([r for r in sub if r["indel_len"] == L])
                                for L in [0] + INDEL_LENGTHS},
        }
    json_path.write_text(json.dumps(summary, indent=2))
    md_path = out_dir / f"runtime_{stamp}.md"
    md_path.write_text(report(summary, read_id, csv_path.name, json_path.name), encoding="utf-8")

    print(f"{read_id} reads, window W = {2 * opts.half + 1}, vote lambda = {opts.lam}, repeat = {opts.repeat}")
    print(f"{'version':>16s} {'total':>9s} {'vote':>8s} {'SW':>8s} {'speedup':>8s} {'same score':>11s} "
          f"{'same diagonal':>14s} {'vote comparisons':>17s}")
    for name in names:
        o = summary["variants"][name]["overall"]
        same_d0 = f"{o['same_d0_as_full_vote_share']:.1%}" if "same_d0_as_full_vote_share" in o else "-"
        print(f"{name:>16s} {o['total_us']:>8.0f}u {o['vote_us']:>7.0f}u {o['align_us']:>7.0f}u "
              f"{o['speedup_vs_ordinary']:>7.2f}x {o['same_score_share']:>10.1%} {same_d0:>14s} {o['vote_comparisons']:>17d}")
    print("saved", csv_path, "and", json_path.name, "and", md_path.name)


def report(summary, n_reads, csv_name, json_name):
    """The run's tables as Markdown: speed and agreement, overall and by indel length."""
    s, variants = summary["settings"], summary["variants"]
    names = list(variants)
    out = [
        f"# Runtime: ordinary SW vs attention ({summary['created']})",
        "",
        f"{n_reads} simulated 150 bp reads, window W = {s['window_cells']} (w = {s['half_width']}), "
        f"vote lambda = {s['vote_lambda']}, fastest of {s['repeat']} run(s), seed {s['seed']}. "
        f"Python {summary['machine']['python']} on {summary['machine']['system']}. "
        "Times are means per read in microseconds. Python timings are a software baseline, not FPGA numbers.",
        "",
        "## Overall",
        "",
        "| version | total (us) | vote (us) | SW (us) | speedup vs ordinary | same score as ordinary "
        "| same diagonal as full vote | SW cells | vote comparisons |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name in names:
        o = variants[name]["overall"]
        same_d0 = f"{o['same_d0_as_full_vote_share']:.1%}" if "same_d0_as_full_vote_share" in o else "-"
        out.append(f"| {name} | {o['total_us']:,.0f} | {o['vote_us']:,.0f} | {o['align_us']:,.0f} | "
                   f"{o['speedup_vs_ordinary']:.2f}x | {o['same_score_share']:.1%} | {same_d0} | "
                   f"{o['cells']:,} | {o['vote_comparisons']:,} |")

    lengths = list(variants[names[0]]["by_indel_length"])
    out += ["", "## Time by indel length (us, speedup in brackets)", "",
            "| indel length (bp) | " + " | ".join(names) + " |",
            "|---:|" + "---:|" * len(names)]
    for L in lengths:
        cells = []
        for name in names:
            b = variants[name]["by_indel_length"][L]
            cells.append(f"{b['total_us']:,.0f}" if name == "ordinary"
                         else f"{b['total_us']:,.0f} ({b['speedup_vs_ordinary']:.2f}x)")
        out.append(f"| {'none' if L == '0' else L} | " + " | ".join(cells) + " |")

    out += ["", "## Same score as ordinary SW, by indel length", "",
            "| indel length (bp) | " + " | ".join(n for n in names if n != "ordinary") + " |",
            "|---:|" + "---:|" * (len(names) - 1)]
    for L in lengths:
        out.append(f"| {'none' if L == '0' else L} | " + " | ".join(
            f"{variants[n]['by_indel_length'][L]['same_score_share']:.1%}" for n in names if n != "ordinary") + " |")

    out += ["", f"Per-read data: `{csv_name}`. All numbers: `{json_name}`.", ""]
    return "\n".join(out)


if __name__ == "__main__":
    main()
