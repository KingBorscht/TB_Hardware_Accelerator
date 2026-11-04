# Proofs for the alignment kernel

This document states and proves the claims behind the proposed architecture (diagonal vote, then a sliding-window Smith-Waterman). Each result is marked as one of:

- **Proved**: a complete proof is given here.
- **Classical**: a standard result, cited.
- **Computed**: an exact number obtained from the reference sequences by exhaustive calculation, with no sampling.
- **Empirical**: supported only by simulated reads. Not a theorem.

## Summary

| # | Claim | Status |
|---|-------|--------|
| 2.1 | $H(i,j)$ is the best score of a local alignment ending at $(i,j)$ | Classical, proof sketched |
| 2.2 | $0 \le H(i,j) \le M\min(i,j)$; scores fit in 9 bits | Proved |
| 3.1 | Neighbouring cells differ by a bounded amount; 3 bits (side) and 2 bits (diagonal) suffice | Proved |
| 4.1 | Score equals $\tfrac{M}{2}(m+n)$ minus a shortest-path cost | Proved |
| 5.1 | The window computes at most $mW$ cells | Proved |
| 5.3 | Window score never exceeds the full score | Proved |
| 5.5 | Window score equals the full score when an optimal path fits in the window | Proved |
| 6.2 | The penalised vote lands on an acceptable diagonal when substitutions are below a margin | Proved |
| 6.3 | Then the window contains the true alignment | Proved |
| 6.4 | Values of the margin for the three loci | Computed |
| 6.5 | Probability bound for the vote under independent substitutions | Proved, given 6.4 |
| 7 | Storage, cell, and cycle counts | Proved (counting) |
| 8.1 | The direct vote makes exactly $mn$ base comparisons | Proved |
| 8.2 | The vote can be computed in $O((m+n)\log(m+n))$ | Classical |
| 9.2, 9.3 | Wavefront update and comparison bounds | Proved; algorithm correctness is classical |
| 10.1 | Segment-vote band computes at most $mW + 2S_{\text{seg}}G$ cells | Proved |
| 11.1 | Expected cost with a gate | Proved |
| 12 | Claims that are not proved | Listed |

## 1. Notation

- Read $a = a_1 \dots a_m$, reference window $b = b_1 \dots b_n$.
- Scores: match $M > 0$, mismatch $x \le M$, gap $g < 0$ per gapped base. The project uses $M = 2$, $x = -1$, $g = -2$, $m = 150$, $n \le 180$.
- $s(a_i, b_j) = M$ if $a_i = b_j$, otherwise $x$.
- Grid nodes are pairs $(i, j)$ with $0 \le i \le m$, $0 \le j \le n$. A node with $i, j \ge 1$ is also called a cell.
- The diagonal of a node is $j - i$.
- A **path** is a sequence of steps between nodes. A diagonal step $(i-1, j-1) \to (i, j)$ scores $s(a_i, b_j)$. A down step $(i-1, j) \to (i, j)$ and a right step $(i, j-1) \to (i, j)$ each score $g$. The score of a path is the sum of its steps. The empty path scores 0. A path is a local alignment of the substrings it spans.
- $W = 2w + 1$ is the window width and $w$ its half-width. $d_0$ is the diagonal chosen by the vote.

## 2. Smith-Waterman

**Definition.** $H(i, 0) = H(0, j) = 0$, and for $i, j \ge 1$

$$H(i,j) = \max\{\,0,\; H(i-1,j-1) + s(a_i,b_j),\; H(i-1,j) + g,\; H(i,j-1) + g\,\}.$$

$S = \max_{i,j} H(i,j)$.

**Lemma 2.1.** $H(i,j)$ is the maximum score over all paths that end at node $(i,j)$, the empty path included.

*Proof sketch.* Induction on $i + j$. A non-empty path ending at $(i,j)$ ends with one of the three steps, and the rest of it is a path ending at the corresponding predecessor. Maximising over the last step and, by the induction hypothesis, over the rest gives the three non-zero terms. The empty path gives the term 0. On the boundary every non-empty path consists of gap steps only and scores below 0. $\square$

**Proposition 2.2.** $0 \le H(i,j) \le M \min(i,j)$. Hence $0 \le S \le M\min(m,n)$, and an unsigned score needs $\lceil \log_2(M\min(m,n) + 1) \rceil$ bits.

*Proof.* The lower bound is the term 0. For the upper bound use Lemma 2.1. A diagonal step increases both coordinates by one, so a path ending at $(i,j)$ contains at most $\min(i,j)$ diagonal steps. Each scores at most $M$, and every other step scores $g < 0$. $\square$

With the project values, $S \le 300$ and 9 bits suffice.

## 3. Differences between neighbouring cells

**Proposition 3.1.** For all $1 \le i \le m$ and $1 \le j \le n$:

$$g \le H(i,j) - H(i,j-1) \le M - g \qquad \text{(a)}$$

$$g \le H(i,j) - H(i-1,j) \le M - g \qquad \text{(b)}$$

$$x \le H(i,j) - H(i-1,j-1) \le M \qquad \text{(c)}$$

*Proof.* **Lower bounds.** The recurrence gives $H(i,j) \ge H(i,j-1) + g$, $H(i,j) \ge H(i-1,j) + g$ and $H(i,j) \ge H(i-1,j-1) + s(a_i,b_j) \ge H(i-1,j-1) + x$.

**Upper bounds.** First extend (a) and (b) to the boundary. For $i = 0$, $H(0,j) - H(0,j-1) = 0 \le M - g$, and for $j = 0$, $H(i,0) - H(i-1,0) = 0 \le M - g$. The lower bounds extend in the same way, because $0 \ge g$.

Now use induction on $i + j$. Fix $i, j \ge 1$ and assume the three upper bounds hold at every node with a smaller coordinate sum. Consider which term attains the maximum in the recurrence.

*Case 0: $H(i,j) = 0$.* All three differences are at most 0, because $H \ge 0$ everywhere, and $0 \le M \le M - g$.

*Case 1: $H(i,j) = H(i-1,j-1) + s(a_i,b_j)$.* Then $H(i,j) \le H(i-1,j-1) + M$, which is (c). The extended lower bound (b) at $(i, j-1)$ gives $H(i-1,j-1) \le H(i,j-1) - g$, so $H(i,j) \le H(i,j-1) + M - g$, which is (a). The extended lower bound (a) at $(i-1,j)$ gives $H(i-1,j-1) \le H(i-1,j) - g$, so $H(i,j) \le H(i-1,j) + M - g$, which is (b).

*Case 2: $H(i,j) = H(i-1,j) + g$.* Then (b) holds because $g \le M - g$. By the induction hypothesis (a) at $(i-1,j)$, or its boundary form, $H(i-1,j) \le H(i-1,j-1) + M - g$, so $H(i,j) \le H(i-1,j-1) + M$, which is (c). Combining this with $H(i-1,j-1) \le H(i,j-1) - g$ gives (a).

*Case 3: $H(i,j) = H(i,j-1) + g$.* Symmetric to case 2, with the roles of $i$ and $j$ exchanged. $\square$

**Corollary 3.2.** A side-by-side difference (a) or (b) takes one of $M - 2g + 1$ integer values and needs $\lceil \log_2(M - 2g + 1) \rceil$ bits. A same-diagonal difference (c) takes one of $M - x + 1$ values and needs $\lceil \log_2(M - x + 1) \rceil$ bits.

With the project values these are the ranges $-2 \dots +4$ (3 bits) and $-1 \dots +2$ (2 bits). The simulated reads showed exactly these ranges (`reference/ideas_study.py`, `reference/window_study.py`).

**Caveat.** Proposition 3.1 concerns the full matrix. In a window, a cell at the edge has a neighbour outside the window that is read as 0, and the bound does not apply to that pair. The simulation measured only pairs with both cells inside the window.

## 4. Score as a shortest-path cost

For a path $P$ from $(i_0, j_0)$ to $(i_1, j_1)$ let $\mu$, $\nu$, $\gamma$ be its numbers of matches, mismatches and gap steps. Let

$$u = i_0 + (m - i_1) + j_0 + (n - j_1)$$

be the number of bases of $a$ and $b$ that $P$ leaves unaligned. Define the cost

$$c(P) = (M - x)\,\nu + \left(\tfrac{M}{2} - g\right)\gamma + \tfrac{M}{2}\,u .$$

All three coefficients are non-negative.

**Proposition 4.1.** For every path $P$, $\;\operatorname{score}(P) = \tfrac{M}{2}(m+n) - c(P)$. Therefore

$$S = \tfrac{M}{2}(m + n) - \min_P c(P).$$

*Proof.* A diagonal step advances both coordinates and a gap step advances one, so $(i_1 - i_0) + (j_1 - j_0) = 2(\mu + \nu) + \gamma$, and hence $m + n = 2(\mu+\nu) + \gamma + u$. Substituting,

$$\tfrac{M}{2}(m+n) - c(P) = M\mu + M\nu + \tfrac{M}{2}\gamma + \tfrac{M}{2}u - (M-x)\nu - \tfrac{M}{2}\gamma + g\gamma - \tfrac{M}{2}u = M\mu + x\nu + g\gamma ,$$

which is the score of $P$. The map $P \mapsto c(P)$ is therefore a decreasing function of the score over the same set of paths, the empty path included, so the maximum score corresponds to the minimum cost. $\square$

With the project values the costs are 0 (match), 3 (mismatch), 3 (gap base), 1 (unaligned base), and $S = m + n - c_{\min}$. This is the form a time-coded ("first arrival") circuit would compute. The identity was also checked on 400 random pairs.

## 5. The sliding window

**Definition.** The window is the set of cells $\mathcal{W} = \{(i,j) : 1 \le i \le m,\ 1 \le j \le n,\ \lvert j - i - d_0 \rvert \le w\}$. $H_{\mathcal W}$ is defined by the recurrence of section 2 on cells of $\mathcal W$ and is 0 at every other node. $S_{\mathcal W} = \max H_{\mathcal W}$.

**Proposition 5.1 (work).** $\lvert \mathcal W \rvert \le mW$, so the window computes at most a fraction $W/n$ of the $mn$ cells of the full matrix.

*Proof.* For each $i$ the condition $\lvert j - i - d_0\rvert \le w$ allows at most $2w + 1 = W$ values of $j$. $\square$

For $m = 150$, $n = 180$, $W = 17$ this is 2,550 against 27,000 cells, 9.4%.

Call a path **admissible** if every node on it except possibly its start node is a cell of $\mathcal W$.

**Lemma 5.2.** For a cell $(i,j) \in \mathcal W$, $H_{\mathcal W}(i,j)$ is the maximum score over admissible paths ending at $(i,j)$, the empty path included.

*Proof.* Induction on $i + j$. Let $p$ be a predecessor of $(i,j)$ and $\sigma$ the score of the step from $p$.

If $p \notin \mathcal W$, the recurrence uses $H_{\mathcal W}(p) = 0$ and the candidate is $\sigma$. That is the score of the one-step path from $p$, which is admissible because $p$ is its start node. No longer admissible path passes through $p$, because $p$ would then be a non-start node outside $\mathcal W$.

If $p \in \mathcal W$, the candidate is $H_{\mathcal W}(p) + \sigma$. By the induction hypothesis $H_{\mathcal W}(p)$ is the best score of an admissible path ending at $p$, and appending the step keeps it admissible. Conversely, an admissible path whose last step comes from $p$ is an admissible path ending at $p$ followed by that step.

Taking the maximum over the three predecessors and the empty path gives the recurrence. $\square$

**Proposition 5.3.** $S_{\mathcal W} \le S$.

*Proof.* By Lemmas 2.1 and 5.2, $S_{\mathcal W}$ is a maximum over a subset of the paths whose maximum is $S$. $\square$

**Lemma 5.4.** The diagonals of the nodes of a path form a set of consecutive integers $\{k_{\text{lo}}, \dots, k_{\text{hi}}\}$, and $k_{\text{hi}} - k_{\text{lo}} \le \gamma$, the number of gap steps.

*Proof.* A diagonal step leaves the diagonal unchanged, a right step raises it by 1 and a down step lowers it by 1. A sequence of integers that changes by at most 1 at each step takes every value between its minimum and maximum, and its range is at most the number of steps that change it. $\square$

**Proposition 5.5 (exactness).** Suppose $S > 0$ and some optimal path $P^\ast$ has diagonal range $[k_{\text{lo}}, k_{\text{hi}}]$ with

$$k_{\text{hi}} - w \;\le\; d_0 \;\le\; k_{\text{lo}} + w .$$

Then $S_{\mathcal W} = S$. In particular this holds if $d_0$ is one of the diagonals of $P^\ast$ and $P^\ast$ has at most $w$ gap steps.

*Proof.* An optimal path with positive score does not begin with a gap step, since deleting that step would raise the score by $-g > 0$. So $P^\ast$ begins with a diagonal step, and every node after its start node has $i \ge 1$ and $j \ge 1$, that is, it is a cell. The hypothesis says every diagonal $k$ of $P^\ast$ satisfies $\lvert k - d_0\rvert \le w$, so those cells lie in $\mathcal W$ and $P^\ast$ is admissible. By Lemma 5.2, $S_{\mathcal W} \ge \operatorname{score}(P^\ast) = S$, and Proposition 5.3 gives equality.

For the last sentence, if $d_0 \in [k_{\text{lo}}, k_{\text{hi}}]$ and $\gamma \le w$, Lemma 5.4 gives $k_{\text{hi}} - k_{\text{lo}} \le w$, hence $k_{\text{hi}} - w \le k_{\text{lo}} \le d_0 \le k_{\text{hi}} \le k_{\text{lo}} + w$. $\square$

If $S = 0$ then $S_{\mathcal W} = 0$ as well, by $0 \le S_{\mathcal W} \le S$.

Proposition 5.5 is conditional. Whether the vote produces a $d_0$ that satisfies its hypothesis is the subject of section 6.

## 6. The diagonal vote

In this section indices start at 0. The read is $r_0 \dots r_{m-1}$ and the window is $b_0 \dots b_{n-1}$. Diagonal $k$ pairs $r_i$ with $b_{i+k}$, which is the same diagonal numbering as before.

**Definitions.** For $-(m-1) \le k \le n-1$ let

- $I(k) = \{i : 0 \le i < m,\ 0 \le i + k < n\}$ and $N(k) = \lvert I(k)\rvert$, the number of compared positions. $N(k)$ depends only on $k$, $m$, $n$.
- $V_r(k) = \lvert\{i \in I(k) : r_i = b_{i+k}\}\rvert$, the number of matches.
- $Z_r(k) = V_r(k) - \lambda\,(N(k) - V_r(k)) = (1+\lambda)V_r(k) - \lambda N(k)$ for a fixed penalty $\lambda \ge 0$.

The vote returns a diagonal $d_0$ that maximises $Z_r$. $\lambda = 0$ counts matches. $\lambda = 1$ is matches minus mismatches.

**Read model.** Let $E$ be the reference genome around the window, with $b_j = E[\omega + j]$. A **template** $t_0 \dots t_{m-1}$ is given by a start $p$ and at most one indel of length $g$ at read position $q$:

- no indel: $t_i = E[p + i]$;
- deletion: $t_i = E[p+i]$ for $i < q$ and $t_i = E[p + g + i]$ for $i \ge q$;
- insertion: $t_i = E[p+i]$ for $i < q$, $t_i$ is unconstrained for $q \le i < q + g$, and $t_i = E[p - g + i]$ for $i \ge q + g$.

Let $F$ be the set of unconstrained positions (empty unless there is an insertion). A read $r$ **follows the template with at most $e$ substitutions** if $r_i = t_i$ for all $i \notin F$ except at most $e$ positions.

The **true diagonals** are $k_1 = p - \omega$ for the part before the indel, and $k_2 = k_1$ (no indel), $k_1 + g$ (deletion) or $k_1 - g$ (insertion) for the part after it. Write $\kappa_{\text{lo}} = \min(k_1,k_2)$, $\kappa_{\text{hi}} = \max(k_1,k_2)$. The **true alignment** $T$ pairs each read base with the reference base it was copied from, wherever that base lies in the window. All its diagonals lie in $[\kappa_{\text{lo}}, \kappa_{\text{hi}}]$.

A diagonal $d$ is **acceptable** if the window centred on it contains every true diagonal:

$$\mathcal A = \{d : \kappa_{\text{hi}} - w \le d \le \kappa_{\text{lo}} + w\}.$$

For the template alone define the pessimistic and optimistic counts

$$V^-(k) = \lvert\{i \in I(k)\setminus F : t_i = b_{i+k}\}\rvert, \qquad V^+(k) = V^-(k) + \lvert I(k) \cap F\rvert ,$$

and $Z^\pm(k) = (1+\lambda)V^\pm(k) - \lambda N(k)$. The **margin** of the template is

$$\Delta_\lambda = \max_{k \in \mathcal A} Z^-(k) \;-\; \max_{k \notin \mathcal A} Z^+(k).$$

$\Delta_\lambda$ depends only on $E$, the window, $m$, $w$, $\lambda$ and the template parameters $(p, q, g, \text{type})$. It can be computed before any read is seen.

**Lemma 6.1.** If $r$ follows the template with at most $e$ substitutions, then for every $k$

$$V^-(k) - e \;\le\; V_r(k) \;\le\; V^+(k) + e .$$

*Proof.* Split $I(k)$ into three sets. On positions outside $F$ where $r_i = t_i$, the read and the template give the same indicator. On the at most $e$ substituted positions each indicator changes by at most 1. On positions in $F$ the read contributes between 0 and 1 each, which is between the contribution counted in $V^-$ (none) and in $V^+$ (all). $\square$

**Theorem 6.2 (vote).** If $r$ follows the template with at most $e$ substitutions and

$$2\,(1+\lambda)\,e \;<\; \Delta_\lambda ,$$

then every maximiser of $Z_r$ is acceptable.

*Proof.* Let $k_A \in \mathcal A$ attain $\max_{k\in\mathcal A} Z^-(k)$. By Lemma 6.1, $Z_r(k_A) \ge Z^-(k_A) - (1+\lambda)e$. For any $k \notin \mathcal A$, Lemma 6.1 and the definition of the margin give

$$Z_r(k) \le Z^+(k) + (1+\lambda)e \le Z^-(k_A) - \Delta_\lambda + (1+\lambda)e < Z^-(k_A) - (1+\lambda)e \le Z_r(k_A).$$

So no diagonal outside $\mathcal A$ maximises $Z_r$. $\square$

**Corollary 6.3 (the window contains the true alignment).** Under the hypothesis of Theorem 6.2, the window centred on $d_0$ contains every diagonal of the true alignment $T$, and

$$\operatorname{score}(T) \;\le\; S_{\mathcal W} \;\le\; S .$$

*Proof.* $d_0 \in \mathcal A$ means $d_0 - w \le \kappa_{\text{lo}}$ and $\kappa_{\text{hi}} \le d_0 + w$, so every diagonal of $T$ is within $w$ of $d_0$. Take $T$ restricted to the grid, starting at its first paired base. Its nodes after the start are cells on those diagonals, so it is admissible, and Lemma 5.2 gives $\operatorname{score}(T) \le S_{\mathcal W}$. The second inequality is Proposition 5.3. $\square$

This is the guarantee that matters for calling mutations: the window cannot lose the alignment the read actually came from. It does not say $S_{\mathcal W} = S$; see section 12.

### 6.4 Margins for the three loci (Computed)

`reference/vote_margin.py` evaluates $\Delta_\lambda$ exactly for $m = 150$, $w = 8$, every read start for which the read covers the hotspot centre with at least 20 bp on each side (111 starts per locus), every indel position, every length $g = 1,\dots,8$ and both indel types: 260,739 templates per locus, 782,217 in total. Windows: rpoB 180 bp, katG and inhA 120 bp, from NC_000962.3. The window coordinates are provisional.

Minimum margin, and the largest $e$ with $2(1+\lambda)e < \Delta$:

| Locus | Case | $\lambda = 0$: min $\Delta$ | tolerated $e$ | $\lambda = 1$: min $\Delta$ | tolerated $e$ |
|-------|------|------:|------:|------:|------:|
| rpoB | no indel | 56 | 27 | 106 | 26 |
| rpoB | deletion $\le 8$ | 5 | 2 | 12 | 2 |
| rpoB | insertion $\le 8$ | −4 | none | 10 | 2 |
| katG | no indel | 32 | 15 | 78 | 19 |
| katG | deletion $\le 8$ | 1 | 0 | 9 | 2 |
| katG | insertion $\le 8$ | −9 | none | 4 | 0 |
| inhA | no indel | 36 | 17 | 78 | 19 |
| inhA | deletion $\le 8$ | 5 | 2 | 14 | 3 |
| inhA | insertion $\le 8$ | −4 | none | 9 | 2 |

Share of all templates for which the vote is guaranteed with at least $e$ substitutions:

| $e$ | 0 | 1 | 2 | 3 | 5 | 8 |
|-----|---|---|---|---|---|---|
| $\lambda = 0$ | 99.91% | 99.83% | 99.70% | 99.50% | 98.79% | 96.55% |
| $\lambda = 1$ | 100% | > 99.99% | > 99.99% | 99.95% | 98.99% | 92.61% |

These shares weight every template equally. They are not probabilities over real reads.

What the numbers say:

1. **With the plain match count ($\lambda = 0$) the vote is not guaranteed.** For some insertion templates the margin is negative, so a wrong diagonal can win even with no substitution at all. The cause is that a diagonal on which more of the read overlaps the window accumulates chance matches over more positions. The sequence is GC-rich and about 40% of unrelated positions match.
2. **With $\lambda = 1$ every template has a positive margin.** A read with no substitutions always votes for an acceptable diagonal, for every start and every single indel of up to 8 bp at all three loci.
3. **Reads without an indel are safe by a wide margin** under either statistic: at least 15 substitutions are tolerated.
4. **Indel reads have thin worst-case margins.** A few templates tolerate no substitution or only two. The bound is adversarial (it assumes every substitution removes a match from the true diagonal and adds one to the same wrong diagonal), so it is conservative.

Cost of $\lambda = 1$ in hardware: $Z(k) = 2V(k) - N(k)$, and $N(k)$ is a function of the shift only, so it is one counter and one subtraction.

**Corollary 6.5 (independent substitutions).** Suppose each of the $m$ read bases is substituted independently with probability $\varepsilon$, and the template tolerates $e^\ast$ substitutions. Then

$$\Pr[\text{the vote is not acceptable}] \;\le\; \Pr[\operatorname{Bin}(m,\varepsilon) > e^\ast] = \sum_{i > e^\ast} \binom{m}{i}\varepsilon^i(1-\varepsilon)^{m-i}.$$

*Proof.* By Theorem 6.2 the vote can fail only if more than $e^\ast$ substitutions occur. $\square$

For $m = 150$:

| $e^\ast$ | applies to | $\varepsilon = 0.005$ | $\varepsilon = 0.01$ | $\varepsilon = 0.02$ |
|------|------------|------|------|------|
| 19 | every read without an indel, $\lambda = 1$ | $1.9\times10^{-22}$ | $1.0\times10^{-16}$ | $3.1\times10^{-11}$ |
| 3 | 99.95% of templates, $\lambda = 1$ | $7.1\times10^{-3}$ | $6.5\times10^{-2}$ | $3.5\times10^{-1}$ |

The second row is weak because the bound is adversarial, not because the vote fails that often. Here $\varepsilon$ must include true SNPs as well as sequencing errors.

**Empirical check.** On 16,500 simulated reads (SNPs, 0.5% errors, indels up to 8 bp), the $W = 17$ window disagreed with full SW on 1 read with either statistic, although the two statistics chose different diagonals on 143 reads. So the change to $\lambda = 1$ buys a guarantee; it showed no measurable difference in simulation.

## 7. Storage, cells and cycles (Proved by counting)

| Quantity | Full SW, systolic | Window SW | Diagonal vote, streaming |
|----------|-------------------|-----------|--------------------------|
| DP cells in hardware | $\min(m,n)$ = 150 | $W$ = 17 | none ($m$ two-bit comparators) |
| Cycles per read | $m + n - 1$ = 329 | $m$ = 150 | $m + n - 1$ = 329 |
| Traceback bits | $2mn$ = 54,000 | $2mW$ = 5,100 | none |

- Without pipelining, one read takes $(m + n - 1) + m = 479$ cycles through the vote and the window. That is more than the 329 cycles of a full systolic array.
- With the two stages pipelined, throughput is set by the vote: one read per $m + n - 1$ cycles, the same as the full array.
- For the DP part alone, cells times cycles falls from about $\min(m,n)(m+n)$ to $mW$.

The design therefore reduces area and storage. It does not reduce the number of clock cycles per read.

## 8. Cost of the vote

**Proposition 8.1.** Computing $V(k)$ for every $k$ by direct comparison uses exactly $mn$ base comparisons.

*Proof.* $\sum_k N(k)$ counts the pairs $(i, j)$ with $0 \le i < m$, $0 \le j < n$, each once, on its diagonal $k = j - i$. $\square$

So in a sequential model the direct vote is $\Theta(mn)$, the same order as full SW. Its advantage there is the constant: a two-bit comparison against a DP cell update. With machine words of $\omega$ bits the comparisons pack into $O(mn/\omega)$ word operations.

**Proposition 8.2 (Classical).** All $V(k)$ can be computed in $O((m+n)\log(m+n))$ time.

*Sketch.* $V(k) = \sum_{c \in \{A,C,G,T\}} \sum_i [r_i = c]\,[b_{i+k} = c]$. For each letter the inner sum is a cross-correlation of two 0/1 vectors, computed for all $k$ at once with a fast Fourier transform. Four transforms suffice. This is the match-counting technique of Fischer and Paterson (1974). $\square$

Hence vote plus window is $O((m+n)\log(m+n) + mW)$ in a sequential model, which is below $\Theta(mn)$ when $W = o(n)$. The hardware does not use this method, so this bound is not a claim about the hardware.

## 9. The wavefront (event-driven) alternative

This section concerns Myers's algorithm for two strings of lengths $N_a$ and $N_b$. $D$ is the minimum number of single-base insertions and deletions that turn one into the other; a substitution counts as 2. Put the strings on a grid with coordinates $(x, y)$ and number diagonals by $k = x - y$. A horizontal step raises $x$ and $k$; a vertical step raises $y$ and lowers $k$; a diagonal step along a match is free.

The algorithm works in rounds $d = 0, 1, 2, \dots$. In round $d$ it visits the diagonals $k = -d, -d+2, \dots, d$. For each it takes the starting point

$$x_{\text{start}}(d,k) = \max\{\,x(d-1,k-1) + 1,\; x(d-1,k+1)\,\}$$

and then slides along matches as far as possible, giving $x(d,k)$. It stops when a slide reaches $(N_a, N_b)$.

**Lemma 9.1 (Classical, Myers 1986).** $x(d,k)$ is the largest $x$ on diagonal $k$ reachable from the origin with exactly $d$ non-diagonal steps.

*Sketch.* Any $d$-step path to diagonal $k$ makes its last non-diagonal step from diagonal $k-1$ or $k+1$ after $d-1$ steps. Replacing its prefix by the furthest-reaching $(d-1)$-path on that diagonal does not move the step's landing point backwards, and sliding along matches afterwards only moves it forwards. Induct on $d$. $\square$

It follows that the algorithm stops in round $D$.

**Proposition 9.2 (updates).** The number $U$ of (round, diagonal) updates satisfies

$$\frac{D(D+1)}{2} + 1 \;\le\; U \;\le\; \frac{(D+1)(D+2)}{2}.$$

*Proof.* Round $d$ has $d + 1$ diagonals. Rounds $0, \dots, D-1$ run in full, $\sum_{d=0}^{D-1}(d+1) = D(D+1)/2$ updates. Round $D$ makes at least 1 and at most $D+1$. $\square$

**Proposition 9.3 (comparisons).** The number $C$ of base comparisons satisfies $C \le U + (2D+1)\min(N_a, N_b)$.

*Proof.* Each slide ends at its first mismatch or at the edge of the grid, so there is at most one failed comparison per update: at most $U$ in total.

A successful comparison advances the slide by one point on its diagonal. We show no point is passed twice. Let $(x, y)$ be where the slide on diagonal $k$ stopped in round $d - 2$. A vertical step from it lands on diagonal $k - 1$ with the same $x$, so $x(d-1,k-1) \ge x$. A horizontal step lands on diagonal $k+1$ at $x + 1$, so $x(d-1,k+1) \ge x+1$. If a step would leave the grid, the other one still gives the bound, or the start falls outside the grid and the slide makes no comparison. Therefore $x_{\text{start}}(d,k) \ge x + 1$: each slide begins beyond where the previous slide on the same diagonal stopped.

At most $2D + 1$ diagonals are ever visited, and a diagonal has at most $\min(N_a, N_b)$ points. $\square$

**Corollary 9.4.** The worst-case work is $O(D\min(N_a,N_b) + D^2)$.

Both bounds were checked on 3,000 random pairs with no violation. For highly similar sequences the observed work is close to $N + D^2$, because almost all matches lie on one diagonal. That is an empirical observation with an expected-case analysis in the literature. It is not a worst-case bound.

## 10. The segment-vote ("attention") band

Split the read into consecutive segments of length at most $S_{\text{seg}}$. Segment $t$ votes for a diagonal $d_t$. Row $i$ in segment $t$ computes the diagonals from $\min(d_{t-1},d_t,d_{t+1}) - w$ to $\max(d_{t-1},d_t,d_{t+1}) + w$, where a missing neighbour is omitted. Let $G = \sum_t \lvert d_{t+1} - d_t \rvert$.

**Proposition 10.1.** The number of cells computed is at most $mW + 2S_{\text{seg}}G$.

*Proof.* A row in segment $t$ computes at most $R_t + W$ cells, where $R_t$ is the range of $\{d_{t-1}, d_t, d_{t+1}\}$. The range of three numbers is at most the sum of the two consecutive absolute differences, $R_t \le \lvert d_t - d_{t-1}\rvert + \lvert d_{t+1} - d_t\rvert$. Summing over the at most $S_{\text{seg}}$ rows of each segment and over all segments,

$$\text{cells} \le mW + S_{\text{seg}}\sum_t R_t \le mW + 2S_{\text{seg}}G ,$$

since each difference $\lvert d_{t+1} - d_t\rvert$ appears in at most two of the $R_t$. $\square$

With no indel $G = 0$ and the count is at most $mW$. With one indel of length $g$, $G = g$. In the worst case the band covers the whole matrix.

This bounds the work only. The segment-vote band missed 7-13% of 5-12 bp indels in simulation, all within one segment length of the end of the overlap, and there is no correctness theorem for it.

## 11. The gate

**Proposition 11.1.** If a read passes the gate with probability $\rho$, the vote costs $C_{\text{vote}}$ and the window costs $C_{\text{SW}}$, then the expected cost per read is

$$\mathbb E[T] = C_{\text{vote}} + (1 - \rho)\,C_{\text{SW}}, \qquad C_{\text{vote}} \le \mathbb E[T] \le C_{\text{vote}} + C_{\text{SW}} .$$

*Proof.* Every read pays for the vote. The window is run with probability $1 - \rho$. Linearity of expectation. $\square$

$\rho$ is an empirical quantity. It depends on the threshold, which is not yet chosen, and on the reads.

## 12. What is not proved

1. **$S_{\mathcal W} = S$ in general.** Corollary 6.3 gives $\operatorname{score}(T) \le S_{\mathcal W} \le S$. Full SW may find a higher-scoring alignment that leaves the window, and nothing here excludes that. In simulation it happened in about 1 read in 16,500.
2. **The vote for reads outside the model.** Theorem 6.2 covers one indel of length at most $w$ and substitutions. Two or more indels, chimeric reads, and reads from elsewhere in the genome are not covered.
3. **The vote with $\lambda = 0$.** The margin is negative for some templates, so there is no guarantee.
4. **The margins on the final windows.** Section 6.4 used provisional window coordinates and forward-strand reads. The numbers must be recomputed when the Bio focus fixes the coordinates and the strand decision is made.
5. **The gate.** There is no theorem about which reads may skip SW. In simulation 3.6% of indel reads had an ungapped score equal to the full SW score because the indel was within 14 bp of the end of the overlap.
6. **The segment-vote band and the adaptive band.** Only the work bound of Proposition 10.1. The adaptive band was worse than a fixed band in simulation.
7. **Typical-case cost of the wavefront method.** Section 9 proves worst-case bounds only.
8. **Anything about speed or power of real hardware.** Section 7 counts cells and cycles. Clock frequency, routing and power have not been measured.

## References

- Smith, T. F., & Waterman, M. S. (1981). Identification of common molecular subsequences. *Journal of Molecular Biology, 147*(1), 195-197.
- Myers, E. W. (1986). An O(ND) difference algorithm and its variations. *Algorithmica, 1*, 251-266.
- Fischer, M. J., & Paterson, M. S. (1974). String-matching and other products. In *Complexity of Computation*, SIAM-AMS Proceedings 7, 113-125.

## Scripts

| Result | Script |
|--------|--------|
| Margins of section 6.4 | `reference/vote_margin.py` |
| Difference ranges of section 3 | `reference/ideas_study.py`, `reference/window_study.py` |
| Window and vote on simulated reads | `reference/band_study.py` |
| Wavefront counts of section 9 | `reference/ideas_study.py` |
