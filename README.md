# Erdős 128: sparse halves of triangle-free graphs

**Author name: komeiji shiki. Research draft for review.**

This repository publishes the manuscript and reproducibility snapshot of
5 September 2026 for a candidate general upper bound

$$b(G) \le \frac{2543}{100000}n^2 = 0.02543n^2,$$

where $G$ is a finite triangle-free graph and $b(G)$ is the least number
of edges induced by a set of $\lfloor n/2\rfloor$ vertices.
**The conjectured bound $n^2/50=0.02n^2$ remains unresolved.**

- [Read the manuscript (PDF)](output/pdf/sparse-halves.pdf)
- [LaTeX source](paper/sparse-halves.tex)
- [Download the complete review snapshot (ZIP)](output/erdos128-02543-review.zip)
- [AI-assisted proof audit (Chinese)](notes/proof-audit-2026-09-05.md)
- [Numerical verification report](results/paper-verification.json)

## Status and attribution

AI assistance was substantial, including mathematical exploration,
programming, drafting and a subsequent audit. The author has confirmed
review of the mathematical arguments, the applicability of cited results
and the computational checks. This is an author self-declaration. The
manuscript has not been accepted by a journal, and research priority has
not been established. See [the author-status declaration](results/author-review-status.json).
Computer-generated reports do not establish whether a human performed a
review; their human-review flags describe the scope of those checks.

The argument combines a degree-sensitive averaging criterion, sparse-set
extension and neighbourhoods anchored across a maximum cut. It uses the
max-cut theorems of Balogh, Clemen and Lidický and earlier work by Razborov.
The four-cycle certificate and the original neighbourhood-perturbation
construction are due to Amir Sarid; see the manuscript and
[certificate provenance](certificates/README.md).

The PDF uses the official E-JC 12-point A4 layout, theorem styles and MSC
05C35 classification. The official style file is included without changes.
Submission-only metadata omits unassigned publication dates, volume and DOI.
The author declaration and public materials URL are included in the paper.
The mathematical argument and proof inputs are unchanged.

## Reproduce the numerical checks

Python 3.11 or later is sufficient. No third-party package is needed for:

```text
python -S src/verify_paper.py
python -S -m unittest discover -s tests -v
```

The replay checks five exact rational interval trees with two separately
implemented checkers, covering 8,565 nodes and 4,285 leaves. The four-cycle
checker enumerates all 32,768 labelled six-vertex graphs and verifies the
certificate on the 38 triangle-free isomorphism classes with 34 positive
rational square terms. The three tests included here exercise current
certificates and failure controls. The full local research workspace's
17-test result is recorded separately; its earlier finite-graph experiments
are not assumptions of this manuscript.

Optional symbolic checking of 12 polynomial identities uses SymPy 1.14.0:

```text
python -m pip install sympy==1.14.0
python src/audit_algebra.py
```

The source of the checks and their exact scope are documented in the paper.
The quoted max-cut and large-independence theorems are external inputs.

## Snapshot integrity

`manifest.json` records SHA-256 hashes of the current repository payload.
The ZIP contains its own manifest for the corresponding review archive. Running
checks rewrites result reports, including timings, so those output hashes
may change after a successful replay without changing the proof inputs.
The `.gitattributes` file preserves the byte sequences during Git checkout.

## Sources

- [Erdős problem 128](https://www.erdosproblems.com/128)
- [Balogh, Clemen and Lidický, Max Cuts in Triangle-free Graphs](https://arxiv.org/abs/2103.14179v1)
- [Razborov, More about sparse halves in triangle-free graphs](https://doi.org/10.1070/SM9615)
- [Sarid's pinned source and certificate](https://github.com/aimir/erdos-128-sparse-halves/tree/ba85f319ebb33ec2194ef4c2929a12916329e1e7)

中文说明：这是署名 komeiji shiki 的 AI 辅助研究草稿及计算复核材料。
当前候选上界为 0.02543n²，原题尚未解决。作者已声明完成核验；稿件已按 E-JC 官方模板排版，尚未经期刊同行评审。
