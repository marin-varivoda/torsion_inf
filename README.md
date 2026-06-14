# torsion_inf

Magma code in this repository is meant to be run from `src/` working directory. Paths for `load` command in Magma are written relative to that path. 

The `src/DT/` directory holds code that uses [CurveArith](https://github.com/nt-lib/CurveArith) (Derickx–Terao) to compute lower bounds on the gonality of modular curves $X_1(m, mn)$. To run all of these in parallel (20 jobs in total) and save logs under `out/DT/`, use:

```sh
./scripts/run_DT.sh
```

### Rough runtime totals

All computations were performed on the **Mordell** workstation at the Department of Mathematics, University of Zagreb, with an **AMD EPYC 9175F 16-Core Processor** and **384 GB of RAM**.

Parallel batches like DT computation use the slowest job. Multi-worker Magma jobs use wall clock from the final `elapsed=` line (not Magma's `Total time` since it tracks cpu time). The $X_1(37)$ run is logged under `src/X1_37/`.

| Category | Time |
|---|---:|
| Cover search ($X_1(53)$ + $X_1(57)$ + $X_1(43)$) | ~377,707 s (~105 h) |
| $X_1(32)$ Magma search | ~16,355 s (~4.5 h) |
| $X_1(37)$ full run | ~7,097 s (~2 h) |
| DT (20 parallel jobs; slowest: $X_1(2,30)$ lb) | ~14,046 s (~3.9 h) |
| Everything else in `out/` | ~240 s |

## Data Sources and Attribution

To ensure this repository is as self-contained as possible, we include copies of the necessary data (models for $X_1(N)$, $X_1(M, N)$) and third-party libraries in the `third-party` directory. This keeps the code functional even if the external sources become unavailable.

### Optimized equations for $X_1(M, N)$
The curve models in the `third_party/DS17_X1MN_models` directory were obtained from the [Optimized equations for X_1(m,mn)](https://math.mit.edu/~drew/X1mn.html) database. 

The method for obtaining them is detailed in their paper:

> **Maarten Derickx and Andrew V. Sutherland**, *Torsion subgroups of elliptic curves over quintic and sextic number fields*, Proceedings of the AMS, **145** (2017), 4233-4245.

---

### Optimized equations for $X_1(N)$
The curve models in the `third_party/SvH_X1N_models` directory were obtained from the [Optimized equations for X_1(N)](https://math.mit.edu/~drew/X1_optcurves.html) database.

These optimized equations $f(x,y)=0$ for $X_1(N)$ are joint work by Andrew V. Sutherland and Mark van Hoeij.

---

## Third-Party Libraries

### mdmagma
- **Authors:** Maarten Derickx and Andrew V. Sutherland
- **Source Repository:** [nt-lib/mdmagma](https://github.com/nt-lib/mdmagma)


### CurveArith
- **Authors:** Maarten Derickx and Kenji Terao
- **Source Repository:** [nt-lib/CurveArith](https://github.com/nt-lib/CurveArith)

> **Maarten Derickx and Kenji Terao**, *Computing class groups and gonalities of algebraic curves over finite fields*, ([arXiv 2026](https://arxiv.org/abs/2602.17417))

### Modular
- **Author:** David Zywina
- **Source Repository:** [davidzywina/Modular](https://github.com/davidzywina/Modular)

> **David Zywina**, *Classification of modular curves with low gonality*, ([preprint](https://pi.math.cornell.edu/~zywina/papers/lowgonality.pdf)).


### mdsage
- **Author:** Maarten Derickx
- **Source Repository:** [koffie/mdsage](https://github.com/koffie/mdsage)
- **Documentation:** [mdsage documentation](https://koffie.github.io/mdsage/doc/html/)

`mdsage` is a SageMath package for computations with modular forms, modular symbols, and modular curves. It builds on top of SageMath and includes additional functionality not yet available in Sage.
