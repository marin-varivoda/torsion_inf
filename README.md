# torsion_inf

To ensure this repository is as self-contained as possbile, we include copies of the necessary data (models for $X1(N)$, $X1(M, N)$) and third-party libraries in the `third-party` directory. This keeps the code functional even if the external sources become unavailable.

## Data Sources and Attribution

### Optimized equations for $X_1(M, N)$
The curve models in the `third_party/DS17_X1MN_models` directory were obtained from the [Optimized equations for X_1(m,mn)](https://math.mit.edu/~drew/X1mn.html) database. 

We thank the authors Maarten Derickx and Andrew V. Sutherland for constructing these models. The method for obtaining them is detailed in their paper:

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