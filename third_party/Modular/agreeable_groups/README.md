# Agreeable groups of higher genus

For a non-CM elliptic curve E over $\mathbb{Q}$, the image of the adelic Galois representation of E is an open subgroup $G_E$ of $GL_2(\widehat{\mathbb{Z}})$.  It is [conjectured ](https://arxiv.org/abs/2206.14959) that the index $[GL_2(\widehat{\mathbb{Z}}):G_E]$ must be one of the following:  2, 4, 6, 8, 10, 12, 16, 20, 24, 30, 32, 36, 40, 48, 54, 60, 72, 80, 84, 96, 108, 112, 120, 128, 144, 160, 182, 192, 200, 216, 220, 224, 240, 288, 300, 
336, 360, 384, 480, 504, 576, 768, 864, 1152, 1200, 1296, 1536, 2736.

If the conjecture fails, then we obtain an "unexpected" rational point on a modular curve $X_G$ of genus at least 2 for a certain open group $G$.

We can take $G$ so that it satisfies one of the following:


1. $G$ has prime power level; this includes the groups in [RSZB](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/ell-adic-images-of-galois-for-elliptic-curves-over-mathbb-q-and-an-appendix-with-john-voight/D5BC92F9949B387570A7D764635B6AC8) as well as normalizers of nonsplit Cartan subgroups modulo $\ell$ for primes $\ell>17$.
    
2. $G$ is one of the finite number of groups given in the file `groups.m`.

This list was made so that one could start determining the rational points of these modular curves $X_G$.

-----------------------------


This follows from a previous computation  [project](https://github.com/davidzywina/OpenImage).   The file `agreeable.dat` contains the information from this early work.


The groups are produced by the file `find_groups.m` (it took around 75 minutes on my machine; your times may vary).  We also look at a few local obstructions to remove some groups $G$ for which $X_G(\mathbb{Q})$ is empty.  

We make use of our faster modular curve code and the significantly faster group theory code of Drew Sutherland.
