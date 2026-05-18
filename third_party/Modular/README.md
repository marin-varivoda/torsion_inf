# Modular

This `Magma` code is for computing modular forms and models of modular curves.  It is an improvement of some [previous code](https://github.com/davidzywina/OpenImage) that was used to obtain a computational version of Serre's open image theorem.  

Details on the theory behind this code can be found in the preprint [Classification of modular curves with low gonality](https://pi.math.cornell.edu/~zywina/papers/lowgonality.pdf).  The computations in that paper can be found in the directory `classification`.

To load the code, type the following in `Magma` from the base directory of the repository (or with appropriate pathname otherwise):

        AttachSpec("Modular.spec");
        
Given a subgroup `G` of `GL(2,Z/NZ)`, the code

        X:=CreateModularCurveRec(G);
creates the corresponding modular curve `X=X_G`.  You could also use
        
        X:=CreateModularCurveRec(N,gens);
where `gens` is a sequence of generators of `G`.

Much of the computations will be defined over a number field ``K:=X`KG``.  When `G` has full determinant, the field `K` is just the rationals.   The code is designed to work well when the degree `[K:Q]` is small.  Note that one has to be a little careful concerning whether bases are over `Q` or `K`.  We will now suppose that `K=Q`.

For an integer weight `k>1`, the following computes a basis of the space `M_{k,G}` of modular forms from the paper.

        X:=FindModularForms(k,X);
A basis for the modular forms is given by ``X`F``.   Note that each modular form is given as a sequence consisting of the q-expansion at each cusp. There is a slight difference with the previous version: the q-expansion is now given in terms of ``qw:=q^(1/w)``, where `w` is the width of the cusp being considered.

Once the modular forms have been computed, the following computes the cusp forms:

        X:=FindCuspForms(X);
with the basis given by ``X`F0``. There is also a function `FindModularFormsWithVanishingConditions` where one can impose other vanishing conditions at the cusps.

When `X` has genus at least 3, we can compute the image of the canonical map:

        X:=FindCanonicalModel(X);
the equations are given in ``X`psi``. When the curve is geometrically hyperelliptic, the equations will cut out a genus 0 curve.

The following computes a model of the curve:

        X:=FindModelOfXG(X);
It will use the canonical model when possible.  Otherwise it will use what the LMFDB calls an *embedded model*.
The equations for the model are found in ``X`psi``.   The polynomials in ``X`psi`` will all be homogenous of degree 2 or 3 (except in the case where `X` is not geometrically hyperelliptic and has genus 3; it is then given as a smooth plane quartic).

A significant change over the previous version is how the precision of the q-expansions is handled.  For example, none of the above functions require precision as input!  
The following function increases the precision of your modular forms

        X:=IncreaseModularFormPrecision(X,prec),
where `prec` indicates the desired precision.   This will replace the modular forms in ``X`F`` and ``X`F0``, if they exist, with the same modular forms but with appropriately increased precision.  


