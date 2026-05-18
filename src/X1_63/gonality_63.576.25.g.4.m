AttachSpec("../third_party/CurveArith/CurveArith.spec");
AttachSpec("../third_party/Modular/Modular.spec");
load "X1_63/63.576.25.g.4.m";

// Compute the F2-gonality of the intermediate modular curve with 
// LMFDB label 63.576.25.g.4.
//
// This is an intermediate curve in the tower
//     X_1(63) -> 63.576.25.g.4 -> X_1(21).
//
// The LMFDB data states that the Jacobian of this curve has analytic rank 0,
// which implies the algebraic rank is also 0, so gonality computations are
// very useful since they give a tight bound on density degree. In fact,
// a lower bound gon_2(63.576.25.g.4) >= 10 implies that only finitely many points 
// of degree <= 9 exist on this intermediate curve over Q. Since X_1(63) has 
// a degree 3 map to it, the same  finiteness conclusion follows for
// degree <= 9 points on X_1(63).

lmfdb_label := "63.576.25.g.4";

// Use Zywina's Modular package to construct a model from the LMFDB subgroup
// generators.
X := CreateModularCurveRec(level, gens);
X := FindModelOfXG(X);

R := Universe(X`psi); // Extract the polynomial ring where equations live
C := Curve(ProjectiveSpace(R), X`psi); // Construct the modular curve

// Reduce at a small prime not dividing the level
p := 2;
C_fin := Curve(Reduction(C, p));

assert IsNonsingular(C_fin);
assert Genus(C_fin) eq g;

printf "Initialized modular curve %o and reduced it over a finite field with %o elements.\n", lmfdb_label, p;
printf "Genus of the reduction: %o\n", Genus(C_fin);
printf "Computing the gonality of the reduction... (this might take a moment)\n";

T := Time();
gonality := Gonality(C_fin);

printf "Calculation complete. Gonality of the reduction: %o\n", gonality;

exit;