AttachSpec("../third_party/CurveArith/CurveArith.spec");
load "../third_party/DS17_X1MN_models/X1_5_10.txt";

// Construct a curve from model
A<u,v> := AffineSpace(F, 2);
C := Curve(A, X);
C_proj := ProjectiveClosure(C);

// Base change to Q(zeta_5)
K<z> := CyclotomicField(5);
f_to_k := hom< F -> K | z^2 + z^3 >; // Map 'a' (root of x^2 - x - 1) to (z^2 + z^3) in Q(zeta_5)
assert a^2 + a - 1 eq 0;
assert f_to_k(a)^2 + f_to_k(a) - 1 eq 0; // Verify that (z^2+z^3) is the root of the same polynomial as a

C_proj_K := BaseChange(C_proj, f_to_k);

// Find a prime ideal lying over p in the field of definition Q(zeta_5)
p := 11;
O := MaximalOrder(K);
prime_ideal := Decomposition(O, p)[1,1];
res_field, res_map := ResidueClassField(prime_ideal);

// Reduce the curve over finite field
C_fin := Curve(Reduction(C_proj_K, prime_ideal));

printf "Initialized the model for X_1(%o, %o) and reduced it over a finite field with %o elements.\n", m, m*n, #res_field;
printf "Calculating gonality... (this might take a moment)\n";

T := Time();
gonality := Gonality(C_fin);

printf "Gonality for the curve X_1(%o, %o) reduced over a finite field of size %o is: %o\n", m, m*n, #res_field, gonality;
printf "Total calculation time: %o seconds\n", Time(T);

exit;
