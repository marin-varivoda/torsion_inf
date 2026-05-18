AttachSpec("../third_party/CurveArith/CurveArith.spec");
load "../third_party/DS17_X1MN_models/X1_3_18.txt";

// Construct a curve from model
A<u,v> := AffineSpace(F, 2);
C := Curve(A, X);
C_proj := ProjectiveClosure(C);

// Find a prime ideal lying over p in the field of definition F = Q(zeta_3)
p := 7;
O := MaximalOrder(F);
prime_ideal := Decomposition(O, p)[1,1];
res_field, res_map := ResidueClassField(prime_ideal);

// Reduce the curve over finite field
C_fin := Curve(Reduction(C_proj, prime_ideal));


printf "Initialized the model for X_1(%o, %o) and reduced it over a finite field with %o elements.\n", m, m*n, #res_field;
printf "Calculating gonality... (this might take a moment)\n";

T := Time();
gonality := Gonality(C_fin);

printf "Gonality for the curve X_1(%o, %o) reduced over a finite field of size %o is: %o\n", m, m*n, #res_field, gonality;
printf "Total calculation time: %o seconds\n", Time(T);
