AttachSpec("../third_party/CurveArith/CurveArith.spec");
load "../third_party/DS17_X1MN_models/X1_2_24.txt";

// Construct a curve from model
A<u,v> := AffineSpace(F, 2);
C := Curve(A, X);
C_proj := ProjectiveClosure(C);

// Reduce the curve over finite field
p := 7; // Note that gon_5 = 7, while gon_7 = 8
C_fin := Curve(Reduction(C_proj, p));

printf "Initialized the model for X_1(%o, %o) and reduced it over a finite field with %o elements.\n", m, m*n, p;
printf "Computing the exact gonality of X_1(%o, %o)... (this might take a moment)\n", m, m*n;

T := Time();
gonality := Gonality(C_fin);

printf "Calculation complete. Gonality: %o\n", gonality; // 8
printf "Total calculation time: %o seconds\n", Time(T);

exit;