AttachSpec("../third_party/CurveArith/CurveArith.spec");
load "../third_party/DS17_X1MN_models/X1_2_36.txt";

// Construct a curve from model
A<u,v> := AffineSpace(F, 2);
C := Curve(A, X);
C_proj := ProjectiveClosure(C);

// Reduce the curve over finite field
p := 5;
C_fin := Curve(Reduction(C_proj, p));

printf "Initialized the model for X_1(%o, %o) and reduced it over a finite field with %o elements.\n", m, m*n, p;
printf "Checking if reduced X_1(%o, %o) has a function of degree at most 9... (this might take a moment)\n", m, m*n;

T := Time();
has_deg_9_function := HasFunctionOfDegreeAtMost(C_fin, 9); // Function from CurveArith by Derickx-Terao

printf "Calculation complete.  Has a function of degree at most 9: %o\n", has_deg_9_function;
printf "Total calculation time: %o seconds\n", Time(T);

exit;
