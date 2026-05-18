AttachSpec("../third_party/CurveArith/CurveArith.spec");
AttachSpec("../third_party/mdmagma/v2/mdmagma.spec");

N := 48;

C := MDX1(N, Rationals());
C_proj := Curve(C);

// Reduce the curve over finite field
p := 5;
C_fin := Curve(Reduction(C_proj, p));

printf "Initialized the model for X_1(%o) and reduced it over a finite field with %o elements.\n", N, p;
printf "Checking if reduced X_1(%o) has a function of degree at most 9... (this might take a moment)\n", N;

T := Time();
has_deg_9_function := HasFunctionOfDegreeAtMost(C_fin, 9); // Function from CurveArith by Derickx-Terao

printf "Calculation complete. Has a function of degree at most 9: %o\n", has_deg_9_function;
printf "Total calculation time: %o seconds\n", Time(T);

exit;
