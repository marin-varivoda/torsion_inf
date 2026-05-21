AttachSpec("/home/mvarivoda/mdmagma/v2/mdmagma.spec");
SetSeed(1337);


N := 32;

C := MDX1(N, Rationals());
C_proj := Curve(C);

// Reduce the curve over finite field
p := 5;
C_fin := Curve(Reduction(C_proj, p));

print "Starting class group computation";
cg, cgToDiv, divToCg := ClassGroup(C_fin);

tcg := TorsionSubgroup(cg);

print "Torsion subgroup size:", #tcg;
