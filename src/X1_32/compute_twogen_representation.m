QQ := Rationals();
A<x, y> := AffineSpace(QQ, 2);

// Sutherland's FFc32.txt model from https://math.mit.edu/~drew/X1_altcurves.html
// f := x^11*y^5 - 2*x^11*y^4 + 2*x^11*y^3 - 5*x^10*y^6 + 7*x^10*y^5 - 6*x^10*y^4 - 
//     2*x^10*y^3 + 10*x^9*y^7 - 5*x^9*y^6 + 3*x^9*y^5 + 5*x^9*y^4 + 6*x^9*y^3 - 
//     4*x^9*y^2 - 10*x^8*y^8 - 10*x^8*y^7 + 5*x^8*y^6 - x^8*y^5 - 14*x^8*y^4 + 
//     6*x^8*y^3 + 4*x^8*y^2 + 5*x^7*y^9 + 20*x^7*y^8 - 5*x^7*y^6 + 5*x^7*y^4 - 
//     7*x^7*y^3 - 6*x^7*y^2 + 3*x^7*y - x^6*y^10 - 13*x^6*y^9 - 12*x^6*y^8 + 
//     20*x^6*y^6 - 4*x^6*y^5 + 3*x^6*y^4 + 3*x^6*y^3 + x^6*y^2 - 3*x^6*y + 
//     3*x^5*y^10 + 11*x^5*y^9 + 7*x^5*y^8 - 10*x^5*y^7 - 15*x^5*y^6 - 13*x^5*y^5 +
//     29*x^5*y^4 - 20*x^5*y^3 + 9*x^5*y^2 + x^5*y - x^5 - 3*x^4*y^10 - 5*x^4*y^9 -
//     6*x^4*y^8 + 14*x^4*y^7 + 25*x^4*y^6 - 33*x^4*y^5 + 16*x^4*y^4 - 12*x^4*y^3 +
//     3*x^4*y^2 + x^4*y + x^3*y^10 + 4*x^3*y^9 - 12*x^3*y^7 - 2*x^3*y^6 + 
//     7*x^3*y^5 + 8*x^3*y^4 - 4*x^3*y^3 - 3*x^3*y^2 + x^3*y - 2*x^2*y^9 + 
//     8*x^2*y^7 - 7*x^2*y^6 + 9*x^2*y^5 - 24*x^2*y^4 + 24*x^2*y^3 - 9*x^2*y^2 + 
//     x^2*y + x*y^8 - 11*x*y^6 + 25*x*y^5 - 25*x*y^4 + 14*x*y^3 - 5*x*y^2 + x*y - 
//     y^6 + 5*y^5 - 10*y^4 + 10*y^3 - 5*y^2 + y;

f := x^4*y - x^3*y^3 - x^3*y + x^2*y^4 + x^2*y - x^2 - x*y^4 + x*y^3 - x*y^2 + x*y + 
    y^3 - 2*y^2 + y;

C := ProjectiveClosure(Curve(A, f));
F := FunctionField(C);

r := F ! ((x^2*y - x*y + y - 1) / (x^2*y-x));
s := F ! ((x*y - y + 1) / (x*y));

b := r * s * (r-1);
c := s * (r-1);

a1 := 1-c;
a2 := -b;
a3 := -b;
a4 := 0;
a6 := 0;

b2 := a1^2 + 4 * a2;
b4 := 2*a4 + a1*a3;
b6 := a3^2 + 4*a6;
b8 := a1^2 * a6  + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2;

Ebc_disc := -b2^2 * b8 - 8*b4^3 -27*b6^2 + 9 * b2 * b4 * b6;

printf "Successfully loaded curve. Computing cusps...\n";

cusps := Support(Divisor(Ebc_disc));

printf "Found %o cusps.\n", #cusps;
printf "Computing TwoGenerators representation for each cusp...\n";

twogen := [];
for cusp in cusps do
    g1, g2 := TwoGenerators(cusp);
    Append(~twogen, [g1, g2]);
end for;

printf "Finished computing TwoGenerators representations.\n";
printf "Result:\n";

print twogen;
