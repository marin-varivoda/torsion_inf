QQ := Rationals();
A<x, y> := AffineSpace(QQ, 2);

// Sutherland's FFc32.txt model from https://math.mit.edu/~drew/X1_altcurves.html
f := x^11*y^5 - 2*x^11*y^4 + 2*x^11*y^3 - 5*x^10*y^6 + 7*x^10*y^5 - 6*x^10*y^4 - 
    2*x^10*y^3 + 10*x^9*y^7 - 5*x^9*y^6 + 3*x^9*y^5 + 5*x^9*y^4 + 6*x^9*y^3 - 
    4*x^9*y^2 - 10*x^8*y^8 - 10*x^8*y^7 + 5*x^8*y^6 - x^8*y^5 - 14*x^8*y^4 + 
    6*x^8*y^3 + 4*x^8*y^2 + 5*x^7*y^9 + 20*x^7*y^8 - 5*x^7*y^6 + 5*x^7*y^4 - 
    7*x^7*y^3 - 6*x^7*y^2 + 3*x^7*y - x^6*y^10 - 13*x^6*y^9 - 12*x^6*y^8 + 
    20*x^6*y^6 - 4*x^6*y^5 + 3*x^6*y^4 + 3*x^6*y^3 + x^6*y^2 - 3*x^6*y + 
    3*x^5*y^10 + 11*x^5*y^9 + 7*x^5*y^8 - 10*x^5*y^7 - 15*x^5*y^6 - 13*x^5*y^5 +
    29*x^5*y^4 - 20*x^5*y^3 + 9*x^5*y^2 + x^5*y - x^5 - 3*x^4*y^10 - 5*x^4*y^9 -
    6*x^4*y^8 + 14*x^4*y^7 + 25*x^4*y^6 - 33*x^4*y^5 + 16*x^4*y^4 - 12*x^4*y^3 +
    3*x^4*y^2 + x^4*y + x^3*y^10 + 4*x^3*y^9 - 12*x^3*y^7 - 2*x^3*y^6 + 
    7*x^3*y^5 + 8*x^3*y^4 - 4*x^3*y^3 - 3*x^3*y^2 + x^3*y - 2*x^2*y^9 + 
    8*x^2*y^7 - 7*x^2*y^6 + 9*x^2*y^5 - 24*x^2*y^4 + 24*x^2*y^3 - 9*x^2*y^2 + 
    x^2*y + x*y^8 - 11*x*y^6 + 25*x*y^5 - 25*x*y^4 + 14*x*y^3 - 5*x*y^2 + x*y - 
    y^6 + 5*y^5 - 10*y^4 + 10*y^3 - 5*y^2 + y;

C := ProjectiveClosure(Curve(A, f));
F := FunctionField(C);

printf "Successfully loaded curve.\n";
printf "Loading cusps from TwoGenerators representation...\n";

// Coerce x and y to F before loading qcusps_twogen32.m
x := F ! x;
y := F ! y;

load "../src/X1_32/qcusps_twogen_32.m";
cusps_from_twogen := [Place(twogen_pair) : twogen_pair in qcusps_twogen];
assert #cusps_from_twogen eq 9;

printf "Computing cusps from scratch...\n";

// Compute cusps from scratch for comparison
r := (x^2*y - x*y + y - 1) / (x^2*y-x);
s := (x*y - y + 1) / (x*y);

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

// TODO elaborate on DvH base
f2 := b/(16*b^2+(1-20*c-8*c^2)*b+c*(c-1)^3);
f3 := b;

function cuspSignatureF23(cusp)
    return [Valuation(f2, cusp), Valuation(f3, cusp)];
end function;

cuspSignatures := [cuspSignatureF23(c) : c in cusps_from_twogen];

print cuspSignatures;


Ebc_disc := -b2^2 * b8 - 8*b4^3 -27*b6^2 + 9 * b2 * b4 * b6;

cusps := Support(Divisor(Ebc_disc));
assert #cusps eq 9; // X1(32) has 9 cusp places

// Check that cusps eq cusps_from_twogen
assert Set(cusps) eq Set(cusps_from_twogen);

printf "Loaded and computed cusps match. Test passed successfully!";
