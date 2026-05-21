p := 3;
A<x, y> := AffineSpace(GF(p), 2);

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


// Coerce x and y to F before loading qcusps_twogen32.m
x := F ! x;
y := F ! y;
load "X1_32/qcusps_twogen_32.m";
qcusp_num := #qcusps_twogen;
printf "Successfully loaded twogen representation for %o rational cusp places.\n", #qcusps_twogen;


// Compute cusps over F3
r := (x^2*y - x*y + y - 1) / (x^2*y-x);
s := (x*y - y + 1) / (x*y);

b := r * s * (r-1);
c := s * (r-1);

a1 := 1-c;
a2 := -b;
a3 := -b;
a4 := 0;
a6 := 0;
C := ProjectiveClosure(Curve(A, f));
F := FunctionField(C);


b2 := a1^2 + 4 * a2;
b4 := 2*a4 + a1*a3;
b6 := a3^2 + 4*a6;
b8 := a1^2 * a6  + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2;

Ebc_disc := -b2^2 * b8 - 8*b4^3 -27*b6^2 + 9 * b2 * b4 * b6;

cusps := Support(Divisor(Ebc_disc));
printf "Computed total of %o cusp places over finite field with %o elements.\n", #cusps, p;


/**
 *  For a cusp place on X1(32) over F_p, returns the index i in [1..#qcusps_twogen] of
 *  the unique rational cusp place that "lies over it", in the sense that the divisor obtained
 *  by reducing the rational cusp place over F_p is >= cusp.
 */
function LabelCusp(cusp)
    labels := [];

    for i in [1..#qcusps_twogen] do
        g1 := qcusps_twogen[i][1];
        g2 := qcusps_twogen[i][2];



        if Valuation(g1, cusp) gt 0 and Valuation(g2, cusp) gt 0 then
            Append(~labels, i);
        end if;
    end for;

    print labels;
    return labels[1];

    error "Unable to label cusp. This shouldn't mathematically be possible, check for code errors.";
end function;

// Compute divisors obtained by reducing rational cusp places modulo p
NULL_DIVISOR := cusps[1] - cusps[1];
qcusps_reduced := [NULL_DIVISOR : x in [1..qcusp_num]];

for cusp in cusps do
    label := LabelCusp(cusp);
    cusp_deg := Degree(cusp);
    qcusps_reduced[label] := qcusps_reduced[label] + cusp;
end for;
printf "Computed divisors over F%o obtained by reducing rational cusp places.\n", p;


cg, mapCgToDiv, mapDivToCg := ClassGroup(C);

// Find the subgroup of all divisor classes that can be obtained by reducing a cusp-supported rational divisor
qcuspsCgImage := mapDivToCg(qcusps_reduced);
qcuspsCgSubgrp, mapQCuspSubgrpToCg := sub<cg | qcuspsCgImage>;


places := [Places(C, i) : i in [1..9]];
places_count := [#places[i] : i in [1..9]];
assert places_count eq [12, 4, 8, 9, 16, 140, 312, 901, 1976];

function LoadCover(filename, places)
    lines := Split(Read(filename));
    T := [];

    for line in lines do
        if #line eq 0 then
            continue;
        end if;

        D := places[1][1] - places[1][1];

        // Example:
        // P_1_2,P_1_2,P_3_7
        // becomes ["P", "1", "2", "P", "1", "2", "P", "3", "7"].
        parts := Split(line, ",_");
        assert #parts mod 3 eq 0;

        for i in [1..#parts by 3] do
            assert parts[i] eq "P";

            deg := StringToInteger(parts[i + 1]);
            idx := StringToInteger(parts[i + 2]);

            D +:= places[deg][idx];
        end for;

        Append(~T, D);
    end for;

    return T;
end function;

divisors_to_check := LoadCover("X1_32/deg9_divisors_cover.txt", places);

printf "Loaded %o covering divisors.\n", #divisors_to_check;
assert #divisors_to_check eq 22606;
assert &and[Degree(D) le 18 : D in divisors_to_check];


has_func_of_deg_at_most_7 := false;
has_func_of_deg9 := false;
deg8_func_pole_divisors := {};

num_tasks := #divisors_to_check;
start_time := Realtime();
report_interval := Maximum(1, Floor(num_tasks / 200)); // roughly every 0.5%

printf "Searching RR spaces for %o divisors...\n", num_tasks;

SetColumns(200);

for idx -> D in divisors_to_check do
    R, m := RiemannRochSpace(D);

    for f in R do
        g := m(f);
        deg_f := Degree(g);

        if deg_f gt 0 and deg_f le 7 then
            has_func_of_deg_at_most_7 := true;

        elif deg_f eq 9 then
            has_func_of_deg9 := true;

        elif deg_f eq 8 then
            pole_divisor := Denominator(Divisor(g));
            Include(~deg8_func_pole_divisors, pole_divisor);
        end if;
    end for;

    if idx mod report_interval eq 0 or idx eq num_tasks then
        elapsed := Realtime() - start_time;
        eta := (num_tasks - idx) * (elapsed / idx);
        pct := (idx * 100.0) / num_tasks;

        printf "%o/%o (%o%%) | elapsed=%os | ETA=%os | deg8 pole divisors=%o | <=7=%o | deg9=%o\n",
            idx, num_tasks, RealField(4)!pct,
            Floor(elapsed), Floor(eta),
            #deg8_func_pole_divisors,
            has_func_of_deg_at_most_7,
            has_func_of_deg9;
    end if;
end for;

printf "Finished RR search.\n";
printf "has_func_of_deg_at_most_7 = %o\n", has_func_of_deg_at_most_7;
printf "has_func_of_deg9 = %o\n", has_func_of_deg9;
printf "#deg8_func_pole_divisors = %o\n", #deg8_func_pole_divisors;

save "x1_32.ws";