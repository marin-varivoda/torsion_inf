AttachSpec("../third_party/mdmagma/v2/mdmagma.spec");
SetSeed(1337);
SetColumns(200); // for prettier output without line wrapping

N := 32;
p := 3;
C := MDX1(N, GF(p));
FC := FunctionField(Curve(C));
g := 17; // genus of X1(32)

cusps := Cusps(C);
assert #cusps eq 18;

// These are the divisors obtained by reducing cusp places of X1(32) over Q modulo p
qcusps_reduced := CuspOrbitsQ(C);

// there are 17 rational cusps: 12 of deg 1, 2 of deg 2, 2 of degree 4 and 1 of degree 8
// we do some sanity checks here
assert #qcusps_reduced eq 17; 
assert #[D : D in qcusps_reduced | Degree(D) eq 1] eq 12;
assert #[D : D in qcusps_reduced | Degree(D) eq 2] eq 2;
assert #[D : D in qcusps_reduced | Degree(D) eq 4] eq 2;
assert #[D : D in qcusps_reduced | Degree(D) eq 8] eq 1;


cg, mapCgToDiv, mapDivToCg := ClassGroup(C);

// Find the subgroup of all divisor classes that can be obtained by reducing a cusp-supported rational divisor
qcuspsCgImage := mapDivToCg(qcusps_reduced);
qcuspsCgSubgrp, mapQCuspSubgrpToCg := sub<cg | qcuspsCgImage>;

// Calculate all places of degree up to 9 on C
places := [Places(Curve(C), i) : i in [1..9]];
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

        // Example: P_1_2,P_1_2,P_3_7 becomes ["P", "1", "2", "P", "1", "2", "P", "3", "7"].
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

// We load a set of divisors that covers all degree 9 effective divisors. In other words,
// if D >= 0 is a divisor of degree 9, then there exists E \in divisors_to_check such that D <= E
divisors_to_check := LoadCover("X1_32/deg9_divisors_cover.txt", places);

printf "Loaded %o covering divisors.\n", #divisors_to_check;
assert #divisors_to_check eq 22606;
// to make RR space search fast we want degree of each divisor in covering set to not be much greater than genus
assert &and[Degree(D) le 18 : D in divisors_to_check]; 



// Perform RR search
// - We will show that all non-constant functions f of degree <= 9, have degree 8
// - We find all possible pole divisors of degree 8 function
has_func_of_deg_at_most_7 := false;
has_func_of_deg9 := false;
deg8_func_pole_divisors := {};

printf "Searching RR spaces for %o divisors...\n", #divisors_to_check;

num_tasks := #divisors_to_check;
start_time := Realtime();
report_interval := Maximum(1, Floor(num_tasks / 200)); // roughly every 0.5%

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


// Verify that all degree 8 divisors can be obtained (up to equivalence) by reducing a cusp-supported rational divisor
non_obtainable_divisor_exists := false;
for D in deg8_func_pole_divisors do
    D_class := mapDivToCg(D);
    if not D_class in qcuspsCgSubgrp then
        non_obtainable_divisor_exists := true;
        error "Unexpectedly found degree 8 function pole divisor which can't be obtained by reducing a cusp-supported rational divisor.\n";
        break;
    end if;
end for;
if not non_obtainable_divisor_exists then
    printf "Successfully verified that all divisors in deg8_func_pole_divisors can be obtained by reducing a cusp-supported rational divisor.\n";
end if;

deg8_func_pole_divisors_list := Setseq(deg8_func_pole_divisors);
assert #deg8_func_pole_divisors eq 56;

min_lcm_degree := 10000;
for i in [1..#deg8_func_pole_divisors_list] do
    for j in [(i+1)..#deg8_func_pole_divisors_list] do
        A := deg8_func_pole_divisors_list[i];
        B := deg8_func_pole_divisors_list[j];

        D := LCM(A, B);  // smallest divisor with D >= A and D >= B
        assert D ge A;
        assert D ge B;

        min_lcm_degree := Min(min_lcm_degree, Degree(D));
    end for;
end for;

assert min_lcm_degree gt 9;
printf "Successfully verified that for all distinct A, B in deg8_func_pole_divisors, Degree(LCM(A,B)) > 9.\n";
printf "Minimum Degree(LCM(A,B)) value found was %o.\n", min_lcm_degree;

exit;