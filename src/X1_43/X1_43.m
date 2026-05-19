AttachSpec("../third_party/mdmagma/v2/mdmagma.spec");
SetSeed(1337);

N := 43;
p := 2;
C := MDX1(N, GF(p));
FC := FunctionField(Curve(C));
g := 57; // genus of X1(43)

cusps := Cusps(C);
cusps_count := #cusps;
cusps_degree := [Degree(cusp) : cusp in cusps];

plc1  := Places(Curve(C),  1);
plc7  := Places(Curve(C),  7);
plc8  := Places(Curve(C),  8);
plc9  := Places(Curve(C),  9);
plc10 := Places(Curve(C), 10);
plc11 := Places(Curve(C), 11);

plc1sum := &+ plc1;

// Note that there are 21 cusps of degree 1 and 3 cusps of degree 7
// [ 21, 0, 0, 0, 0, 0, 6, 42, 105, 84, 126 ]
print("The number of places of X1(43)(F2) of degrees 1..11 is:");
[#plc1] cat [#Places(Curve(C), i) : i in [2..6]] cat [#plc7, #plc8, #plc9, #plc10, #plc11];


// S := { D : D >=0, deg D = 18, and supp D has at least 7 different divisors of degree 1 }
// We build a list of divisors which dominates S up to diamond operator action.

// This means that for every s \in S, there should exist t \in divisors_to_check and a 
// diamond operator <a> such that <a> s <= t
divisors_to_check := [];


// 1) Cover divisors in S which have a place of degree 10 or 11 in support 

/**
 *  Takes a list of places which is closed under diamond operator action,
 *  partitions it into orbits and returns a list where there is exactly one representative
 *  for each diamond orbit.
 */
function getDiamondOrbitRepresentatives(place_list)
    diamond_orbits := [];
    rep := [];
    for plc in place_list do
        alredy_accounted_for := false;

        // check if plc already is in some orbit
        for diamond_orbit in diamond_orbits do
            if plc in diamond_orbit then
                alredy_accounted_for := true;
                break;
            end if;
        end for;

        if not alredy_accounted_for then
            Append(~diamond_orbits, DiamondOrbit(C, plc));
            Append(~rep, plc);
        end if;
    end for;

    return rep;
end function;

/* 1.1) Cover divisors in S with degree 11 place in support */
rep11 := getDiamondOrbitRepresentatives(plc11); // 6 elements
for P in rep11 do
    D := plc1sum + P;
    Append(~divisors_to_check, D);
end for;

/* 1.2) Cover divisors in S with degree 10 place in support */
rep10 := getDiamondOrbitRepresentatives(plc10); // 4 elements
for P in rep10 do
    D := 2 * plc1sum + P;
    Append(~divisors_to_check, D);
end for;


// 2) Cover divisors in S which have a place of degree 8 or 9 in support

// Observe that a divisor of this class has the form A + P where P is degree 8 or 9 place
// and A >= 0 is supported on places of degree 1.

// By using diamond operators which act transitively on places of degree 1, we can
// WLOG assume that plc1[1] has highest multiplicity in A, counting argument shows that this
// multiplicity is at most 4. 

// Also it is impossible that A has two different places with multiplicity >= 3 since that would 
// mean the degree of the entire A + P is at least >= (3+3+1+1+1+1+1) + (8) = 19 which is impossible.  

// Therefore A + P <= 2*plc1sum + 2*plc1[1] + P.

// Finally, we store all the divisors which can appear on the right side of above inequality:

for P in plc9 do
    D := 2*plc1sum + 2*plc1[1] + P;
    Append(~divisors_to_check, D);
end for;

for P in plc8 do
    D := 2*plc1sum + 2*plc1[1] + P;
    Append(~divisors_to_check, D);
end for;


// 3) Cover divisors in S which have a place of degree 7 in support

// The form is A + P again where P is a degree 7 place and A >= 0 is supported on places of degree 1.
// WLOG assume that plc1[1] has highest multiplicity in A, counting argument shows that this
// multiplicity is at most 5. 

// Moreover, at most two places in A can have multiplicity >= 3, and in that case both have 
// exactly the multiplicity 3. Note that (3+3+1+1+1+1+1) + (7) = 18

// Therefore A + P <= 2*plc1sum + 3*plc1[1] + plc1[i] + plc7[j], for some i,j 

for i in [2..#plc1] do
    for j in [1..#plc7] do
        D := 2*plc1sum + 3*plc1[1] + plc1[i] + plc7[j];
        Append(~divisors_to_check, D);
    end for;
end for;


// 4) Cover divisors in S which only have degree 1 places in support

// We use a separate C++ file to generate the cover in this case.
// The file X1_43/case4_cover.txt has one divisor per line:
//      a1 a2 ... a21
// representing a1*plc1[1] + ... + a21*plc1[21].

function LoadCase4Cover(filename, plc1)
    assert #plc1 eq 21;

    FP := Open(filename, "r");
    ret := [];
    line_no := 0;

    while true do
        line := Gets(FP);

        if IsEof(line) then
            break;
        end if;

        line_no +:= 1;

        // Skip empty lines, just in case.
        if #line eq 0 then
            continue;
        end if;

        coeffs := StringToIntegerSequence(line);

        assert #coeffs eq #plc1;
        assert &and[c ge 0 : c in coeffs];

        D := &+[ coeffs[i] * plc1[i] : i in [1..#plc1] ];

        assert Degree(D) le 57;

        Append(~ret, D);
    end while;

    return ret;
end function;

case4_cover := LoadCase4Cover("X1_43/case4_cover.txt", plc1);
assert #case4_cover eq 694; // make sure the parsed size matches C++ output

printf "Loaded case 4 cover of size: %o\n", #case4_cover; 

divisors_to_check := divisors_to_check cat case4_cover;

/**
 *  Returns true/ false depending on whether or not RR Space has a non-constant function
 *  of degree at most 18
 */
function RRSpaceHasFuncOfDegAtMost18(D)
    R, m := RiemannRochSpace(D);
    for f in R do
        g := m(f);
        deg_f := Degree(g);
        if deg_f gt 0 and deg_f le 18 then
            return true;
        end if;
    end for;
    return false;
end function;
