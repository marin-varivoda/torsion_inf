AttachSpec("../third_party/mdmagma/v2/mdmagma.spec");
SetSeed(1337);

N := 57;
p := 2;
C := MDX1(N, GF(p));
FC := FunctionField(Curve(C));
g := 85; // genus of X1(57)

// 18 cusps of degree 1, 9 cusps of degree 2, 2 cusps of degree 9, 1 cusp of degree 18
cusps := Cusps(C);
cusps_count := #cusps;
cusps_degree := [Degree(cusp) : cusp in cusps];

plc1  := Places(Curve(C),  1);
plc2  := Places(Curve(C),  2);
plc6  := Places(Curve(C),  6);
plc7  := Places(Curve(C),  7);
plc8  := Places(Curve(C),  8);
plc9  := Places(Curve(C),  9);
plc10 := Places(Curve(C), 10);
plc11 := Places(Curve(C), 11);
plc12 := Places(Curve(C), 12);

plc1sum := &+ plc1;
plc2sum := &+ plc2;
plc6sum := &+ plc6;

// [ 18, 9, 0, 0, 0, 2, 36, 18, 4, 216, 126, 381 ]
print("The number of places of X1(57)(F2) of degrees 1..12 is:");
[#plc1, #plc2] cat [#Places(Curve(C), i) : i in [3..5]] cat [#plc6, #plc7, #plc8, #plc9, #plc10, #plc11, #plc12];


/**
 *  Takes a list of places which is closed under diamond operator action,
 *  partitions it into orbits and returns a list where there is exactly one representative
 *  for each diamond orbit.
 */
function getDiamondOrbitRepresentatives(place_list)
    diamond_orbits := [];
    rep := [];
    for plc in place_list do
        already_accounted_for := false;

        // check if plc already is in some orbit
        for diamond_orbit in diamond_orbits do
            if plc in diamond_orbit then
                already_accounted_for := true;
                break;
            end if;
        end for;

        if not already_accounted_for then
            Append(~diamond_orbits, DiamondOrbit(C, plc));
            Append(~rep, plc);
        end if;
    end for;

    return rep;
end function;

// S := { D : D >=0, deg D = 18, and supp D has at least 6 different divisors of degree 1 }
// We build a list of divisors which dominates S up to diamond operator action.

// This means that for every s \in S, there should exist t \in divisors_to_check and a 
// diamond operator <a> such that <a> s <= t
divisors_to_check := [];


// 1) Cover divisors in S which have a place of degree 10, 11 or 12 in support

// Divisor of this class has the form D = A + B + P where P is the degree 10, 11 or 12 place, 
// A>=0 is supported on places of degree 1, and B>=0 is supported on places of degree 2.
// Note that A has to have at least 6 elements in support, thus the largest multiplicity appearing
// in A is at most 3. Indeed, (3+1+1+1+1+1) + (0) + (10) = 18

// therefore A <= 3*plc1sum and similarly B <= plc2sum, therefore D <= 3*plc1sum + plc2sum + P

// WLOG we apply diamond operators to D so that P is one of the 12+7+23=42 orbit representatives calculated below:

rep10 := getDiamondOrbitRepresentatives(plc10); // 12 elements
rep11 := getDiamondOrbitRepresentatives(plc11); //  7 elements
rep12 := getDiamondOrbitRepresentatives(plc12); // 23 elements

rep_10_11_12 := rep10 cat rep11 cat rep12;
for P in rep_10_11_12 do
    D := 3*plc1sum + plc2sum + P;
    Append(~divisors_to_check, D);
end for;


// 2) Cover divisors in S which have a place of degree 9 in support

// Divisor of this class has the form D = A + B + P where P is the degree 9 place, 
// A>=0 is supported on places of degree 1, and B>=0 is supported on places of degree 2.

// WLOG we apply diamond operators to D so that plc1 has the highest multiplicity in A

// Note that A has to have at least 6 elements in support, thus the largest multiplicity appearing
// in A is at most 4. Indeed, (4+1+1+1+1+1) + (0) + (9) = 18. Moreover, A has at most one place
// with multiplicity >= 3 (since otherwise we would have deg A >= 3+3+1+1+1+1=10 which is too much).

// Therefore A <= 2*plc1 + 2*plc1sum

// Moreover 18 = degA + deg B + deg P >= 6 + deg B + 9, which implies 3 >= deg B. Since deg B is even,
// this means 2 >= deg B. Therefore, B <= plc2sum

// Finally, we have: D = A + B + P <= 2*plc1[1] + 2*plc1sum + plc2sum + P

// By varying P, we get four different possible divisors on RHS but since degree isn't too big, we cover these
// four possibilities by grouping them into pairs and adding covering divisors with two different degree 9 places:

Append(~divisors_to_check, 2*plc1[1] + 2*plc1sum + plc2sum + plc9[1] + plc9[2]);
Append(~divisors_to_check, 2*plc1[1] + 2*plc1sum + plc2sum + plc9[3] + plc9[4]);


// 3) Cover divisors in S which have a place of degree 8 in support

// First, note that the 18 degree 8 places belong to the same diamond orbit!
rep8 := getDiamondOrbitRepresentatives(plc8); // 1 element
assert #rep8 eq 1;

// Again, we write D = A + B + P where P now has degree 8.
// We apply the diamond operator to D so that P = rep8[1]!

//  3.1) deg B = 4 
// We have D = A + B + P <= plc1sum + 2*plc2sum + rep8[1]
Append(~divisors_to_check, plc1sum + 2*plc2sum + rep8[1]);

//  3.2) deg B = 2
// We have deg A = 8, therefore A <= 3*plc1sum. Thus, D <= 3*plc1sum + plc2sum + rep8[1]
Append(~divisors_to_check, 3*plc1sum + plc2sum + rep8[1]);

//  3.3) deg B = 0
// We have deg A = 10. Let plc1[i] be the place with highest multiplicity in A. Then A <= 2*plc1[i] + 3*plc1sum.
// This is because 5 is the highest possible multiplicity in A, and at most one multiplicity can be >= 4 for the degree
// of A to not exceed 10.
// Therefore, D <= 2*plc1[i] + 3*plc1sum + rep8[1]

// We group the possible plc1[i] into two blocks of nine to keep divisors_to_check small.
Append(~divisors_to_check, 3*plc1sum + rep8[1] + 2*plc1[1] + 2*plc1[2] + 2*plc1[3] 
                                               + 2*plc1[4] + 2*plc1[5] + 2*plc1[6] 
                                               + 2*plc1[7] + 2*plc1[8] + 2*plc1[9] );

Append(~divisors_to_check, 3*plc1sum + rep8[1] + 2*plc1[10] + 2*plc1[11] + 2*plc1[12]
                                               + 2*plc1[13] + 2*plc1[14] + 2*plc1[15]
                                               + 2*plc1[16] + 2*plc1[17] + 2*plc1[18] );



// 4) Cover divisors in S which have a place of degree 7 in support

// First, note that the 36 degree 7 places partition into two diamond orbits.
rep7 := getDiamondOrbitRepresentatives(plc7); // 2 elements
assert #rep7 eq 2;
rep7sum := &+ rep7;

// Again, we write D = A + B + P where P is degree 7.
// We apply the diamond operator to D so that P \in {rep7[1], rep7[2]}!

//  4.1) deg B = 4 
// We have D = A + B + P <= 2*plc1sum + 2*plc2sum + rep7sum
Append(~divisors_to_check, 2*plc1sum + 2*plc2sum + rep7sum);

//  4.2) deg B = 2
// We have deg A = 9. Let plc1[i] be the place with highest multiplicity in A. Then A <= 2*plc1[i] + 2*plc1sum.
// This is because 4 is the highest possible multiplicity in A, and at most one multiplicity can be >= 3 for the degree
// of A to not exceed 9.
// We have D = A + B + P <= 2*plc1[i] + 2*plc1sum + plc2sum + rep7sum

Append(~divisors_to_check, 2*plc1sum + plc2sum + rep7sum + 2*plc1[1] + 2*plc1[2] + 2*plc1[3] 
                                                         + 2*plc1[4] + 2*plc1[5] + 2*plc1[6] 
                                                         + 2*plc1[7] + 2*plc1[8] + 2*plc1[9] );

Append(~divisors_to_check, 2*plc1sum + plc2sum + rep7sum + 2*plc1[10] + 2*plc1[11] + 2*plc1[12]
                                                         + 2*plc1[13] + 2*plc1[14] + 2*plc1[15]
                                                         + 2*plc1[16] + 2*plc1[17] + 2*plc1[18] );

//  4.3) deg B = 0
// We have deg A = 11. Let plc1[i] be the place with highest multiplicity in A. Then A <= 3*plc1[i] + 3*plc1sum.
// This is because 6 is the highest possible multiplicity in A, and at most one multiplicity can be >= 4 for the degree
// of A to not exceed 11 (note that 4+4+1+1+1+1=12).
// Therefore, D <= 3*plc1[i] + 3*plc1sum + rep7sum

Append(~divisors_to_check, 3*plc1sum + rep7sum + 3*plc1[1] + 3*plc1[2] + 3*plc1[3] 
                                               + 3*plc1[4] + 3*plc1[5] + 3*plc1[6] );

Append(~divisors_to_check, 3*plc1sum + rep7sum + 3*plc1[7]  + 3*plc1[8]  + 3*plc1[9] 
                                               + 3*plc1[10] + 3*plc1[11] + 3*plc1[12] );

Append(~divisors_to_check, 3*plc1sum + rep7sum + 3*plc1[13] + 3*plc1[14] + 3*plc1[15] 
                                               + 3*plc1[16] + 3*plc1[17] + 3*plc1[18] );



// 5) Cover divisors in S which have a place of degree 6 in support

// Divisor of this class has the form D = A + B + C where A>=0, B>=0, C>=0 are effective divisors
// supported on places of degree 1, degree 2 and degree 6 respectively. We consider all the possible
// 3-tuples (deg A, deg B, deg C) and do the case split accordingly:

//  5.1)    (6, 0, 12)
//  In this case we have D <= plc1sum + 2*plc6sum
Append(~divisors_to_check, plc1sum + 2*plc6sum);


//  5.2)    (6, 6, 6)
//  In this case we have D <= plc1sum + 3*plc2sum + plc6sum
Append(~divisors_to_check, plc1sum + 3*plc2sum + plc6sum);


//  5.3)    (8, 4, 6)
//  Apply the diamond operator to D so that plc1[1] has highest multiplicity in A. We then have
//  A <= plc1[1] + 2*plc1sum (highest multiplicity at most 3 and only one place with multiplicity 3 possible)
//  Therefore, we have D <= plc1[1] + 2*plc1sum + 2*plc2sum + plc6sum
Append(~divisors_to_check, plc1[1] + 2*plc1sum + 2*plc2sum + plc6sum);


//  5.4)    (10, 2, 6)
//  Apply the diamond operator to D so that plc1[1] has highest multiplicity in A. We then have
//  A <= 2*plc1[1] + 3*plc1sum (see case 3.3 for reasoning)
//  Therefore, we have D <= 2*plc1[1] + 3*plc1sum + plc2sum + plc6sum
Append(~divisors_to_check, 2*plc1[1] + 3*plc1sum + plc2sum + plc6sum);


//  5.5)    (12, 0, 6)
//  Apply the diamond operator to D so that plc1[1] has highest multiplicity in A. We then have
//  A <= 3*plc1[1] + 4*plc1sum (highest multiplicity at most 7 and only one place with multiplicity >=5 possible)
Append(~divisors_to_check, 3*plc1[1] + 4*plc1sum + plc6sum);
