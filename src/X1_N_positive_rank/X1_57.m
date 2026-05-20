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


// 6) Cover divisors in S which are only supported on places of degree 1 and degree 2

// Divisor of this class has the form D = A + B where A>=0 and B>=0 are effective divisors
// supported on places of degree 1 and degree 2 respectively. We consider all the possible
// pairs (deg A, deg B) and do the case split accordingly:

//  6.1)    (6, 12)
//  In this case we have A <= plc1sum. The divisor B has total multiplicity 6 on degree 2 places,
//  therefore B <= 3*plc2sum + 3*plc2[i], where plc2[i] is a place of highest multiplicity in B.
//  We have D <= plc1sum + 3*plc2sum + 3*plc2[i]
Append(~divisors_to_check, plc1sum + 3*plc2sum + 3*plc2[1] + 3*plc2[2] + 3*plc2[3]);
Append(~divisors_to_check, plc1sum + 3*plc2sum + 3*plc2[4] + 3*plc2[5] + 3*plc2[6]);
Append(~divisors_to_check, plc1sum + 3*plc2sum + 3*plc2[7] + 3*plc2[8] + 3*plc2[9]);


//  6.2)    (8, 10)
//  Apply the diamond operator to D so that plc1[1] has highest multiplicity in A. We then have
//  A <= plc1[1] + 2*plc1sum. The divisor B has total multiplicity 5 on degree 2 places,
//  therefore B <= 2*plc2sum + 3*plc2[i], where plc2[i] is a place of highest multiplicity in B.
//  We have D <= plc1[1] + 2*plc1sum + 2*plc2sum + 3*plc2[i]
Append(~divisors_to_check, plc1[1] + 2*plc1sum + 2*plc2sum + 3*plc2[1] + 3*plc2[2]);
Append(~divisors_to_check, plc1[1] + 2*plc1sum + 2*plc2sum + 3*plc2[3] + 3*plc2[4]);
Append(~divisors_to_check, plc1[1] + 2*plc1sum + 2*plc2sum + 3*plc2[5] + 3*plc2[6]);
Append(~divisors_to_check, plc1[1] + 2*plc1sum + 2*plc2sum + 3*plc2[7] + 3*plc2[8]);
Append(~divisors_to_check, plc1[1] + 2*plc1sum + 2*plc2sum + 3*plc2[9]);


//  6.3)    (10, 8)
//  Apply the diamond operator to D so that plc1[1] has highest multiplicity in A. We then have
//  A <= 3*plc1[1] + 2*plc1sum + plc1[i] for some i. This is because 5 is the highest possible multiplicity in A,
//  and apart from plc1[1] at most one place can have multiplicity >=3.
//  Similarly, B <= 2*plc2sum + 2*plc2[j] for some j.
//  We have D <= 3*plc1[1] + 2*plc1sum + 2*plc2sum +  plc1[i] + 2*plc2[j] <= 3*plc1[1] + 2*plc1sum + 2*plc2sum + plc1block_xxx + 2*plc2block_yyy

//  We group the possible plc1[i] into blocks of size six, and the possible plc2[j] into blocks of size two.

plc1block_1_6   := plc1[1]  + plc1[2]  + plc1[3]  + plc1[4]  + plc1[5]  + plc1[6];
plc1block_7_12  := plc1[7]  + plc1[8]  + plc1[9]  + plc1[10] + plc1[11] + plc1[12];
plc1block_13_18 := plc1[13] + plc1[14] + plc1[15] + plc1[16] + plc1[17] + plc1[18];

plc2block_1_2 := plc2[1] + plc2[2];
plc2block_3_4 := plc2[3] + plc2[4];
plc2block_5_6 := plc2[5] + plc2[6];
plc2block_7_8 := plc2[7] + plc2[8];
plc2block_9   := plc2[9];

for R in [plc1block_1_6, plc1block_7_12, plc1block_13_18] do
    for Q in [plc2block_1_2, plc2block_3_4, plc2block_5_6, plc2block_7_8, plc2block_9] do
        Append(~divisors_to_check, 3*plc1[1] + 2*plc1sum + R + 2*plc2sum + 2*Q);
    end for;
end for;

// Note that if we had iterated over (i,j) pairs, we would need to add 17*9=153 more divisors. With this pair of blocks approach
// we get away with only adding 3*5=15 new divisors, while keeping the degree < 90.


//  6.4)    (12, 6)
//  Apply the diamond operator to D so that plc1[1] has highest multiplicity in A. We then have
//  A <= 4*plc1[1] + 3*plc1sum + plc1[i] for some i. This is because 7 is the highest possible multiplicity in A,
//  and apart from plc1[1] at most one place can have multiplicity >= 4.
//  Similarly, B <= plc2sum + 2*plc2[j] for some j.
//  We have D <= 4*plc1[1] + 3*plc1sum + plc1[i] + plc2sum + 2*plc2[j] <= 4*plc1[1] + 3*plc1sum + plc1block_xxx + plc2sum + 2*plc2[j]

//  We group the possible plc1[i] into two blocks of size nine, and we iterate over all plc2[j].
plc1block_1_9   := plc1[1]  + plc1[2]  + plc1[3]  + plc1[4]  + plc1[5]  + plc1[6]  + plc1[7]  + plc1[8]  + plc1[9];
plc1block_10_18 := plc1[10] + plc1[11] + plc1[12] + plc1[13] + plc1[14] + plc1[15] + plc1[16] + plc1[17] + plc1[18];

for R in [plc1block_1_9, plc1block_10_18] do
    for Q in plc2 do
        Append(~divisors_to_check, 4*plc1[1] + 3*plc1sum + R + plc2sum + 2*Q);
    end for;
end for;


//  6.5)    (14, 4)
//  Apply the diamond operator to D so that plc1[1] has highest multiplicity in A. We then have
//  A <= 6*plc1[1] + 3*plc1sum + 2*plc1[i] for some i. This is because 9 is the highest possible multiplicity in A, 
//  and apart from plc1[1] at most one place can have multiplicity >= 4.
//  Similarly, B <= plc2sum + plc2[j] for some j.
//  We have D <= 6*plc1[1] + 3*plc1sum + 2*plc1[i] + plc2sum + plc2[j] <= 6*plc1[1] + 3*plc1sum + 2*plc1block_xxx + plc2sum + plc2block_yyy

//  We group the possible plc1[i] into blocks of size two, and the possible plc2[j] into blocks of size three.

plc1block_1_2   := plc1[1]  + plc1[2];
plc1block_3_4   := plc1[3]  + plc1[4];
plc1block_5_6   := plc1[5]  + plc1[6];
plc1block_7_8   := plc1[7]  + plc1[8];
plc1block_9_10  := plc1[9]  + plc1[10];
plc1block_11_12 := plc1[11] + plc1[12];
plc1block_13_14 := plc1[13] + plc1[14];
plc1block_15_16 := plc1[15] + plc1[16];
plc1block_17_18 := plc1[17] + plc1[18];

plc2block_1_3 := plc2[1] + plc2[2] + plc2[3];
plc2block_4_6 := plc2[4] + plc2[5] + plc2[6];
plc2block_7_9 := plc2[7] + plc2[8] + plc2[9];

for R in [plc1block_1_2, plc1block_3_4, plc1block_5_6, plc1block_7_8, plc1block_9_10, plc1block_11_12, plc1block_13_14, plc1block_15_16, plc1block_17_18] do
    for Q in [plc2block_1_3, plc2block_4_6, plc2block_7_9] do
        Append(~divisors_to_check, 6*plc1[1] + 3*plc1sum + 2*R + plc2sum + Q);
    end for;
end for;


//  6.6)    (16, 2)
//  Apply the diamond operator to D so that plc1[1] has highest multiplicity in A. We then have
//  A <= 7*plc1[1] + 4*plc1sum + 2*plc1[i] for some i. This is because 11 is the highest possible multiplicity in A,
//  and apart from plc1[1] at most one place can have multiplicity greater than 4.
//  Then, B = plc2[j] for some j.
//  We have D <= 7*plc1[1] + 4*plc1sum + 2*plc1[i] + plc2[j] <= 7*plc1[1] + 4*plc1sum + 2*plc1block_xxx + plc2[j]

//  We group the possible plc1[i] into blocks of size four, and we iterate over all plc2[j].

plc1block_1_4   := plc1[1]  + plc1[2]  + plc1[3]  + plc1[4];
plc1block_5_8   := plc1[5]  + plc1[6]  + plc1[7]  + plc1[8];
plc1block_9_12  := plc1[9]  + plc1[10] + plc1[11] + plc1[12];
plc1block_13_16 := plc1[13] + plc1[14] + plc1[15] + plc1[16];

for R in [plc1block_1_4, plc1block_5_8, plc1block_9_12, plc1block_13_16, plc1block_17_18] do
    for Q in plc2 do
        Append(~divisors_to_check, 7*plc1[1] + 4*plc1sum + 2*R + Q);
    end for;
end for;


//  6.7)    (18, 0)
//  Apply the diamond operator to D so that plc1[1] has highest multiplicity in A. We then have
//  A <= 9*plc1[1] + 4*plc1sum + 3*plc1[i] + 3*plc1[j] for some i,j. This is because 13 is the highest
//  possible multiplicity in A, and apart from plc1[1] at most two places can have multiplicity >= 5. Also
//  the highest possible multiplicity for a place which is not plc1[1] is 7 (in case of 7+7+1+1+1+1=18).

//  We iterate over all possible pairs of plc1[i], plc1[j].
for i in [2..17] do
    for j in [i+1..18] do
        Append(~divisors_to_check, 9*plc1[1] + 4*plc1sum + 3*plc1[i] + 3*plc1[j]);
    end for;
end for;





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

/**
 *  Parallel implementation for RRSpaceHasFuncOfDegAtMost18. The idea is that 
 *  following two calls produce the same result:
 *  
 *  results := ParallelMapRRSpaceHasFuncOfDegAtMost18(divisors_to_check);
 *  results := [RRSpaceHasFuncOfDegAtMost18(D) : D in divisors_to_check];
 *  
 *  but ParallelMapRRSpaceHasFuncOfDegAtMost18 will do this work in parallel
 *  therefore being faster in practice on multicore systems.
 */
function ParallelMapRRSpaceHasFuncOfDegAtMost18(task_inputs)
    CORE_COUNT := 32;
    MEMORY_PER_WORKER := 16 * 10^9; // 0 for unlimited

    num_tasks := #task_inputs;
    if num_tasks eq 0 then
        return [];
    end if;

    worker_count := Minimum(CORE_COUNT, num_tasks);

    server_socket := Socket(: LocalHost := "localhost");
    host, port := Explode(SocketInformation(server_socket));

    // Start worker processes. Each worker receives a task index, computes the
    // corresponding result, and sends the pair "task_idx,result" back.
    for worker_idx in [1..worker_count] do
        child_pid := Fork();

        if child_pid eq 0 then
            SetMemoryLimit(MEMORY_PER_WORKER); 
            client_socket := Socket(host, port);
            
            // Ask for the first task. Index 0 means no completed task yet.
            Write(client_socket, "0,false");

            while true do
                _ := WaitForIO([client_socket]);
                msg := Read(client_socket);
                
                task_idx := -1;
                try
                    task_idx := eval msg;
                catch e
                    break;
                end try;

                // The parent sends -1 when there are no more tasks.
                if task_idx eq -1 then
                    break;
                end if;

                task_input := task_inputs[task_idx];
                task_result := RRSpaceHasFuncOfDegAtMost18(task_input);
                result_str := Sprintf("%o", task_result);

                try
                    Write(client_socket, Sprintf("%o,%o", task_idx, result_str));
                catch e
                    break;
                end try;
            end while;

            // Terminate the child process to prevent normal cleanup on exit. Because fork()
            // means processes share open file descriptors, letting the child quit normally
            // would affect the file descriptors of the parent, possibly interfering with the interpreter.
            System(Sprintf("kill -9 %o", Getpid()));
            quit;
 
        end if;
    end for;

    worker_sockets := [];
    for worker_idx in [1..worker_count] do
        Append(~worker_sockets, WaitForConnection(server_socket));
    end for;

    results_by_task_idx := AssociativeArray();
    next_task_idx := 1;
    tasks_completed := 0;

    start_time := Realtime();
    report_interval := Maximum(1, Floor(num_tasks / 200)); // roughly every 0.5%
    printf "Starting parallel computation with %o tasks across %o workers...\n", num_tasks, worker_count;

    while tasks_completed lt num_tasks do
        if #worker_sockets eq 0 then
            error "Parallel computation failed: no worker sockets remain.";
        end if;

        ready_sockets := WaitForIO(worker_sockets);

        for sock in ready_sockets do
            try
                msg := Read(sock);
            catch e
                error "Parallel computation failed while reading from a worker.";
            end try;
            
            task_idx := -1;
            task_result := false;

            try
                comma_idx := Index(msg, ",");
                if comma_idx le 1 or comma_idx ge #msg then
                    error "Invalid worker message.";
                end if;

                task_idx := eval msg[1..comma_idx-1];
                result_str := msg[comma_idx+1..#msg];

                if task_idx lt 0 or task_idx gt num_tasks then
                    error "Invalid task index.";
                end if;

                if result_str eq "true" then
                    task_result := true;
                elif result_str eq "false" then
                    task_result := false;
                else
                    error "Invalid task result.";
                end if;
            catch e
                error "Parallel computation failed: worker returned an invalid message.";
            end try;

            if task_idx gt 0 then
                results_by_task_idx[task_idx] := task_result;
                tasks_completed +:= 1;

                if tasks_completed mod report_interval eq 0 or tasks_completed eq num_tasks then
                    elapsed := Realtime() - start_time;
                    core_time_per_task := (elapsed * worker_count) / tasks_completed;
                    eta := (num_tasks - tasks_completed) * (elapsed / tasks_completed);
                    pct := (tasks_completed * 100.0) / num_tasks;
                    
                    // We round the duration values to int so the output fits on the same line
                    printf "%o/%o (%o%%) | elapsed=%os | ETA=%os | core-time/task=%os\n",
                            tasks_completed, num_tasks, RealField(4)!pct,
                            Floor(elapsed), Floor(eta), Floor(core_time_per_task);
                end if;
            end if;

            if next_task_idx le num_tasks then
                try
                    Write(sock, Sprintf("%o", next_task_idx));
                catch e
                    error "Parallel computation failed while sending a task to a worker.";
                end try;
                next_task_idx +:= 1;
            else
                try
                    Write(sock, "-1");
                catch e
                    error "Parallel computation failed while stopping a worker.";
                end try;

                // This worker has no more tasks to receive. Remove its socket from the active set
                Exclude(~worker_sockets, sock);
            end if;
        end for;
    end while;

    WaitForAllChildren();
    
    return [ results_by_task_idx[i] : i in [1..num_tasks] ];
end function;

printf "Searching RR spaces for %o divisors... (this might take a moment)\n", #divisors_to_check;

T := Time();
results := ParallelMapRRSpaceHasFuncOfDegAtMost18(divisors_to_check);

num_true := #[res : res in results | res];
num_false := #results - num_true;

printf "Computation finished. %o divisors produced a non-constant function of degree <= 18 (true), and %o did not (false).\n", num_true, num_false;
printf "Total calculation time: %o seconds\n", Time(T);

exit;