AttachSpec("/home/mvarivoda/mdmagma/v2/mdmagma.spec");


N := 32;

C := MDX1(N, Rationals());
C_proj := Curve(C);

// Reduce the curve over finite field
p := 3;
C_fin := Curve(Reduction(C_proj, p));

printf "Initialized the model for X_1(%o) and reduced it over a finite field with %o elements.\n", N, p;
printf "Checking if reduced X_1(%o) has a function of degree at most 9... (this might take a moment)\n", N;

plcs := [];
for i in [1..9] do
    printf "Starting Places(Curve(C), %o) [%o/12]...\n", i, i;
    t0 := Cputime();

    Append(~plcs, Places(Curve(C_fin), i));

    printf "Finished Places(Curve(C), %o) [%o/12] in %o seconds. Found %o places.\n",
        i, i, Cputime(t0), #plcs[#plcs];
end for;


Div := DivisorGroup(Curve(C_fin));
T_divisors := [];

filename := "covering_set_T_deg8.txt";
printf "\nOpening %o to load the dominating set T...\n", filename;
t0 := Cputime();

// Open file for reading
F := Open(filename, "r");

count := 0;
line := Gets(F);

while not IsEof(line) do
    if #line gt 0 then
        tokens := Split(line, ",");
        
        // Build the divisor by summing the places
        d := Div ! 0;
        for token in tokens do
            // Parse strings like "P_9_1976" -> parts = ["P", "9", "1976"]
            parts := Split(token, "_");
            deg := StringToInteger(parts[2]);
            idx := StringToInteger(parts[3]);
            
            d +:= Div ! plcs[deg][idx];
        end for;
        
        Append(~T_divisors, d);
    end if;
    
    count +:= 1;
    if count mod 5000 eq 0 then
        printf "  Parsed %o divisors...\n", count;
    end if;
    
    line := Gets(F);
end while;

delete F; // Close the file

printf "Successfully loaded %o divisors into Magma in %o seconds.\n", #T_divisors, Cputime(t0);


function RRSpaceHasFuncOfDegExactly9(D)
    R, m := RiemannRochSpace(D);
    for f in R do
        g := m(f);
        deg_f := Degree(g);
        if deg_f eq 9 then
            return true;
        end if;
    end for;
    return false;
end function;

picGrp, mapPicToDiv, mapDivToPic := ClassGroup(C_fin);
print "Successfully computed class group for curve over finite field";

function RRSpaceDeg8Numerators(D)
    numerators := [];
    R, m := RiemannRochSpace(D);
    for f in R do
        g := m(f);
        deg_f := Degree(g);
        if deg_f eq 8 then
            num_divisor := Numerator(Divisor(g));
            pic_grp_el := mapDivToPic(num_divisor);
            Append(~numerators, pic_grp_el);
        end if;
    end for;

    return Set(numerators);
end function;


results := [];
total := #T_divisors;

// Calculate how many iterations make up 0.5%
step_size := Max(1, Floor(total * 0.005));

// Start the wall-clock timer
t_start := Realtime();

for i := 1 to total do
    D := T_divisors[i];
    res_set := RRSpaceDeg8Numerators(D);
    for el in res_set do
        Append(~results, el);
    end for;
    
    // Check if we hit the 0.5% threshold
    if i mod step_size eq 0 then
        elapsed := Realtime(t_start);
        
        // Calculate the rate and ETA
        time_per_iter := elapsed / i;
        rem_iters := total - i;
        eta := time_per_iter * rem_iters;
        
        // Calculate the number of unique elements found so far
        unique_count := #Set(results);
        
        // Print dynamically on the same line
        printf "Progress: %o / %o | Unique: %o | Elapsed: %o s | ETA: %o s    \n", 
            i, total, (i * 100.0) / total, unique_count, elapsed, eta;
    end if;
end for;

// Final print out to preserve the last line and show the final unique count
printf "\nDone! Total compute time: %.1f s\n", Realtime(t_start);
printf "Total unique Riemann-Roch spaces found: %o\n", #Set(results);

//results := @@PARALLEL_SOCKET_MAP[ T_divisors | RRSpaceDeg8Numerators ]@@;

print results;

save "x1_32_nodeg9.ws";
