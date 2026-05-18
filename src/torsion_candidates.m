load "common.m";

/**
 *  Returns a finite list of pairs (M, N) that is a superset of \Phi^\infty(d).
 */
function TorsionInfCandidates(d)
    printf "Computing candidate torsion groups for extensions of degree %o over Q.\n", d;

    // List of all M such that EulerPhi(M) | d
    S := MWithEulerPhiDividing(d);

    /**
     * 1) Crude Abramovich Bound
     * Make a finite list of pairs (M, N) such that 325/2^16 * 6/pi^2 * M * N^2 <= 2 * d
     */
    step1_candidates := [];
    for M in S do
        N := M;
        abramovich_lb := Abramovich_lb(Gamma1_M_N_index_lb(M, N));

        while abramovich_lb le (2*d+0.5) do
            step1_candidates := Append(step1_candidates, [M, N]);

            N := N+M;
            abramovich_lb := Abramovich_lb(Gamma1_M_N_index_lb(M, N));
        end while;
    end for;
    printf "Pass 1 complete. Candidates remaining: %o\n", #step1_candidates;


    /**
     * 2) Keep (M, N) if exact Abramovich lower bound <= 2 * d / phi(M). 
     * If rk J(X_1(M,N)) is known to be 0, tighten criterion to <= d / phi(M).
     */
    step2_candidates := [];
    for MN_pair in step1_candidates do
        M := MN_pair[1];
        N := MN_pair[2];

        exact_bound := Abramovich_lb(Gamma1_M_N_index(M, N));
        phi_M := EulerPhi(M);
        
        // Default upper limit
        upper_limit := 2 * d / phi_M;
        
        // Tighten criterion if rank is known to be 0
        if IsX1MNRankKnownZero(M, N) then
            upper_limit := d / phi_M;
        end if;
        
        // This is an exact comparison with rational numbers so we don't need to
        // guard against floating-point fuzziness.
        if exact_bound le upper_limit then
            step2_candidates := Append(step2_candidates, MN_pair);
        end if;
    end for;
    printf "Pass 2 complete. Candidates remaining: %o\n", #step2_candidates;

    return step2_candidates;
end function;

/**
 * Prints a list of (M, N) candidates grouped nicely by M.
 */
procedure PrintCandidatesNicely(d, cands)
    printf "--- Candidates for degree %o (Total: %o) ---\n", d, #cands;
    
    // Extract unique M values
    M_set := { pair[1] : pair in cands };
    M_seq := SetToSequence(M_set);
    Sort(~M_seq);
    
    for M in M_seq do
        // Find all N associated with this M
        N_seq := [ pair[2] : pair in cands | pair[1] eq M ];
        Sort(~N_seq);
        
        // Format the output
        printf "  M = %2o : N in { ", M;
        for i in [1..#N_seq] do
            printf "%o", N_seq[i];
            if i lt #N_seq then
                printf ", ";
            end if;
        end for;
        printf " }\n";
    end for;
    printf "\n";
end procedure;

cand7 := TorsionInfCandidates(7);
cand7_final := [];
for pair in cand7 do
    M := pair[1];
    N := pair[2];
    
    if M eq 1 and not HasInfDegree7Points_X1N(N) then
        // printf "Removing [ %o, %o ] (Excluded by DvH14).\n", M, N;
    else
        Append(~cand7_final, pair);
    end if;
end for;
printf "Additional pass complete. Removed (1, N) pairs excluded by DvH14. Candidates remaining: %o\n", #cand7_final;
PrintCandidatesNicely(7, cand7_final);


cand8 := TorsionInfCandidates(8);
cand8_final := [];
for pair in cand8 do
    M := pair[1];
    N := pair[2];
    
    if M eq 1 and not HasInfDegree8Points_X1N(N) then
        // printf "Removing [ %o, %o ] (Excluded by DvH14).\n", M, N;
    else
        Append(~cand8_final, pair);
    end if;
end for;
printf "Additional pass complete. Removed (1, N) pairs excluded by DvH14. Candidates remaining: %o\n", #cand8_final;
PrintCandidatesNicely(8, cand8_final);

cand9_final := TorsionInfCandidates(9);
PrintCandidatesNicely(9, cand9_final);
