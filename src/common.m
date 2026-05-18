/**
 * Constructs the congruence subgroup Gamma1(M, N), defined as the 
 * intersection of Gamma1(N) and the principal congruence subgroup Gamma(M).
 */
function Gamma1_M_N(M, N)
    assert IsDivisibleBy(N, M);

    grp1 := Gamma1(N);
    grp2 := CongruenceSubgroup(M);

    return Intersection(grp1, grp2);
end function;

/**
 * Finds the index [PSL(2,Z) : Gamma1(N)].
 * Returns the same value but is faster alternative to Index(Gamma1(N)) Magma call .
 */
function Gamma1_index(N)
    assert N ge 1;
    
    // Calculate the index in SL(2, Z)
    sl2_index := N^2;
    for p in PrimeDivisors(N) do
        sl2_index := (sl2_index div p^2) * (p^2 - 1);
    end for;
    
    // Project to PSL(2, Z)
    if N le 2 then
        return sl2_index;
    else
        return sl2_index div 2;
    end if;
end function;

/**
 * Finds the index [PSL(2,Z) : Gamma1(M, N)]. 
 * Returns the same value but is faster alternative to Index(Gamma1_M_N(M, N)) Magma call.
 */
function Gamma1_M_N_index(M, N)
    assert IsDivisibleBy(N, M);

    return M * Gamma1_index(N);
end function;

/**
 * Returns the Abramovich lower bound for the C-gonality of the modular curve X_G,
 * given the index of the modular subgroup G.
 */
function Abramovich_lb(congr_subgrp_index)
    return (325/2^15) * congr_subgrp_index;
end function;


/**
 *  Use the fact that infinite product over primes of (p^2-1) / p^2 is exactly 6 / pi^2
 *  to establish a lower bound. Crucially Gamma1_index_lb(N) is strictly increasing in N.  
 */
function Gamma1_index_lb(N)
    assert N ge 1;
    
    sl2_index_lb := N^2 * 6 / Pi(RealField())^2;
    
    return sl2_index_lb / 2;
end function;


/**
 * For fixed M this is a strictly increasing function in N that serves as a strict lower bound
 * for Gamma1_M_N.
 */
function Gamma1_M_N_index_lb(M, N)
    assert IsDivisibleBy(N, M);

    return M * Gamma1_index_lb(N);
end function;


/**
 * Checks if the rank of J(X_1(M, N)) over Q(zeta_M) is known to be zero.
 * * Sources:
 * For M = 1 and M = 2: 
 * Derickx, Etropolski, van Hoeij, Morrow, and Zureick-Brown, 
 * "Sporadic cubic torsion", Algebra & Number Theory 15 (2021), 
 * no. 7, Theorem 3.1.
 * * For M >= 3: 
 * M. Derickx and A. V. Sutherland, "Torsion subgroups of elliptic 
 * curves over quintic and sextic number fields", Proc. Amer. Math. Soc. 
 * 145 (2017), no. 10, Theorem 4.1.
 * * Note: A return value of false does not mean the rank is strictly positive, 
 * only that it is not covered by these explicit lists.
 */
function IsX1MNRankKnownZero(M, N)
    assert IsDivisibleBy(N, M);
    n := N div M;
    
    if M eq 1 then
        // Derived from S0 minus {63, 80, 95, 104, 105, 126, 144}
        rank_zero_N := {1..36} join {38..42} join {44..52} join 
            {54, 55, 56, 59, 60, 62, 64, 66, 68, 69, 70, 71, 72, 75, 76, 78, 
             81, 84, 87, 90, 94, 96, 98, 100, 108, 110, 119, 120, 132, 140, 
             150, 168, 180};
        return n in rank_zero_N;
    elif M eq 2 then
        // Derived from S1
        rank_zero_M2_N := {1..21} join {24, 25, 26, 27, 30, 33, 35, 42, 45};
        return n in rank_zero_M2_N;
    elif M eq 3 and n le 10 then
        return true;
    elif M eq 4 and n le 6 then
        return true;
    elif M eq 5 and n le 4 then
        return true;
    elif M eq 6 and n le 5 then
        return true;
    else
        return false;
    end if;
end function;

/**
 * Checks if X_1(N) has infinitely many points of degree = 7 over Q,
 * which is equivalent to (1, N) \in \Phi^\infty(7).
 *  
 * Based on M. Derickx and M. van Hoeij, "Gonality of the modular curve X1(N)",
 * J. Algebra 417 (2014), 52-71, Theorem 3.
 */
function HasInfDegree7Points_X1N(N)
    assert N ge 1;
    
    // 1..24, 26, 27, 28, 30
    if N le 24 then
        return true;
    elif N eq 26 or N eq 27 or N eq 28 or N eq 30 then
        return true;
    else 
        return false;
    end if;
end function;


/**
 * Checks if X_1(N) has infinitely many points of degree = 8 over Q,
 * which is equivalent to (1, N) \in \Phi^\infty(8).
 *  
 * Based on M. Derickx and M. van Hoeij, "Gonality of the modular curve X1(N)",
 * J. Algebra 417 (2014), 52-71, Theorem 3.
 */
function HasInfDegree8Points_X1N(N)
    assert N ge 1;
    
    if N le 28 then
        return true;
    elif N eq 30 or N eq 32 or N eq 36 then
        return true;
    else
        return false;
    end if;
end function;


/**
 * Returns a sorted list of all integers M such that EulerPhi(M) divides k.
 * Uses the lower bound EulerPhi(M) >= sqrt(M/2) to limit the search 
 * space to M <= 2 * k^2.
 * 
 * See for lower bound on EulerPhi(M): 
 * https://math.stackexchange.com/questions/527946/prove-that-phin-geq-sqrtn-2
 */
function MWithEulerPhiDividing(k)
    // The theoretical maximum M we need to check
    max_M := 2 * k^2;
    
    return [M : M in [1..max_M] | IsDivisibleBy(k, EulerPhi(M))];
end function;
