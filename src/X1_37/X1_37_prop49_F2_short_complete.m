/*
    X1_37_prop49_F2_short_complete.m

    Compact first-run version of the F_2, degree-9 Proposition 4.9
    translate test for X_1(37).

    The main line is intentionally close to the X0(112)/X0(117) scripts:

        1. build X_1(37) over F_2;
        2. load the hard-coded quotient map X_1(37) -> X_0(37)^+;
        3. get a verified order-5 class D0 from two quotient fibers;
        4. test all effective degree-9 divisors D against D+iD0, i=1,...,4.

    The safeguards from the completed run are retained:

        - optional reconstruction check of S37,T37 from the universal curve;
        - exact Pic^0 order-5 check for the quotient generator;
        - duplicate filtering for candidate fibers;
        - characteristic-2 degree-18 branch doubling candidates;
        - lazy fast-linear W0 test with recursive residue-coordinate flattening.

    Normal use:

        ./run_first.sh 16 16

    One worker:

        magma -b NCHUNKS:=16 CHUNK_INDEX:=0 \
              X1_37_prop49_F2_short_complete.m
*/

///////////////////////////////////////////////////////////////////////////
// Options.
///////////////////////////////////////////////////////////////////////////

function BoolOptionValue(value, option_name)
    if Type(value) eq BoolElt then
        return value;
    end if;
    if Type(value) eq RngIntElt then
        if value in {0, 1} then
            return value eq 1;
        end if;
        error Sprintf("Option %o must be Boolean-like.", option_name);
    end if;
    if Type(value) eq MonStgElt then
        if value in {"true", "True", "TRUE", "1", "yes", "Yes", "YES"} then
            return true;
        elif value in {"false", "False", "FALSE", "0", "no", "No", "NO"} then
            return false;
        end if;
        error Sprintf("Option %o must be true/false/1/0/yes/no.", option_name);
    end if;
    error Sprintf("Option %o must be Boolean-like; got type %o.", option_name, Type(value));
end function;

function IntegerOptionValue(value, option_name)
    if Type(value) eq RngIntElt then
        return value;
    end if;
    if Type(value) eq MonStgElt then
        try
            return StringToInteger(value);
        catch err
            error Sprintf("Option %o must be an integer; got %o.", option_name, value);
        end try;
    end if;
    error Sprintf("Option %o must be an integer; got type %o.", option_name, Type(value));
end function;

if not assigned RUN_TRANSLATE_TEST then
    RUN_TRANSLATE_TEST := true;
else
    RUN_TRANSLATE_TEST := BoolOptionValue(RUN_TRANSLATE_TEST, "RUN_TRANSLATE_TEST");
end if;

if not assigned VERIFY_MAP_BY_RECONSTRUCTION then
    VERIFY_MAP_BY_RECONSTRUCTION := false;
else
    VERIFY_MAP_BY_RECONSTRUCTION := BoolOptionValue(VERIFY_MAP_BY_RECONSTRUCTION, "VERIFY_MAP_BY_RECONSTRUCTION");
end if;

if not assigned VERIFY_DELTA_ORDER then
    VERIFY_DELTA_ORDER := true;
else
    VERIFY_DELTA_ORDER := BoolOptionValue(VERIFY_DELTA_ORDER, "VERIFY_DELTA_ORDER");
end if;

if not assigned VERIFY_MODEL_BASIC then
    VERIFY_MODEL_BASIC := true;
else
    VERIFY_MODEL_BASIC := BoolOptionValue(VERIFY_MODEL_BASIC, "VERIFY_MODEL_BASIC");
end if;

if not assigned QUIET_SETUP then
    QUIET_SETUP := false;
else
    QUIET_SETUP := BoolOptionValue(QUIET_SETUP, "QUIET_SETUP");
end if;

if not assigned F2_CACHED_SCAN_TARGETS_ONLY then
    F2_CACHED_SCAN_TARGETS_ONLY := true;
else
    F2_CACHED_SCAN_TARGETS_ONLY := BoolOptionValue(F2_CACHED_SCAN_TARGETS_ONLY, "F2_CACHED_SCAN_TARGETS_ONLY");
end if;

if not assigned VERIFY_FAST_LINEAR_W0_AGAINST_RR then
    VERIFY_FAST_LINEAR_W0_AGAINST_RR := false;
else
    VERIFY_FAST_LINEAR_W0_AGAINST_RR := BoolOptionValue(VERIFY_FAST_LINEAR_W0_AGAINST_RR, "VERIFY_FAST_LINEAR_W0_AGAINST_RR");
end if;

if not assigned FAST_LINEAR_W0_VERIFY_LIMIT then
    FAST_LINEAR_W0_VERIFY_LIMIT := 0;
else
    FAST_LINEAR_W0_VERIFY_LIMIT := IntegerOptionValue(FAST_LINEAR_W0_VERIFY_LIMIT, "FAST_LINEAR_W0_VERIFY_LIMIT");
end if;

if not assigned NCHUNKS then
    NCHUNKS := 1;
else
    NCHUNKS := IntegerOptionValue(NCHUNKS, "NCHUNKS");
end if;

if not assigned CHUNK_INDEX then
    CHUNK_INDEX := 0;
else
    CHUNK_INDEX := IntegerOptionValue(CHUNK_INDEX, "CHUNK_INDEX");
end if;

if not assigned VerboseEvery then
    VerboseEvery := 1000;
else
    VerboseEvery := IntegerOptionValue(VerboseEvery, "VerboseEvery");
end if;

if not assigned MaxSurvivors then
    MaxSurvivors := 1;
else
    MaxSurvivors := IntegerOptionValue(MaxSurvivors, "MaxSurvivors");
end if;

// Verified in the completed F_2 run.  Set both to 0 to search all pairs.
if not assigned CACHED_GENERATOR_I then
    CACHED_GENERATOR_I := 1;
else
    CACHED_GENERATOR_I := IntegerOptionValue(CACHED_GENERATOR_I, "CACHED_GENERATOR_I");
end if;

if not assigned CACHED_GENERATOR_J then
    CACHED_GENERATOR_J := 4;
else
    CACHED_GENERATOR_J := IntegerOptionValue(CACHED_GENERATOR_J, "CACHED_GENERATOR_J");
end if;

if NCHUNKS lt 1 then
    error "NCHUNKS must be positive.";
end if;
if CHUNK_INDEX lt 0 or CHUNK_INDEX ge NCHUNKS then
    error "CHUNK_INDEX must satisfy 0 <= CHUNK_INDEX < NCHUNKS.";
end if;

///////////////////////////////////////////////////////////////////////////
// Small utilities.
///////////////////////////////////////////////////////////////////////////

procedure CheckOrDie(condition, message)
    if not condition then
        error message;
    end if;
end procedure;

procedure PrintSection(title)
    print "";
    print "======================================================================";
    print title;
    print "======================================================================";
end procedure;

function PositivePartOfDivisor(D, C)
    Dpos := DivisorGroup(C)!0;
    for P0 in Support(D) do
        v := Valuation(D, P0);
        if v gt 0 then
            Dpos +:= v*P0;
        end if;
    end for;
    return Dpos;
end function;

function InW0FF(D)
    V, rrmap := RiemannRochSpace(D);
    return Dimension(V) gt 0;
end function;

function IsPrincipalDegreeZero(D)
    if Degree(D) ne 0 then
        return false;
    end if;
    return Dimension(RiemannRochSpace(D)) gt 0;
end function;

///////////////////////////////////////////////////////////////////////////
// X0-style effective-divisor enumeration, streamed by degree type.
///////////////////////////////////////////////////////////////////////////

function DegreeTypesOfDegreeFF(degree, C, place_cache)
    occurring_degrees := [i : i in [1..degree] | #place_cache[i] gt 0];
    old := [ <[<0,0>], degree> ];
    new := old;
    done := [];

    for d0 in Reverse(occurring_degrees) do
        for n in Reverse([1..Floor(degree/d0)]) do
            for degree_type in old do
                for i in [1..Floor(degree_type[2]/(n*d0))] do
                    dt := <degree_type[1] cat [<n,d0> : j in [1..i]],
                           degree_type[2] - n*d0*i>;
                    Append(~new, dt);
                end for;
            end for;
            done cat:= [dt[1][2..#dt[1]] : dt in new | dt[2] eq 0];
            new := [dt : dt in new | dt[2] gt 0];
            old := new;
        end for;
    end for;

    seen := AssociativeArray();
    keys := [];
    for dt in done do
        key := Sprint(Sort([Sprint(pair) : pair in dt]));
        if not IsDefined(seen, key) then
            seen[key] := dt;
            Append(~keys, key);
        end if;
    end for;
    return [seen[key] : key in Sort(keys)];
end function;

function DegreeTypeGroups(degree_type)
    groups := [];
    for t in degree_type do
        mult := t[1];
        degp := t[2];
        found := false;
        for i in [1..#groups] do
            if groups[i][1] eq mult and groups[i][2] eq degp then
                groups[i] := <mult, degp, groups[i][3] + 1>;
                found := true;
                break;
            end if;
        end for;
        if not found then
            Append(~groups, <mult, degp, 1>);
        end if;
    end for;
    return groups;
end function;

function CountDivisorsFromGroups(groups, place_cache)
    used_by_degree := AssociativeArray();
    ndivs := 1;

    for g in groups do
        degp := g[2];
        count := g[3];
        used := IsDefined(used_by_degree, degp) select used_by_degree[degp] else 0;
        available := #place_cache[degp] - used;
        if available lt count then
            return 0;
        end if;
        ndivs *:= Binomial(available, count);
        used_by_degree[degp] := used + count;
    end for;

    return ndivs;
end function;

///////////////////////////////////////////////////////////////////////////
// Fast W0 membership.
///////////////////////////////////////////////////////////////////////////

function ResidueCoordinatesOverBaseFF(c, Fbase)
    if Type(c) eq RngIntElt then
        return [Fbase!c];
    end if;
    if Parent(c) cmpeq Fbase then
        return [Fbase!c];
    end if;
    ok, c_base := IsCoercible(Fbase, c);
    if ok then
        return [c_base];
    end if;

    coords := [Fbase | ];
    for a in Eltseq(c) do
        coords cat:= ResidueCoordinatesOverBaseFF(a, Fbase);
    end for;
    return coords;
end function;

function LocalExpansionCoefficientsFF(h, P0, nterms)
    if nterms eq 0 then
        return [];
    end if;

    t := LocalUniformizer(P0);
    rem := h;
    coeffs := [];
    for i in [1..nterms] do
        c := Evaluate(rem, P0);
        if c cmpeq Infinity() then
            error "LocalExpansionCoefficientsFF: input is not regular at P0.";
        end if;
        Append(~coeffs, c);
        rem := (rem - Lift(c, P0))/t;
    end for;
    return coeffs;
end function;

function FastRowsForPlace(rr_functions, P0, e, rmax, Fbase)
    ncols := #rr_functions;
    if ncols eq 0 or rmax eq 0 then
        return [];
    end if;

    coeffs_by_basis := [];
    t := LocalUniformizer(P0);
    for f in rr_functions do
        Append(~coeffs_by_basis, LocalExpansionCoefficientsFF(f*t^e, P0, rmax));
    end for;

    rows := [];
    for n in [1..rmax] do
        for coord in [1..Degree(P0)] do
            row := [Fbase | ];
            for j in [1..ncols] do
                coords := ResidueCoordinatesOverBaseFF(coeffs_by_basis[j][n], Fbase);
                Append(~row, coord le #coords select coords[coord] else Fbase!0);
            end for;
            Append(~rows, row);
        end for;
    end for;

    return rows;
end function;

function DegreeTypeCoverDivisorFF(degree_type, C, place_cache)
    max_mult_by_degree := AssociativeArray();
    for t in degree_type do
        mult := t[1];
        degp := t[2];
        if (not IsDefined(max_mult_by_degree, degp)) or mult gt max_mult_by_degree[degp] then
            max_mult_by_degree[degp] := mult;
        end if;
    end for;

    H := DivisorGroup(C)!0;
    for degp in Keys(max_mult_by_degree) do
        for P0 in place_cache[degp] do
            H +:= max_mult_by_degree[degp]*P0;
        end for;
    end for;
    return H;
end function;

function BuildFastW0Data(C, H, delta)
    delta_plus := PositivePartOfDivisor(delta, C);
    delta_minus := PositivePartOfDivisor(-delta, C);
    G := H + delta_plus;
    Rmax := H + delta_minus;

    t_rr := Cputime();
    V, rrmap := RiemannRochSpace(G);
    rr_time := Cputime(t_rr);
    basis := Basis(V);
    rr_functions := [rrmap(b) : b in basis];

    places := [];
    rmaxs := [];
    eval_shifts := [];
    for P0 in Support(Rmax) do
        rmax := Valuation(Rmax, P0);
        if rmax gt 0 then
            Append(~places, P0);
            Append(~rmaxs, rmax);
            Append(~eval_shifts, Valuation(G, P0));
        end if;
    end for;

    return <#basis, places, rmaxs, rr_functions, eval_shifts,
            AssociativeArray(), Dimension(V), Degree(G), Degree(Rmax), rr_time>;
end function;

function FastW0RowCount(data)
    return &+[data[3][i]*Degree(data[2][i]) : i in [1..#data[2]]];
end function;

function FastW0DataWithCache(data, row_cache)
    return <data[1], data[2], data[3], data[4], data[5],
            row_cache, data[7], data[8], data[9], data[10]>;
end function;

function CachedRowsForFastW0Place(data, i, Fbase)
    row_cache := data[6];
    if IsDefined(row_cache, i) then
        return data, row_cache[i];
    end if;

    rows := FastRowsForPlace(data[4], data[2][i], data[5][i], data[3][i], Fbase);
    row_cache[i] := rows;
    return FastW0DataWithCache(data, row_cache), rows;
end function;

function FastInW0(D, data, Fbase)
    ncols := data[1];
    if ncols eq 0 then
        return false, data;
    end if;

    jobs := [];
    for i in [1..#data[2]] do
        P0 := data[2][i];
        keep_orders := data[3][i] - Valuation(D, P0);
        if keep_orders lt 0 then
            error "FastInW0: degree-type cover divisor does not contain D.";
        end if;
        keep_rows := keep_orders*Degree(P0);
        if keep_rows gt 0 then
            Append(~jobs, <keep_rows, i>);
        end if;
    end for;

    if #jobs eq 0 then
        return true, data;
    end if;

    rows := [];
    row_count := 0;
    for job in Reverse(Sort(jobs)) do
        keep_rows := job[1];
        data, place_rows := CachedRowsForFastW0Place(data, job[2], Fbase);
        rows cat:= place_rows[1..keep_rows];
        row_count +:= keep_rows;
        if row_count ge ncols then
            M := Matrix(Fbase, #rows, ncols, &cat rows);
            if Rank(M) eq ncols then
                return false, data;
            end if;
        end if;
    end for;

    M := Matrix(Fbase, #rows, ncols, &cat rows);
    return Rank(M) lt ncols, data;
end function;

procedure TestDivisor(~survivors, ~stop_search, ~checked, D,
        deltas, labels, ~fast_data, Fbase, t0)

    if stop_search then
        return;
    end if;

    checked +:= 1;
    survives := true;
    first_bad := 0;

    for j in [1..#deltas] do
        in_w0, data_j := FastInW0(D, fast_data[j], Fbase);
        fast_data[j] := data_j;

        if VERIFY_FAST_LINEAR_W0_AGAINST_RR and checked le FAST_LINEAR_W0_VERIFY_LIMIT then
            rr_in_w0 := InW0FF(D + deltas[j]);
            if rr_in_w0 ne in_w0 then
                error Sprintf("Fast W0 mismatch at divisor %o, delta %o: fast=%o, RR=%o",
                              checked, j, in_w0, rr_in_w0);
            end if;
        end if;

        if not in_w0 then
            survives := false;
            first_bad := j;
            break;
        end if;
    end for;

    if survives then
        Append(~survivors, D);
        printf "SURVIVOR %o in chunk %o after %o tested divisors: %o\n",
               #survivors, CHUNK_INDEX, checked, D;
        if MaxSurvivors gt 0 and #survivors ge MaxSurvivors then
            print "INCONCLUSIVE: stopped after reaching MaxSurvivors.";
            print "CHUNK_RESULT: INCONCLUSIVE";
            stop_search := true;
        end if;
    elif VerboseEvery gt 0 and (checked mod VerboseEvery eq 0) then
        printf "Chunk %o checked %o divisors; last killed by %o; survivors %o; elapsed %o\n",
               CHUNK_INDEX, checked, labels[first_bad], #survivors, Cputime(t0);
    end if;
end procedure;

procedure StreamGroupDivisors(~survivors, ~stop_search, ~checked,
        groups, group_index, chosen_in_group, next_place_index,
        currentD, used_places, place_cache, deltas, labels, ~fast_data,
        Fbase, t0)

    if stop_search then
        return;
    end if;

    if group_index gt #groups then
        TestDivisor(~survivors, ~stop_search, ~checked, currentD,
            deltas, labels, ~fast_data, Fbase, t0);
        return;
    end if;

    mult := groups[group_index][1];
    degp := groups[group_index][2];
    count := groups[group_index][3];

    if chosen_in_group eq count then
        StreamGroupDivisors(~survivors, ~stop_search, ~checked,
            groups, group_index + 1, 0, 1, currentD, used_places,
            place_cache, deltas, labels, ~fast_data, Fbase, t0);
        return;
    end if;

    places := place_cache[degp];
    remaining := count - chosen_in_group;
    max_place_index := #places - remaining + 1;
    if max_place_index lt next_place_index then
        return;
    end if;

    for ip in [next_place_index..max_place_index] do
        P0 := places[ip];
        if P0 in used_places then
            continue;
        end if;
        StreamGroupDivisors(~survivors, ~stop_search, ~checked,
            groups, group_index, chosen_in_group + 1, ip + 1,
            currentD + mult*P0, used_places cat [P0],
            place_cache, deltas, labels, ~fast_data, Fbase, t0);
        if stop_search then
            break;
        end if;
    end for;
end procedure;

function CheckNoTranslateFromDeltasFF(C, d, deltas, labels, quotient_map_degree)
    if #deltas eq 0 or #labels ne #deltas then
        print "CHUNK_RESULT: INCONCLUSIVE";
        return false, [];
    end if;

    if VERIFY_MODEL_BASIC then
        kC, kC_into_C := ExactConstantField(C);
        printf "Exact constant field size: %o\n", #kC;
        try
            printf "genus(C) = %o\n", Genus(C);
        catch err_genus
            print "WARNING: could not compute Genus(C); continuing.";
        end try;
    else
        print "Basic model checks skipped by VERIFY_MODEL_BASIC:=false.";
    end if;

    printf "quotient map degree = %o\n", quotient_map_degree;
    printf "Testing %o nonzero quotient element(s).\n", #deltas;
    printf "Type-chunking: CHUNK_INDEX=%o, NCHUNKS=%o\n", CHUNK_INDEX, NCHUNKS;

    t_places := Cputime();
    place_cache := AssociativeArray();
    for i in [1..d] do
        place_cache[i] := Places(C, i);
    end for;
    printf "Place counts by degree 1..%o: %o\n", d, [#place_cache[i] : i in [1..d]];

    degree_types := DegreeTypesOfDegreeFF(d, C, place_cache);
    assigned_types := [i : i in [1..#degree_types] | (i - 1) mod NCHUNKS eq CHUNK_INDEX];
    printf "Unique degree types: %o; assigned type indices: %o\n",
           #degree_types, assigned_types;
    printf "Place/type precomputation time: %o seconds.\n", Cputime(t_places);

    survivors := [];
    checked := 0;
    t0 := Cputime();
    Fbase := GF(2);

    for typeno in assigned_types do
        degree_type := degree_types[typeno];
        groups := DegreeTypeGroups(degree_type);
        ndivs := CountDivisorsFromGroups(groups, place_cache);
        printf "Degree type %o/%o: %o; groups %o; divisors %o; elapsed %o\n",
               typeno, #degree_types, degree_type, groups, ndivs, Cputime(t0);

        H := DegreeTypeCoverDivisorFF(degree_type, C, place_cache);
        fast_data := [* *];
        dims := [];
        rows := [];
        rr_times := [];
        for j in [1..#deltas] do
            data_j := BuildFastW0Data(C, H, deltas[j]);
            Append(~fast_data, data_j);
            Append(~dims, data_j[7]);
            Append(~rows, FastW0RowCount(data_j));
            Append(~rr_times, data_j[10]);
        end for;
        printf "Fast W0 setup: cover degree %o; dims %o; row maxima %o; RR times %o\n",
               Degree(H), dims, rows, rr_times;

        stop_search := false;
        StreamGroupDivisors(~survivors, ~stop_search, ~checked,
            groups, 1, 0, 1, DivisorGroup(C)!0, [], place_cache,
            deltas, labels, ~fast_data, Fbase, t0);

        if stop_search then
            return false, survivors;
        end if;
    end for;

    printf "Divisors tested in chunk %o: %o\n", CHUNK_INDEX, checked;
    printf "Surviving translates in this chunk: %o\n", #survivors;

    if #survivors eq 0 then
        print "SUCCESS: this chunk contains no translate survivor.";
        print "CHUNK_RESULT: SUCCESS";
        return true, survivors;
    end if;

    print "INCONCLUSIVE: at least one quotient coset survived in this chunk.";
    print "CHUNK_RESULT: INCONCLUSIVE";
    return false, survivors;
end function;

///////////////////////////////////////////////////////////////////////////
// 1. Build X_1(37) over F_2.
///////////////////////////////////////////////////////////////////////////

PrintSection("1. Build X_1(37) over F_2");

p := 2;
F := GF(p);
RF<x> := FunctionField(F);
RY<y> := PolynomialRing(RF);

f37 := (x + 1)^4*y^18
        + (x + 1)^3*(4*x^3 + 15*x^2 + 4*x - 8)*y^17
        + (x + 1)^3*(6*x^5 + 46*x^4 + 71*x^3 - 39*x^2 - 56*x + 28)*y^16
        + (x + 1)^3*(4*x^7 + 58*x^6 + 204*x^5 + 109*x^4 - 300*x^3 - 91*x^2 + 186*x - 56)*y^15
        + (x + 1)^2*(x^10 + 34*x^9 + 265*x^8 + 692*x^7 + 295*x^6 - 1055*x^5 - 745*x^4 + 570*x^3 + 135*x^2 - 205*x + 70)*y^14
        + (x + 1)*(7*x^12 + 130*x^11 + 755*x^10 + 1659*x^9 + 426*x^8 - 3357*x^7 - 3594*x^6 + 612*x^5 + 1760*x^4 + 274*x^3 + 32*x^2 + 6*x - 56)*y^13
        + (21*x^14 + 294*x^13 + 1420*x^12 + 2579*x^11 - 870*x^10 - 9749*x^9 - 11481*x^8 + 664*x^7 + 10483*x^6 + 6962*x^5 - 527*x^4 - 2730*x^3 - 733*x^2 + 322*x + 28)*y^12
        + (35*x^15 + 363*x^14 + 1210*x^13 + 414*x^12 - 6515*x^11 - 14529*x^10 - 6010*x^9 + 16913*x^8 + 23552*x^7 + 4413*x^6 - 13117*x^5 - 8675*x^4 + 2488*x^3 + 2270*x^2 - 494*x - 8)*y^11
        + (26*x^16 + 184*x^15 + 129*x^14 - 2320*x^13 - 8129*x^12 - 7355*x^11 + 12096*x^10 + 31592*x^9 + 16424*x^8 - 20987*x^7 - 29776*x^6 - 1252*x^5 + 14573*x^4 + 1991*x^3 - 3145*x^2 + 412*x + 1)*y^10
        - x*(x^17 + 5*x^16 + 52*x^15 + 501*x^14 + 2023*x^13 + 2681*x^12 - 4715*x^11 - 19581*x^10 - 18001*x^9 + 14660*x^8 + 40593*x^7 + 16842*x^6 - 25005*x^5 - 20440*x^4 + 8210*x^3 + 5975*x^2 - 2525*x + 210)*y^9
        - x*(x^18 + 10*x^17 + 52*x^16 + 172*x^15 + 162*x^14 - 1311*x^13 - 5622*x^12 - 7317*x^11 + 5435*x^10 + 25859*x^9 + 20839*x^8 - 16281*x^7 - 34475*x^6 - 2999*x^5 + 20366*x^4 + 2065*x^3 - 5725*x^2 + 1264*x - 66)*y^8
        + x*(3*x^17 + 32*x^16 + 180*x^15 + 595*x^14 + 789*x^13 - 1585*x^12 - 7639*x^11 - 9008*x^10 + 5183*x^9 + 22272*x^8 + 12561*x^7 - 16256*x^6 - 16253*x^5 + 7183*x^4 + 5741*x^3 - 3053*x^2 + 395*x - 12)*y^7
        - x*(3*x^16 + 39*x^15 + 243*x^14 + 807*x^13 + 1079*x^12 - 1278*x^11 - 6267*x^10 - 6155*x^9 + 4834*x^8 + 13100*x^7 + 1905*x^6 - 9941*x^5 - 1372*x^4 + 3726*x^3 - 978*x^2 + 72*x - 1)*y^6
        + x^2*(x^14 + 22*x^13 + 166*x^12 + 564*x^11 + 678*x^10 - 891*x^9 - 3469*x^8 - 2293*x^7 + 3656*x^6 + 4495*x^5 - 2053*x^4 - 2249*x^3 + 1275*x^2 - 180*x + 6)*y^5
        - 5*x^3*(x^11 + 12*x^10 + 42*x^9 + 36*x^8 - 111*x^7 - 270*x^6 - 37*x^5 + 331*x^4 + 80*x^3 - 182*x^2 + 48*x - 3)*y^4
        + x^4*(10*x^8 + 35*x^7 - 15*x^6 - 265*x^5 - 325*x^4 + 201*x^3 + 321*x^2 - 180*x + 20)*y^3
        + x^5*(26*x^4 + 77*x^3 + 30*x^2 - 72*x + 15)*y^2
        - x^6*(x^3 + 7*x^2 + 12*x - 6)*y
        + x^7;

K<yC> := FunctionField(f37);
xC := K!x;
C := K;

if VERIFY_MODEL_BASIC then
    try
        gK := Genus(K);
        printf "genus X_1(37) over F_2 = %o\n", gK;
        CheckOrDie(gK eq 40, "Unexpected genus for the X_1(37) model.");
    catch err_genus
        print "WARNING: Magma could not compute Genus(K). Continuing.";
    end try;
    Kconst, Kconst_in_K := ExactConstantField(K);
    printf "exact constant field size = %o\n", #Kconst;
    CheckOrDie(#Kconst eq #F, "The exact constant field is larger than F_2.");
else
    print "Basic genus and constant-field checks skipped.";
end if;

function ElementFromCoefficientList(coeffs)
    if #coeffs eq 0 then
        return K!0;
    end if;
    return &+[K!(coeffs[i+1])*yC^i : i in [0..#coeffs-1]];
end function;

///////////////////////////////////////////////////////////////////////////
// 2. Load and verify the quotient map data.
///////////////////////////////////////////////////////////////////////////

PrintSection("2. Load X_0(37)^+ map data");

load "X1_37_X0plus_map_data.m";
CheckOrDie(assigned S37_hard_coeffs and assigned T37_hard_coeffs,
           "X1_37_X0plus_map_data.m must define S37_hard_coeffs and T37_hard_coeffs.");

S37 := ElementFromCoefficientList(S37_hard_coeffs);
T37 := ElementFromCoefficientList(T37_hard_coeffs);
CheckOrDie(S37 ne K!0 and T37 ne K!0, "The loaded quotient map data is zero.");
printf "#S37 coefficients = %o; #T37 coefficients = %o\n", #S37_hard_coeffs, #T37_hard_coeffs;

if VERIFY_MAP_BY_RECONSTRUCTION then
    PrintSection("2a. Reconstruct S37,T37 from the universal curve");

    r_lmfdb := (xC^2*yC + xC*yC + yC - 1)/(xC^2 + xC + yC - 1);
    s_lmfdb := (xC*yC + yC - 1)/(xC + yC - 1);
    a1_lmfdb := s_lmfdb - r_lmfdb*s_lmfdb + 1;
    a2_lmfdb := r_lmfdb*s_lmfdb - r_lmfdb^2*s_lmfdb;
    a3_lmfdb := r_lmfdb*s_lmfdb - r_lmfdb^2*s_lmfdb;

    Euniv := EllipticCurve([a1_lmfdb, a2_lmfdb, a3_lmfdb, K!0, K!0]);
    P37 := Euniv![K!0, K!0, K!1];
    CheckOrDie(P37 ne Euniv!0, "The marked point P=(0,0) is the identity.");
    CheckOrDie(37*P37 eq Euniv!0, "The marked point is not 37-torsion.");
    for i in [1..18] do
        CheckOrDie(i*P37 ne Euniv!0, Sprintf("The marked point has order dividing %o.", i));
    end for;

    R<T> := PolynomialRing(K);
    ker37 := R!1;
    for i in [1..18] do
        ker37 *:= T - (i*P37)[1];
    end for;
    Equot, isog37 := IsogenyFromKernel(Euniv, ker37);

    CheckOrDie(S37 eq jInvariant(Euniv) + jInvariant(Equot),
               "Hard-coded S37 does not match reconstructed j-sum.");
    CheckOrDie(T37 eq jInvariant(Euniv) * jInvariant(Equot),
               "Hard-coded T37 does not match reconstructed j-product.");
    print "MAP_VERIFICATION: SUCCESS";
else
    print "Map reconstruction verification skipped.";
end if;

E37a1 := EllipticCurve([F | 0, 0, 1, -1, 0]);
CheckOrDie(#Points(E37a1) eq 5, "Expected #37.a1(F_2)=5.");

///////////////////////////////////////////////////////////////////////////
// 3. Recover enough verified quotient fibers to get D0.
///////////////////////////////////////////////////////////////////////////

PrintSection("3. Recover quotient fibers from [S37:T37:1]");

expected_quotient_map_degree := 36;
quotient_fibers := [];
quotient_fiber_labels := [];
fiber_seen := AssociativeArray();

function ZeroDivisorOfFunction(f)
    if f eq K!0 then
        return DivisorGroup(K)!0;
    end if;
    return PositivePartOfDivisor(Divisor(f), K);
end function;

function CommonZeroSupport(u, v)
    supp := [];
    Zv := ZeroDivisorOfFunction(v);
    for P0 in Support(ZeroDivisorOfFunction(u)) do
        if Valuation(Zv, P0) gt 0 then
            Append(~supp, P0);
        end if;
    end for;
    return supp;
end function;

function ResidueKey(P0, h)
    try
        kP, red := ResidueClassField(P0);
        return Sprintf("kdeg=%o,res=%o", Degree(kP), red(h));
    catch err
        return "residue_unavailable";
    end try;
end function;

function BranchDivisorsAtLocalPoint(u, v, label)
    groups := AssociativeArray();
    for P0 in CommonZeroSupport(u, v) do
        vu := Valuation(u, P0);
        vv := Valuation(v, P0);
        e := GCD(vu, vv);
        if e gt 0 then
            m := vu div e;
            n := vv div e;
            key := Sprintf("%o | vpair=(%o,%o), residue(u^%o/v^%o)=%o",
                           label, m, n, n, m, ResidueKey(P0, u^n/v^m));
            if not IsDefined(groups, key) then
                groups[key] := DivisorGroup(K)!0;
            end if;
            groups[key] +:= e*P0;
        end if;
    end for;

    keys := Sort([key : key in Keys(groups)]);
    return [groups[key] : key in keys], keys;
end function;

procedure AddFiber(~out_divs, ~out_labels, ~seen, D, label)
    if Degree(D) ne expected_quotient_map_degree then
        return;
    end if;

    key := Sprint(D);
    if IsDefined(seen, key) then
        return;
    end if;
    for Dold in out_divs do
        if Dold eq D then
            return;
        end if;
    end for;

    seen[key] := true;
    Append(~out_divs, D);
    Append(~out_labels, label);
end procedure;

procedure AddLocalFiberCandidates(~out_divs, ~out_labels, ~seen, u, v, label)
    divs, labs := BranchDivisorsAtLocalPoint(u, v, label);
    totalD := DivisorGroup(K)!0;
    for D in divs do
        totalD +:= D;
    end for;

    if not QUIET_SETUP then
        printf "plane point %o: %o branch candidate(s), total branch-degree %o\n",
               label, #divs, Degree(totalD);
    end if;

    for i in [1..#divs] do
        if not QUIET_SETUP then
            printf "    branch %o degree %o: %o\n", i, Degree(divs[i]), labs[i];
        end if;
        AddFiber(~out_divs, ~out_labels, ~seen, divs[i],
            "degree-36 branch: " cat labs[i]);
    end for;

    AddFiber(~out_divs, ~out_labels, ~seen, totalD,
        "degree-36 total plane fiber: " cat label);

    for i in [1..#divs] do
        if 2*Degree(divs[i]) eq expected_quotient_map_degree then
            AddFiber(~out_divs, ~out_labels, ~seen, 2*divs[i],
                "twice degree-18 branch: " cat labs[i]);
        end if;
    end for;
end procedure;

F_elements := [a : a in F];

if F2_CACHED_SCAN_TARGETS_ONLY and CACHED_GENERATOR_I eq 1 and CACHED_GENERATOR_J eq 4 then
    print "F2 cached scan: using the three target points needed for the verified <1,4> pair.";
    AddLocalFiberCandidates(~quotient_fibers, ~quotient_fiber_labels, ~fiber_seen,
        S37, T37, "[0:0:1]");
    AddLocalFiberCandidates(~quotient_fibers, ~quotient_fiber_labels, ~fiber_seen,
        S37, T37 - K!1, "[0:1:1]");
    AddLocalFiberCandidates(~quotient_fibers, ~quotient_fiber_labels, ~fiber_seen,
        1/T37, S37/T37, "[0:1:0]");
else
    for a in F_elements do
        for b in F_elements do
            AddLocalFiberCandidates(~quotient_fibers, ~quotient_fiber_labels, ~fiber_seen,
                S37 - K!a, T37 - K!b, Sprintf("[%o:%o:1]", a, b));
        end for;
    end for;
    for lam in F_elements do
        AddLocalFiberCandidates(~quotient_fibers, ~quotient_fiber_labels, ~fiber_seen,
            1/S37, T37/S37 - K!lam, Sprintf("[1:%o:0]", lam));
    end for;
    AddLocalFiberCandidates(~quotient_fibers, ~quotient_fiber_labels, ~fiber_seen,
        1/T37, S37/T37, "[0:1:0]");
end if;

printf "#degree-36 candidate fibers retained = %o\n", #quotient_fibers;
for i in [1..#quotient_fibers] do
    if not QUIET_SETUP then
        printf "    candidate fiber %o degree %o: %o\n",
               i, Degree(quotient_fibers[i]), quotient_fiber_labels[i];
    end if;
end for;
CheckOrDie(#quotient_fibers ge 2, "Fewer than two degree-36 candidate fibers were found.");

///////////////////////////////////////////////////////////////////////////
// 4. Verify D0 and build D0,2D0,3D0,4D0.
///////////////////////////////////////////////////////////////////////////

PrintSection("4. Verify an order-5 quotient delta");

function VerifyPicardOrder5(D)
    if Degree(D) ne 0 then
        return false, Sprintf("degree is %o, not 0", Degree(D));
    end if;
    if not IsPrincipalDegreeZero(5*D) then
        return false, "5*D is not principal";
    end if;
    for i in [1..4] do
        if IsPrincipalDegreeZero(i*D) then
            return false, Sprintf("%o*D is already principal", i);
        end if;
    end for;
    return true, "exact order 5 in Pic^0";
end function;

found_generator := false;
cyclic_generator_pair := <0,0>;

if CACHED_GENERATOR_I ne 0 or CACHED_GENERATOR_J ne 0 then
    CheckOrDie(CACHED_GENERATOR_I ge 1 and CACHED_GENERATOR_I le #quotient_fibers and
               CACHED_GENERATOR_J ge 1 and CACHED_GENERATOR_J le #quotient_fibers and
               CACHED_GENERATOR_I ne CACHED_GENERATOR_J,
               "Cached generator pair indices are not valid for the retained fibers.");
    D0 := quotient_fibers[CACHED_GENERATOR_J] - quotient_fibers[CACHED_GENERATOR_I];
    generator_label := Sprintf("fiber %o minus fiber %o", CACHED_GENERATOR_J, CACHED_GENERATOR_I);
    if VERIFY_DELTA_ORDER then
        ok_order, msg_order := VerifyPicardOrder5(D0);
    else
        ok_order := true;
        msg_order := "order check skipped because the runner already verified it";
    end if;
    print generator_label cat ": " cat msg_order;
    found_generator := ok_order;
    cyclic_generator_pair := <CACHED_GENERATOR_I, CACHED_GENERATOR_J>;
else
    for i in [1..#quotient_fibers] do
        if found_generator then
            break;
        end if;
        for j in [1..#quotient_fibers] do
            if i eq j then
                continue;
            end if;
            Dtry := quotient_fibers[j] - quotient_fibers[i];
            ok_order, msg_order := VerifyPicardOrder5(Dtry);
            if not QUIET_SETUP then
                printf "fiber %o minus fiber %o: %o\n", j, i, msg_order;
            end if;
            if ok_order then
                D0 := Dtry;
                generator_label := Sprintf("fiber %o minus fiber %o", j, i);
                cyclic_generator_pair := <i,j>;
                found_generator := true;
                break;
            end if;
        end for;
    end for;
end if;

CheckOrDie(found_generator, "No verified order-5 quotient delta was found.");

cyclic_deltas := [i*D0 : i in [1..4]];
cyclic_delta_labels := [Sprintf("%o*(%o)", i, generator_label) : i in [1..4]];
print Sprintf("Using candidate fiber pair <%o,%o> as a verified order-5 generator over F_2.",
              cyclic_generator_pair[1], cyclic_generator_pair[2]);

///////////////////////////////////////////////////////////////////////////
// 5. X0-style translate test.
///////////////////////////////////////////////////////////////////////////

if RUN_TRANSLATE_TEST then
    PrintSection("5. Proposition 4.9-style translate test for W^0_9");
    ok, survivors := CheckNoTranslateFromDeltasFF(C, 9, cyclic_deltas,
        cyclic_delta_labels, expected_quotient_map_degree);
else
    PrintSection("5. Translate test skipped");
    print "Set RUN_TRANSLATE_TEST:=true to enumerate effective degree-9 divisors.";
end if;

///////////////////////////////////////////////////////////////////////////
// Summary.
///////////////////////////////////////////////////////////////////////////

PrintSection("Summary");
printf "RUN_TRANSLATE_TEST        = %o\n", RUN_TRANSLATE_TEST;
printf "VERIFY_MAP_BY_RECONSTRUCTION = %o\n", VERIFY_MAP_BY_RECONSTRUCTION;
printf "VERIFY_DELTA_ORDER        = %o\n", VERIFY_DELTA_ORDER;
printf "VERIFY_MODEL_BASIC        = %o\n", VERIFY_MODEL_BASIC;
printf "NCHUNKS                   = %o\n", NCHUNKS;
printf "CHUNK_INDEX               = %o\n", CHUNK_INDEX;
printf "CACHED_GENERATOR_PAIR     = <%o,%o>\n", CACHED_GENERATOR_I, CACHED_GENERATOR_J;
printf "#quotient_fibers          = %o\n", #quotient_fibers;
printf "quotient degrees          = %o\n", [Degree(D) : D in quotient_fibers];
printf "cyclic_delta_ok           = %o\n", found_generator;
if assigned ok then
    printf "translate_test_ok         = %o\n", ok;
end if;
if assigned survivors then
    printf "#survivors                = %o\n", #survivors;
end if;

print "Finished X1_37_prop49_F2_short_complete.m.";
exit;
