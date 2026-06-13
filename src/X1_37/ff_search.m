/* Generic effective-divisor search for the X1(37) proof. */

function Pos(D, C)
    E := DivisorGroup(C)!0;
    for P in Support(D) do
        v := Valuation(D, P);
        if v gt 0 then E +:= v*P; end if;
    end for;
    return E;
end function;

function DegreeTypes(d, places)
    ds := [i : i in [1..d] | #places[i] gt 0];
    old := [ <[<0,0>], d> ];
    new := old;
    done := [];
    for deg in Reverse(ds) do
        for m in Reverse([1..Floor(d/deg)]) do
            for t in old do
                for r in [1..Floor(t[2]/(m*deg))] do
                    Append(~new, <t[1] cat [<m,deg> : j in [1..r]], t[2] - m*deg*r>);
                end for;
            end for;
            done cat:= [t[1][2..#t[1]] : t in new | t[2] eq 0];
            new := [t : t in new | t[2] gt 0];
            old := new;
        end for;
    end for;
    seen := AssociativeArray(); keys := [];
    for t in done do
        key := Sprint(Sort([Sprint(x) : x in t]));
        if not IsDefined(seen, key) then seen[key] := t; Append(~keys, key); end if;
    end for;
    return [seen[k] : k in Sort(keys)];
end function;

function Groups(t)
    g := [];
    for a in t do
        found := false;
        for i in [1..#g] do
            if g[i][1] eq a[1] and g[i][2] eq a[2] then
                g[i] := <a[1], a[2], g[i][3] + 1>; found := true; break;
            end if;
        end for;
        if not found then Append(~g, <a[1], a[2], 1>); end if;
    end for;
    return g;
end function;

function ResidueCoords(c, F)
    if Type(c) eq RngIntElt then return [F!c]; end if;
    if Parent(c) cmpeq F then return [F!c]; end if;
    ok, b := IsCoercible(F, c);
    if ok then return [b]; end if;
    v := [F | ];
    for a in Eltseq(c) do v cat:= ResidueCoords(a, F); end for;
    return v;
end function;

function LocalCoeffs(h, P, n)
    u := LocalUniformizer(P); rem := h; coeffs := [];
    for i in [1..n] do
        c := Evaluate(rem, P);
        if c cmpeq Infinity() then error "non-regular local expansion"; end if;
        Append(~coeffs, c); rem := (rem - Lift(c, P))/u;
    end for;
    return coeffs;
end function;

function RowsAt(fs, P, shift, n, F)
    rows := [];
    cs := [LocalCoeffs(f*LocalUniformizer(P)^shift, P, n) : f in fs];
    for r in [1..n] do
        for a in [1..Degree(P)] do
            row := [F | ];
            for j in [1..#fs] do
                v := ResidueCoords(cs[j][r], F);
                Append(~row, a le #v select v[a] else F!0);
            end for;
            Append(~rows, row);
        end for;
    end for;
    return rows;
end function;

function CoverDivisor(t, C, places)
    max := AssociativeArray();
    for a in t do
        if not IsDefined(max, a[2]) or a[1] gt max[a[2]] then max[a[2]] := a[1]; end if;
    end for;
    H := DivisorGroup(C)!0;
    for deg in Keys(max) do for P in places[deg] do H +:= max[deg]*P; end for; end for;
    return H;
end function;

function FastData(C, H, delta)
    G := H + Pos(delta, C);
    R := H + Pos(-delta, C);
    t := Cputime(); V, m := RiemannRochSpace(G); rr := Cputime(t);
    fs := [m(b) : b in Basis(V)];
    ps := []; ns := []; shifts := [];
    for P in Support(R) do
        n := Valuation(R, P);
        if n gt 0 then Append(~ps, P); Append(~ns, n); Append(~shifts, Valuation(G, P)); end if;
    end for;
    return <#fs, ps, ns, fs, shifts, AssociativeArray(), Dimension(V), rr>;
end function;

function ReplaceCache(data, cache)
    return <data[1], data[2], data[3], data[4], data[5], cache, data[7], data[8]>;
end function;

function CachedRows(data, i, F)
    cache := data[6];
    if not IsDefined(cache, i) then cache[i] := RowsAt(data[4], data[2][i], data[5][i], data[3][i], F); end if;
    return ReplaceCache(data, cache), cache[i];
end function;

function FastInW0(D, data, F)
    ncols := data[1];
    if ncols eq 0 then return false, data; end if;
    jobs := [];
    for i in [1..#data[2]] do
        keep := data[3][i] - Valuation(D, data[2][i]);
        if keep lt 0 then error "cover divisor missed D"; end if;
        if keep gt 0 then Append(~jobs, <keep*Degree(data[2][i]), i>); end if;
    end for;
    if #jobs eq 0 then return true, data; end if;
    rows := []; nrows := 0;
    for job in Reverse(Sort(jobs)) do
        data, r := CachedRows(data, job[2], F);
        rows cat:= r[1..job[1]]; nrows +:= job[1];
        if nrows ge ncols and Rank(Matrix(F, #rows, ncols, &cat rows)) eq ncols then
            return false, data;
        end if;
    end for;
    return Rank(Matrix(F, #rows, ncols, &cat rows)) lt ncols, data;
end function;

procedure TestD(~survivors, ~stop, ~checked, D, deltas, labels, ~data, F, t0)
    checked +:= 1; bad := 0;
    for j in [1..#deltas] do
        ok, dj := FastInW0(D, data[j], F); data[j] := dj;
        if not ok then bad := j; break; end if;
    end for;
    if bad eq 0 then
        Append(~survivors, D);
        printf "SURVIVOR %o in chunk %o after %o tested divisors: %o\n", #survivors, CHUNK_INDEX, checked, D;
        print "INCONCLUSIVE: stopped after first survivor."; print "CHUNK_RESULT: INCONCLUSIVE"; stop := true;
    elif checked mod 1000 eq 0 then
        printf "Chunk %o checked %o divisors; last killed by %o; survivors %o; elapsed %o\n",
               CHUNK_INDEX, checked, labels[bad], #survivors, Cputime(t0);
    end if;
end procedure;

procedure Stream(~survivors, ~stop, ~checked, groups, i, chosen, nextp, D, used, places, deltas, labels, ~data, F, t0)
    if stop then return; end if;
    if i gt #groups then TestD(~survivors, ~stop, ~checked, D, deltas, labels, ~data, F, t0); return; end if;
    mult := groups[i][1]; deg := groups[i][2]; need := groups[i][3];
    if chosen eq need then
        Stream(~survivors, ~stop, ~checked, groups, i+1, 0, 1, D, used, places, deltas, labels, ~data, F, t0); return;
    end if;
    last := #places[deg] - (need - chosen) + 1;
    for k in [nextp..last] do
        P := places[deg][k];
        if P notin used then
            Stream(~survivors, ~stop, ~checked, groups, i, chosen+1, k+1, D + mult*P,
                   used cat [P], places, deltas, labels, ~data, F, t0);
            if stop then break; end if;
        end if;
    end for;
end procedure;

function NoTranslate(C, d, deltas, labels)
    printf "Testing %o nonzero translates; chunk %o/%o.\n", #deltas, CHUNK_INDEX, NCHUNKS;
    places := AssociativeArray();
    for i in [1..d] do places[i] := Places(C, i); end for;
    types := DegreeTypes(d, places);
    mine := [i : i in [1..#types] | (i - 1) mod NCHUNKS eq CHUNK_INDEX];
    printf "Place counts by degree: %o\n", [#places[i] : i in [1..d]];
    printf "Degree types: %o; assigned %o\n", #types, mine;
    survivors := []; checked := 0; t0 := Cputime(); F := GF(2);
    for n in mine do
        g := Groups(types[n]); H := CoverDivisor(types[n], C, places);
        printf "Degree type %o/%o: %o; groups %o; elapsed %o\n", n, #types, types[n], g, Cputime(t0);
        data := [* *]; dims := []; times := [];
        for delta in deltas do
            D := FastData(C, H, delta); Append(~data, D); Append(~dims, D[7]); Append(~times, D[8]);
        end for;
        printf "Fast W0 setup: cover degree %o; dims %o; RR times %o\n", Degree(H), dims, times;
        stop := false;
        Stream(~survivors, ~stop, ~checked, g, 1, 0, 1, DivisorGroup(C)!0, [], places, deltas, labels, ~data, F, t0);
        if stop then return false, survivors; end if;
    end for;
    printf "Divisors tested in chunk %o: %o\n", CHUNK_INDEX, checked;
    printf "Surviving translates in this chunk: %o\n", #survivors;
    if #survivors eq 0 then print "SUCCESS: this chunk contains no translate survivor."; print "CHUNK_RESULT: SUCCESS"; return true, survivors; end if;
    print "INCONCLUSIVE: at least one quotient coset survived in this chunk."; print "CHUNK_RESULT: INCONCLUSIVE";
    return false, survivors;
end function;
