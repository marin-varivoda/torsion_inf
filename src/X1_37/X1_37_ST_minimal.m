/*
    Minimal ST-relation proof for X_1(37) over F_2, degree 9.

    verify.m checks the map and the order-5 divisor once.
    worker.m runs one chunk of the effective-divisor search.
*/

function ToInt(a)
    return Type(a) eq RngIntElt select a else StringToInteger(a);
end function;

if not assigned verify then verify := true; end if;
if not assigned do_search then do_search := false; end if;
if not assigned NCHUNKS then NCHUNKS := 1; else NCHUNKS := ToInt(NCHUNKS); end if;
if not assigned CHUNK_INDEX then CHUNK_INDEX := 0; else CHUNK_INDEX := ToInt(CHUNK_INDEX); end if;
assert NCHUNKS ge 1;
assert 0 le CHUNK_INDEX and CHUNK_INDEX lt NCHUNKS;

load "ff_search.m";

print "Build X_1(37) over F_2.";
F := GF(2);
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
if verify then
    assert Genus(K) eq 40;
    kK, kK_to_K := ExactConstantField(K);
    assert #kK eq 2;
end if;

function FromCoeffs(a)
    return #a eq 0 select K!0 else &+[K!a[i+1]*yC^i : i in [0..#a-1]];
end function;

load "X1_37_X0plus_map_data.m";
S37 := FromCoeffs(S37_hard_coeffs);
T37 := FromCoeffs(T37_hard_coeffs);
assert S37 ne K!0 and T37 ne K!0;
printf "#S37 coefficients = %o; #T37 coefficients = %o\n", #S37_hard_coeffs, #T37_hard_coeffs;

if verify then
    print "Verify S37,T37 from the universal elliptic curve.";
    r := (xC^2*yC + xC*yC + yC - 1)/(xC^2 + xC + yC - 1);
    s := (xC*yC + yC - 1)/(xC + yC - 1);
    E := EllipticCurve([s - r*s + 1, r*s - r^2*s, r*s - r^2*s, K!0, K!0]);
    P := E![K!0, K!0, K!1];
    assert P ne E!0;
    assert 37*P eq E!0;
    for i in [1..18] do assert i*P ne E!0; end for;
    R<z> := PolynomialRing(K);
    ker := R!1;
    for i in [1..18] do ker *:= z - (i*P)[1]; end for;
    E2, phi := IsogenyFromKernel(E, ker);
    assert S37 eq jInvariant(E) + jInvariant(E2);
    assert T37 eq jInvariant(E)*jInvariant(E2);
    print "MAP_VERIFICATION: SUCCESS";
end if;

print "Use the quotient image q(S,T)=0.";
q_at_map := S37^38 + S37^32*T37^5 + S37^24*T37^13 + S37^14*T37^20
            + S37^8*T37^29 + S37^8*T37^25 + S37^2*T37^34
            + T37^37 + T37^33;
assert q_at_map eq K!0;
print "ST_RELATION_VERIFICATION: SUCCESS";

FS<S> := FunctionField(F);
PT<T> := PolynomialRing(FS);
q := S^38 + S^32*T^5 + S^24*T^13 + S^14*T^20
     + S^8*T^29 + S^8*T^25 + S^2*T^34 + T^37 + T^33;
Y<TY> := FunctionField(q);
SY := Y!S;
Yplaces := Places(Y, 1);
if verify then
    assert Genus(Y) eq 1;
    assert #Yplaces eq 5;
    E37a1 := EllipticCurve([F | 0, 0, 1, -1, 0]);
    assert #Yplaces eq #Points(E37a1);
end if;
printf "#Places(q(S,T)=0,1) = %o\n", #Yplaces;

function EvalFS(a)
    return Evaluate(Numerator(a), S37)/Evaluate(Denominator(a), S37);
end function;

function EvalY(h)
    c := Eltseq(h);
    return #c eq 0 select K!0 else &+[EvalFS(c[i])*T37^(i-1) : i in [1..#c]];
end function;

function Zero(f)
    return f eq K!0 select DivisorGroup(K)!0 else Pos(Divisor(f), K);
end function;

function GcdDiv(A, B)
    G := DivisorGroup(K)!0;
    for P in Setseq(Seqset(Support(A) cat Support(B))) do
        v := Min(Valuation(A, P), Valuation(B, P));
        if v gt 0 then G +:= v*P; end if;
    end for;
    return G;
end function;

function Pull(g1, g2)
    return GcdDiv(Zero(EvalY(g1)), Zero(EvalY(g2)));
end function;

function Hits(g1, g2)
    return [<i, Valuation(g1, Yplaces[i]), Valuation(g2, Yplaces[i])> :
            i in [1..#Yplaces] | Valuation(g1, Yplaces[i]) gt 0 and Valuation(g2, Yplaces[i]) gt 0];
end function;

g_inf := <1/SY, TY^36/SY^37>;
g_branch := <SY, TY/SY + Y!1>;
for g in [g_inf, g_branch] do
    h := Hits(g[1], g[2]);
    assert #h eq 1;
    assert Min(h[1][2], h[1][3]) eq 1;
end for;

Dinf := Pull(g_inf[1], g_inf[2]);
Dbranch := Pull(g_branch[1], g_branch[2]);
printf "fiber degrees = %o\n", [Degree(Dinf), Degree(Dbranch)];
assert Degree(Dinf) eq 36 and Degree(Dbranch) eq 36;
assert Dinf ne Dbranch;

D0 := Dinf - Dbranch;
function Principal(D)
    return Degree(D) eq 0 and Dimension(RiemannRochSpace(D)) gt 0;
end function;

if verify then
    assert Principal(5*D0);
    for i in [1..4] do assert not Principal(i*D0); end for;
    print "DELTA_ORDER: SUCCESS";
end if;

if do_search then
    deltas := [i*D0 : i in [1..4]];
    labels := [Sprintf("%o*D0", i) : i in [1..4]];
    ok, survivors := NoTranslate(K, 9, deltas, labels);
end if;

print "Finished X1_37_ST_minimal.m.";
exit;
