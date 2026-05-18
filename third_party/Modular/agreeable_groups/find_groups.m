//  nohup magma "find_groups.m" &

function HenselLift(F,pts,p,e : max:=1)
    /*
        Input:
            F : nonempty sequence of polynomials f in Z[x_1,..,x_r] with r>0
            pts:  sequence of a in Z^r such that f(a) is congruent to 0 modulo p^e for all f in F
            p: prime
            e: positive integer

        Output:
            A sequence of points in Z^r that give representatives of the solutions of the equations
            f(x)=0 mod p^(e+1) so that x is congruent modulo p^e to an element of pts.

        Function recurses until solutions are found modulo at least p^max.
    */
    
    if #pts eq 0 then return []; end if;

    r:=#pts[1];
    Pol<[x]>:=PolynomialRing(Integers(),r);

    pts_new:=[];
    for P in pts do
        c:=[Evaluate(f,P):f in F];
        assert &and[a mod p^e eq 0 : a in c]; //check
        c:=Matrix([[a div p^e : a in c]]);
        c:=ChangeRing(c,GF(p));

        b:=[P[i]+x[i]: i in [1..r]];
        A:=Matrix([ [MonomialCoefficient(Evaluate(f,b),x[i]): i in [1..r]]: f in F]);
    
        A:=ChangeRing(A,GF(p));
        flag,v0,N:=IsConsistent(Transpose(A),-c);
        if not flag then
            return [];
        end if;

        for n in N do
            v:=Vector(v0)+n;
            v:=[Integers()!a: a in Eltseq(v)];
            P1:=[P[i]+p^e*v[i]: i in [1..r]];
            pts_new cat:= [P1];
            assert &and [Evaluate(f,P1) mod p^(e+1) eq 0: f in F]; //check
        end for;
    end for;

    if pts_new eq [] then return []; end if;
    if e lt max then
        return HenselLift(F,pts_new,p,e+1: max:=max); //recurse
    end if;

    return pts_new;
end function;


function HasLocalPoints(M,p,n :dim_bound:=Infinity())
    /*
        M:  a modular curve with a projective model defined by polynomials with integer ccoefficients, using these models view as scheme over Z
        p:  a prime
        n:  a positive integer

        Returns true if and only if M(Z/p^n) is nonempty or M base changed to F_p has dimension at least "dim_bound".

        [In particular if false is returned, then M has no Q_p-points]
    */
    if M`psi eq [] then return true; end if; // Curve is P1

    r:=Rank(Parent(M`psi[1]));
    
    R<[x]>:=PolynomialRing(Integers(),r);
    psi:=[R!f: f in M`psi];

    // Use built in Magma functions to find points of M(Z/pZ); this can be slow
    PP:=ProjectiveSpace(GF(p),r-1);
    C:=Scheme(PP,psi);
    dimC:=Dimension(C);
    if dimC ge dim_bound then
        // if C has dimension too large can give up; may be too slow to find F_p-points.
        return true;
    end if;


    ptsC:=Points(C);
    if #ptsC eq 0 then return false; end if;
    if n eq 1 then return true; end if;

    R_<[y]> :=PolynomialRing(Integers(),r-1);
    
    // Run over each point and try to lift them
    for P0 in  ptsC do

        // Dehomogenize before trying to lift
        P1:=Eltseq(P0);
        j:=[i: i in [1..#P1] | P1[i] ne 0][1];
        P1:=[P1[i]/P1[j]: i in [1..#P1] | i ne j];
        P1:=[Integers()!a: a in P1];
        b:=[y[i]: i in [1..j-1]] cat [1] cat [y[i]: i in [j..r-1]];
        psi1:=[Evaluate(f,b): f in psi];

        lifts:=HenselLift(psi1,[P1],p,1 : max:=n);
        if #lifts ne 0 then
            return true;
        end if;
    end for;
      
    return false;
end function;



" ";
"Loading groups from `OpenImage` and initializing modular curves.";

AttachSpec("../Modular.spec");

// We load the groups G from the "OpenImage" repository and initialize a modular curve.
I:=Open("agreeable.dat", "r"); 
X:=AssociativeArray();
repeat
	b,y:=ReadObjectCheck(I);
	if b then
		X[y`key]:=CreateModularCurveRec(y`level,y`gens);
        if assigned y`is_agreeable then
            X[y`key]`is_agreeable:=y`is_agreeable;
        end if;
        if assigned y`has_infinitely_many_points then
            X[y`key]`has_infinitely_many_points:=y`has_infinitely_many_points;
        end if;   
        if assigned y`is_entangled then
            X[y`key]`is_unentangled:= not y`is_entangled;
        end if;                 
	end if;
until not b;

" ";
"Finding those entangled groups of genus at least 2 from `OpenImage' that have no obvious obstruction to rational points";
// For those modular curves of genus at least 2, we look for easy obstructions to rational points.
// They were constructed so that they have real points.
S1:=[];
for k in Keys(X) do
    if X[k]`genus le 1 then continue k; end if;
    
    M:=X[k];
    M:=FindModelOfXG(M); // compute a model

    n:=3;
    obstruction:=exists{p: p in PrimeDivisors(2*3*M`N) | HasLocalPoints(M,p,n) eq false};
    if obstruction eq false then
        S1:=S1 cat [k]; // keep track of curves that might have a Q-point.
    end if;
end for;


" ";
"Finding those unentangled groups of genus at least 2 that are implict in `OpenImage'";

/*  
    We now find agreeable groups, up to conjugacy, that are unentangled and whose ell-adic 
    projections give modular curves with infinitely many Q-points.   We now find all such groups G of 
    genus at least 2, up to conjugacy, that are maximal with respect to the property of X_G(Q) being finite.
*/
keys:=[k: k in Keys(X) | X[k]`genus le 1 and X[k]`is_agreeable and X[k]`has_infinitely_many_points and X[k]`is_unentangled];

pairs:=[];
S:={**};
count:=0;
for k1 in keys do
    count:=count+1;
    //[count,#keys];
    for k2 in keys do
        N1:=X[k1]`N;
        N2:=X[k2]`N;

        // Want N1 and N2 to be >1, relatively prime, with N2 a power of a prime ell that larger than any prime dividing N1
        if GCD(N1,N2) ne 1 then continue k2; end if;
        if N1 eq 1 or N2 eq 1 then continue k2; end if;
        if IsPrimePower(N2) eq false then continue k2; end if;
        _,ell,_:=IsPrimePower(N2);
        if &and[p lt ell : p in PrimeDivisors(N1)] eq false then continue k2; end if;

        N:=N1*N2;
        G1:=GL2Lift(X[k1]`G,N);
        G2:=GL2Lift(X[k2]`G,N);
        G:=G1 meet G2;

        if GL2Genus(G) le 1 then
            // We already know the agreeable groups of genus at most 1
            continue k2;
        end if;

        // We now compute some agreeable groups W that contain those for which
        // G is a maximal agreeable subgroup.   If any of these have genus at least
        // 2 or if they have genus at most 1 and X_W(Q) is finite, then we can
        // ignore this G.
        minimal_overgroups:=MinimalOvergroups(GL2Ambient(N),G);
        agreeable_closures:=[];
        for W in minimal_overgroups do
            W_:=GL2AgreeableClosure(W);
            g_:=GL2Genus(W_);
            if g_ ge 2 then 
                // In this case G is not maximal amongst agreeable groups of genus at least 2
                continue k2;
            end if;
            agreeable_closures cat:= [W_];
        end for;
        for W in agreeable_closures do
            genus:=GL2Genus(W);
            assert genus in {0,1};
            level:=GL2Level(W);
            index:=GL2Index(W);
            GL2:=GL2Ambient(level);
            for k in Keys(X) do
                if X[k]`genus eq genus and X[k]`N eq level and X[k]`index eq index then
                    if IsConjugate(GL2,X[k]`G,W) then
                        continue W;
                    end if;
                end if;
            end for;

            // W is agreeable of genus at most 1 and not in the list X.            
            continue k2;  // have already dealt with these cases
        end for;

        pairs:=pairs cat [[k1,k2]];
        //#pairs;

        S:=S join {*GL2Genus(G)*}; 
        //S;
    end for;
end for;

" ";
"Looking for obvious obstructions to rational points";
// We now look for some obvious local obstructions to rational points
S2:=[];
count:=0;
for u in pairs do
    count:=count+1;

    k1:=u[1];
    k2:=u[2];
    N1:=X[k1]`N;
    N2:=X[k2]`N;
    N:=N1*N2;
    assert GCD(N1,N2) eq 1;
    G1:=GL2Lift(X[k1]`G,N);
    G2:=GL2Lift(X[k2]`G,N);
    G:=G1 meet G2;
    
    g:=GL2Genus(G);
 
    if GL2Genus(G) le 40 then
        // restrict to genus at least 40 for now
        M:=CreateModularCurveRec(G);
        M:=FindModelOfXG(M);
 
        n:=2; 
        dim_bound:=2; // choices
        obstruction:=exists{p: p in PrimeDivisors(2*M`N) | HasLocalPoints(M,p,n: dim_bound:=dim_bound) eq false};

        if obstruction eq false then
            S2:=S2 cat [u];
        end if;
    else
        S2:=S2 cat [u];
    end if;    
    //[#S2,count];
end for;

" ";
"Optimizing groups and putting in a sequence ready for output";

ExtraGroupRec := recformat<N:RngIntElt, gens,label:SeqEnum, LMFDBlabel:MonStgElt, is_unentangled:BoolElt>;

groups:=[];
for k in S1 do
    N:=X[k]`N;

    G:=X[k]`G;
    G:=GL2MinimizeGenerators(G);  // choose smaller set of generators
    gens:=[[Integers()!a : a in Eltseq(g)] : g in Generators(G)];

    r:=rec<ExtraGroupRec  | N:=N, gens:=gens, label:=k, is_unentangled:=false>;
    groups cat:= [r];
end for;

for u in S2 do
    k1:=u[1];
    k2:=u[2];
    N1:=X[k1]`N;
    N2:=X[k2]`N;
    N:=N1*N2;
    
    G1:=GL2Lift(X[k1]`G,N);
    G2:=GL2Lift(X[k2]`G,N);
    G:=G1 meet G2;  
    assert GL2Level(G) eq N1*N2;
    G:=GL2MinimizeGenerators(G);
    gens:=[[Integers()!a : a in Eltseq(g)] : g in Generators(G)];

    label:=["1A0-1a": i in [1..9]];
    for i in [1..9] do
        if k1[i] ne "1A0-1a" then label[i]:=k1[i]; end if;
        if k2[i] ne "1A0-1a" then label[i]:=k2[i]; end if;
    end for;

    r:=rec<ExtraGroupRec  | N:=N, gens:=gens, label:=label, is_unentangled:=true>;
    groups cat:= [r];      
end for;

// Outputs list of groups to a file.  This was copy-pasted into "groups.m".
PrintFileMagma("temp.m", groups);
" ";
"Done.";