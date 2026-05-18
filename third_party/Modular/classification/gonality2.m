/*  
    The following is the list of all congruence subgroups of SL(2,Z), up to conjugacy in GL(2,Z),
    that contain -I, that have genus at least 3 and have gonality 2.  They are given by their 
    Cummins-Pauli label.  (Congruence subgroups of genus 1 and 2 always have gonality 2)

    This code checks that the list indeed has these properties!
*/
gonality_equals_2:=[ "8B3", "10B3", "12C3", "12D3", "12E3", "12F3", "12G3", "12H3", "12K3", 
"12L3", "14A3", "14C3", "14F3", "15F3", "15G3", "16B3", "16C3", "16D3", "16E3", "16F3", 
"16I3", "16J3", "16M3", "16S3", "18A3", "18C3", "18F3", "18G3", "20C3", "20F3", "20G3", 
"20H3", "20I3", "20J3", "20M3", "20O3", "21A3", "21B3", "21D3", "24A3", "24B3", "24C3", 
"24G3", "24I3", "24K3", "24L3", "24M3", "24S3", "24U3", "24V3", "24W3", "28C3", "28E3", 
"30B3", "30G3", "30J3", "30K3", "30L3", "32B3", "32C3", "32D3", "32H3", "32K3", "32M3", 
"33C3", "34B3", "35A3", "36E3", "36F3", "36G3", "39A3", "40D3", "40E3", "40F3", "40I3", 
"41A3", "42E3", "48C3", "48E3", "48F3", "48H3", "48I3", "48J3", "48M3", "50A3", "54A3", 
"60C3", "60D3", "64A3", "96A3", "18B4", "25A4", "25D4", "32B4", "36C4", "42A4", "44B4", 
"47A4", "48C4", "50A4", "50D4", "10A5", "14C5", "16G5", "18A5", "24A5", "24D5", "26A5", 
"30C5", "30F5", "36A5", "36B5", "36H5", "40A5", "42A5", "44B5", "45A5", "45C5", "46A5", 
"48A5", "48E5", "48F5", "48G5", "48H5", "50A5", "50D5", "50F5", "52B5", "54A5", "57A5", 
"58A5", "59A5", "60A5", "96A5", "48A6", "71A6", "32E7", "48N7", "56B7", "64D7", "82B7", 
"96A7", "93A8", "50A9", "50D9", "96B9", "48B11", "72A11", "96B11"];

// loads all the needed modular curve functions
load "classification_functions.m";

// Cummins-Pauli data of congruence subgroups of genus at most 24
cp_data:=CumminsPauliData();

total_time:=Realtime();
count:=0; KG_degrees:=[]; // only used for some numerics in the paper

// list of groups found so far with gonality at most 2
gonality_at_most_2:=[];

for r in cp_data do
    // We check all congruence subgroups in Cummin-Pauli data; this will include all
    // congruence subgroups of genus up to 24 (this will be enough to find all gonality 2 cases)

    if r`genus le 2 then
        // Curves of genus 0 have gonality 1.  Curves of genus 1 or 2 have gonality 2.
        gonality_at_most_2 cat:=[r`name];
        continue r;
    end if;
  
    if (r`level le 226 and r`index gt 191) or (r`level gt 226 and r`index gt 201) then
        // Our gonality bounds show that gonality is not 2
        continue r;
    end if;

    for p in [2,3,5,7] do
        if r`level mod p ne 0 and r`index gt 24*(p^2+1)/(p-1) then
            // Our gonality bounds show that gonality is not 2
            continue r;        
        end if;
    end for;

    for i in r`supergroups do
        s:=cp_data[i];
        // corresponds to a larger congruence subgroup

        // A curve of gonality 2 can only have a nonconstant map to a curve of gonality at most 2.
        // Note that groups in "cp_data" are listed in terms of increasing genus, 
        // then increasing level, then increasing index.
        if s`name notin gonality_at_most_2 then
            continue r;
        end if;

        d1:=r`index div s`index; 
        if s`genus eq 0 and d1 eq 2 then
            // The curve is a degree 2 cover of a modular curve of genus 0,
            // and hence gonality is at most 2.
            gonality_at_most_2 cat:=[r`name];
            continue r;           
        end if;

    end for;

    // We now check the Castelnuovo-Severi inequality using all strictly 
    // larger congruence subgroups.  We first find the set I of these larger 
    // groups up to conjugacy in GL(2,Z).
    I_:=Set(r`supergroups);
    repeat
        I:=I_;
        I_:=I join &join[ Set(cp_data[i]`supergroups): i in I];
    until I eq I_;

    for i in I do
        s:=cp_data[i];
        d1:=r`index div s`index;
		g1:=s`genus;
		g :=r`genus;
        if g gt d1*g1+d1-1 then
            //Castelnuovo-Severi inequality would fail if the gonality was 2
            continue r;
        end if;
    end for;

    /*  
        We now construct an open subgroup G of GL(2,Zhat) such that the intersection of G with 
        SL(2,Zhat) corresponds to the congruence subgroup.

        We find G so that the index of det(G) in Zhat^* is minimal and G has the same level
        as the congruence subgroup.
    */

    N:=r`level; 
    gens:=r`matgens;
    SL2:=SL2Ambient(N);
    GL2:=GL2Ambient(N);   
    H:=sub<SL2|gens>; 
    H`SL:=true;
    H`Index:=r`index;
    H`Genus:=r`genus;
    H`Order:=SL2Size(N) div H`Index;
 
    GG:=Normalizer(GL2,H);
    HH:=SL2Intersection(GG);
    GG`Order:=GL2Order(GG); 

    Q,iota:=quo<GG|H>;
    iotaHH:=iota(HH);
    MM:=[M`subgroup @@ iota: M in Subgroups(Q) | #(iotaHH meet M`subgroup) eq 1];

    min:=Minimum([GL2Index(M): M in MM]);
    MM:=[M: M in MM | GL2Index(M) eq min];
    G:=MM[1];
    assert SL2Intersection(G) eq H;
    assert GL2Level(G) eq N;

    r`name;
    count:=count+1; 
    KG_degrees:=Sort(KG_degrees cat [GL2DeterminantIndex(G)]);

    // Create the modular curve X_G
    X:=CreateModularCurveRec(G); 

    // Checks whether gonality is 2
    if HasGonalityTwo(X) then
        gonality_at_most_2 cat:=[r`name];
    end if;

end for;

// Check that our initial list agrees with the one just computed!
for r in cp_data do
    if r`genus ge 3 then
        assert r`name in gonality_at_most_2 eq r`name in gonality_equals_2;
    end if;
end for;

"Done";
"Total time:",Realtime(total_time);
