/* 
    The following is the list of all congruence subgroups of SL(2,Z), up to conjugacy in GL(2,Z), that contains -I such that
    the corresponding modular curve is isomorphic to a smooth plane quintic.  These groups are given by their Cummins-Pauli label.
        
    This code check that the list indeed has these properties!
*/

smooth_plane_quintics:=["75A6","75D6"];


// We have already computed the congruence subgroups of gonality 2; this is the list with genus > 2
gonality2:=[ "8B3", "10B3", "12C3", "12D3", "12E3", "12F3", "12G3", "12H3", "12K3", 
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

// We have already computed the congruence subgroups of gonality at 3; this is the list with genus > 4

gonality3:=[ "54C5", "16A6", "18A6", "18D6", "24D6", "27A6", "28D6", "28E6", 
    "30C6", "32A6", "36C6", "36H6", "36J6", "36K6", "39A6", "45D6", "54A6", "54B6", "56D6", 
    "64A6", "84A6", "108A6", "27B7", "27C7", "30D7", "42M7", "24A8", "24B8", "36H8", "36I8", 
    "36J8", "36K8", "48A8", "48C8", "48E8", "72F8", "72G8", "84A8", "96A8", "108A8", "108B8", 
    "144A8", "15A10", "36A10", "36C10", "42G10", "72A10", "75A10", "108A10", "108C10", "108A12"];


// loads all the needed modular curve functions
load "classification_functions.m";

// Cummins-Pauli data of congruence subgroups of genus at most 24
cp_data:=CumminsPauliData();

smooth_plane_quintics_found:=[];

for r in cp_data do
    if r`genus ne 6 or r`name in gonality2 cat gonality3 then continue r; end if;

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

    // Create the modular curve X_G
    X:=CreateModularCurveRec(G); 

    b1,b2:=HasGonalityThree(X);
    assert not b1;
    if b2 then  
        smooth_plane_quintics_found cat:= [r`name];
    end if;
end for;

assert Set(smooth_plane_quintics_found) eq Set(smooth_plane_quintics);
