// This file saves the groups of our classification for future use.

load "classification_functions.m";

// Cummins-Pauli data of congruence subgroups of genus at most 24
cp_data:=CumminsPauliData();

// Congruence subgroups of gonality 1
gonality1:=[r`name: r in cp_data | r`genus eq 0];

// Congruence subgroups of gonality 2
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
gonality2:=[r`name: r in cp_data | r`genus in {1,2}] cat gonality2;

// Congruence subgroups of gonality 3
gonality3:=[ "54C5", "16A6", "18A6", "18D6", "24D6", "27A6", "28D6", "28E6", 
    "30C6", "32A6", "36C6", "36H6", "36J6", "36K6", "39A6", "45D6", "54A6", "54B6", "56D6", 
    "64A6", "84A6", "108A6", "27B7", "27C7", "30D7", "42M7", "24A8", "24B8", "36H8", "36I8", 
    "36J8", "36K8", "48A8", "48C8", "48E8", "72F8", "72G8", "84A8", "96A8", "108A8", "108B8", 
    "144A8", "15A10", "36A10", "36C10", "42G10", "72A10", "75A10", "108A10", "108C10", "108A12"];
gonality3:=[r`name: r in cp_data | r`genus in {3,4} and r`name notin gonality2] cat gonality3;

// Congruence subgroups for which the modular curve is isomorphic to a smooth plane quintic
smooth_plane_quintic:=["75A6","75D6"];

// Congruence subgroups that are bielliptic
bielliptic:=    
[ "6A1", "6B1", "6C1", "6D1", "6E1", "6F1", "7A1", "7B1", "7C1", "8A1", "8B1", "8C1", "8D1", "8E1", "8F1", "8G1", "8H1", "8I1", 
  "8J1", "8K1", "9A1", "9B1", "9C1", "9D1", "9E1", "9F1", "9G1", "9H1", "10A1", "10B1", "10C1", "10D1", "10E1", "10F1", "10G1", 
  "10H1", "10I1", "10J1", "10K1", "11A1", "11B1", "11C1", "11D1", "12A1", "12B1", "12C1", "12D1", "12E1", "12F1", "12G1", "12H1", 
  "12I1", "12J1", "12K1", "12L1", "12M1", "12N1",  "12O1", "12P1", "12Q1", "12R1", "12S1", "12T1", "12U1", "12V1", "14A1", "14B1", 
  "14C1", "14D1", "14E1", "14F1", "14G1", "14H1", "15A1", "15B1", "15C1", "15D1", "15E1", "15F1", "15G1", "15H1", "15I1", "16A1", 
  "16B1", "16C1", "16D1", "16E1", "16F1", "16G1", "16H1", "16I1", "16J1", "16K1", "16L1", "16M1", "17A1", "17B1", "17C1", "18A1", 
  "18B1", "18C1", "18D1", "18E1", "18F1", "18G1", "18H1", "18I1", "18J1", "18K1", "19A1", "19B1", "20A1", "20B1", "20C1", "20D1", 
  "20E1", "20F1", "20G1", "20H1", "20I1", "20J1", "21A1", "21B1", "21C1", "21D1", "21E1", "21F1", "22A1", "24A1", "24B1", "24C1", 
  "24D1", "24E1", "24F1", "24G1", "24H1", "24I1", "24J1", "26A1", "26B1", "27A1", "27B1", "27C1", "28A1", "30A1", "30B1", "30C1", 
  "30D1", "32A1", "32B1", "32C1", "32D1", "32E1", "33A1", "36A1", "36B1", "36C1", "39A1", "40A1", "42A1", "42B1", "49A1", "52A1", 
  "8A2", "8B2", "8C2", "9A2", "10A2", "10B2", "10C2", "10D2", "10E2", "10F2", "11A2", "12A2", "12B2", "12C2", "12D2", "12E2", "12F2", 
  "12G2", "12H2", "12I2", "13A2", "14A2", "14C2", "14D2", "14E2", "15A2", "15B2", "15D2", "16A2", "16B2", "16C2", "16D2", "16E2", 
  "16F2", "16G2", "16I2", "16J2", "16K2", "16L2", "18B2", "18C2", "18D2", "18E2", "18G2", "18H2", "18I2", "18L2", "18N2", "18O2", 
  "18P2", "18Q2", "20A2", "20B2", "20C2", "20D2", "20F2", "21A2", "21C2", "22A2", "22C2", "24A2", "24B2", "24C2", "24D2", "24E2", 
  "24F2", "24G2", "24H2", "24I2", "24J2", "24K2", "24L2", "24M2", "24N2", "24O2", "24P2", "24Q2", "26A2", "26B2", "27A2", "27B2", 
  "28A2", "28B2", "28C2", "28D2", "28E2", "28F2", "30A2", "30B2", "30C2", "30D2", "30E2", "30F2", "32A2", "32B2", "32C2", "36A2", 
  "36B2", "36C2", "36D2", "37A2", "38A2", "39A2", "40A2", "42B2", "42C2", "45A2", "48A2", "50A2", "50B2", "54A2", "54B2", "64A2", 
  "78A2", "7A3", "8A3", "8B3", "10A3", "10B3", "10C3", "10D3", "11A3", "12A3", "12B3", "12C3", "12D3", "12E3", "12F3", "12G3", "12H3", 
  "12I3", "12J3", "12K3", "12L3", "12M3", "12N3", "12O3", "12P3", "14B3", "14D3", "14E3", "14F3", "15A3", "15C3", "15D3", "15E3", 
  "15F3", "15G3", "15H3", "15I3", "16A3", "16B3", "16C3", "16D3", "16E3", "16F3", "16G3", "16H3", "16I3", "16J3", "16K3", "16L3", 
  "16M3", "16N3", "16O3", "16P3", "16Q3", "16R3", "16S3", "18A3", "18B3", "18C3", "18D3", "18E3", "18F3", "18G3", "18H3", "18I3", 
  "18J3", "18K3", "20A3", "20B3", "20C3", "20D3", "20E3", "20F3", "20G3", "20H3", "20I3", "20J3", "20K3", "20L3", "20M3", "20N3", 
  "20O3", "20P3", "20Q3", "20R3", "20S3", "20T3", "21A3", "21B3", "21D3", "24A3", "24B3", "24C3", "24D3", "24E3", "24F3", "24G3", 
  "24H3", "24I3", "24J3", "24K3", "24L3", "24M3", "24N3", "24O3", "24P3", "24Q3", "24R3", "24S3", "24T3", "24U3", "24V3", "24W3", 
  "24X3", "24Y3", "24Z3", "24AA3", "24AB3", "24AC3", "26A3", "27A3", "28A3", "28B3", "28C3", "28D3", "28E3", "30A3", "30B3", "30C3", 
  "30D3", "30E3", "30F3", "30G3", "30H3", "30I3", "30J3", "30K3", "30L3", "32A3", "32B3", "32C3", "32D3", "32E3", "32F3", "32G3", 
  "32H3", "32I3", "32J3", "32K3", "32L3", "32M3", "32N3", "32O3", "32P3", "32Q3", "33A3", "33B3", "33C3", "34A3", "34B3", "34C3", 
  "35A3", "36A3", "36B3", "36C3", "36D3", "36E3", "36F3", "36G3", "36H3", "36I3", "36J3", "36K3", "39A3", "40A3", "40B3", "40C3", 
  "40D3", "40E3", "40F3", "40G3", "40H3", "40I3", "40J3", "42A3", "42C3", "42D3", "42E3", "42F3", "43A3", "45B3", "45C3", "45D3", 
  "48A3", "48B3", "48C3", "48D3", "48E3", "48F3", "48G3", "48H3", "48I3", "48J3", "48K3", "48L3", "48M3", "49A3", "50A3", "51A3", 
  "52A3", "52B3", "54A3", "54B3", "54C3", "56A3", "56B3", "56C3", "56D3", "60A3", "60B3", "60C3", "60D3", "64A3", "64B3", "66A3", 
  "72A3", "84A3", "96A3", "9A4", "9B4", "9C4", "10A4", "10B4", "11A4", "12A4", "12C4", "14A4", "14B4", "15B4", "15C4", "15D4", "15E4", 
  "16A4", "16B4", "16C4", "18A4", "18C4", "18D4", "18E4", "18F4", "18G4", "18H4", "18I4", "18J4", "18K4", "18M4", "18N4", "18O4", 
  "18P4", "18Q4", "18S4", "18T4", "20A4", "20B4", "20D4", "20E4", "21A4", "21B4", "21D4", "21E4", "22A4", "22B4", "24A4", "24B4", 
  "24C4", "24D4", "24E4", "24H4", "24J4", "24K4", "24L4", "24O4", "24P4", "24Q4", "24R4", "24S4", "24T4", "26A4", "26B4", "26C4", 
  "27A4", "27B4", "27C4", "27D4", "28A4", "28C4", "28D4", "28E4", "28F4", "29A4", "30C4", "30D4", "30E4", "30F4", "30G4", "30I4", 
  "32A4", "32C4", "33A4", "36A4", "36D4", "36E4", "36F4", "36G4", "36H4", "36I4", "36J4", "36K4", "36L4", "36Q4", "36S4", "37B4", 
  "38A4", "38B4", "39A4", "39B4", "39C4", "40A4", "40B4", "42B4", "42C4", "42F4", "42G4", "42H4", "42I4", "44A4", "44D4", "46A4", 
  "48A4", "48B4", "48E4", "48F4", "48H4", "48I4", "48J4", "50F4", "52A4", "53A4", "54A4", "54B4", "54C4", "54E4", "55A4", "56B4", 
  "56C4", "58A4", "60B4", "61A4", "62A4", "65A4", "70B4", "72A4", "72B4", "72C4", "72D4", "72E4", "74A4", "75A4", "76A4", "78A4", 
  "81A4", "84A4", "84B4", "108A4", "8A5", "11A5", "12A5", "12B5", "12C5", "12D5", "12E5", "14A5", "14B5", "14E5", "14F5", "15B5", 
  "15C5", "15D5", "16A5", "16B5", "16C5", "16D5", "16E5", "16F5", "16H5", "16I5", "16J5", "16K5", "16M5", "16N5", "16O5", "17A5", 
  "20C5", "20D5", "20F5", "20H5", "20I5", "20J5", "20K5", "21B5", "21C5", "21D5", "21E5", "22B5", "22C5", "24B5", "24C5", "24E5", 
  "24F5", "24G5", "24H5", "24I5", "24J5", "24K5", "24L5", "24M5", "24N5", "24O5", "24P5", "24Q5", "24R5", "24S5", "24T5", "24U5", 
  "24V5", "24W5", "24X5", "24Y5", "24Z5", "24AA5", "24AB5", "26C5", "28B5", "28E5", "28F5", "28G5", "28I5", "30A5", "30E5", "30G5", 
  "30H5", "30I5", "30J5", "30L5", "30M5", "30N5", "30O5", "30R5", "30S5", "32A5", "32B5", "32C5", "32D5", "32E5", "32F5", "32G5", 
  "32H5", "32I5", "32J5", "32K5", "32L5", "32M5", "32N5", "32O5", "33A5", "33B5", "34A5", "34B5", "34C5", "34D5", "35A5", "35C5", 
  "36E5", "36F5", "36I5", "36K5", "36L5", "39A5", "39B5", "40B5", "40C5", "40G5", "40H5", "40I5", "40J5", "40L5", "40M5", "40N5", 
  "40O5", "41A5", "42C5", "42E5", "42F5", "42G5", "42H5", "42I5", "45D5", "45F5", "45G5", "45H5", "48B5", "48C5", "48D5", "48I5", 
  "48J5", "51A5", "51B5", "55B5", "56A5", "56D5", "57C5", "63A5", "63B5", "64A5", "64B5", "64C5", "64D5", "65A5", "68A5", "72B5", 
  "74A5", "75A5", "78A5", "78B5", "78C5", "98A5", "18E6", "21C6", "22A6", "22C6", "24B6", "28I6", "36A6", "36B6", "36D6", "42D6", 
  "48D6", "60A6", "63A6", "63D6", "66A6", "76A6", "79A6", "12A7", "12E7", "14B7", "18B7", "18C7", "18F7", "18J7", "20D7", "20J7", 
  "20N7", "21A7", "22A7", "24C7", "24G7", "24M7", "24S7", "27D7", "27E7", "28D7", "30A7", "30Q7", "30R7", "32A7", "32G7", "32H7", 
  "32K7", "35B7", "36A7", "36H7", "40D7", "40K7", "40T7", "40AC7", "42A7", "42B7", "42C7", "42G7", "42H7", "44A7", "48G7", "48H7", 
  "48S7", "48T7", "48W7", "48AD7", "48AE7", "48AH7", "56A7", "60M7", "60N7", "60P7", "60S7", "62A7", "63B7", "64F7", "66B7", "69A7", 
  "70C7", "72C7", "77A7", "78A7", "80K7", "81A7", "83A7", "85A7", "89A7", "72A8", "88A8", "101A8", "14A9", "16E9", "24E9", "24J9", 
  "24AB9", "30I9", "30K9", "30P9", "32E9", "34A9", "40D9", "40F9", "40J9", "42A9", "45A9", "48D9", "48AG9", "55B9", "57A9", "64A9", 
  "80F9", "80H9", "95A9", "18B10", "18D10", "36K10", "42A10", "42E10", "44C10", "46A10", "54B10", "72D10", "80B10", "84B10", 
  "88C10", "90A10", "92A10", "105A10", "114A10", "115A10", "118A10", "48G11", "94A11", "119A11", "131A11", "104B13", "117A13" ];



/* Cummins and Pauli classify congruence subgroups of SL(2,Z) containing -I up to conjugacy in GL(2,Z).
   Each such class of conjugate groups will consist either of 1 or 2 classes of group up to conjugacy in SL(2,Z).
*/

f:=AssociativeArray();
for i in [1..#cp_data] do
    f[cp_data[i]`name]:=i;
end for;

names:=[k: k in Set(gonality1 cat gonality2 cat gonality3 cat smooth_plane_quintic cat bielliptic)];
gonality1_groups:=[];
gonality2_groups:=[];
gonality3_groups:=[];
quintic_groups:=[];
bielliptic_groups:=[];

m:=AssociativeArray();
for k in names do
    r:=cp_data[f[k]]; // congruence subgroup with given Cummins-Pauli name
    assert r`name eq k;

    if r`level eq 1 then 
        m[k]:=1; 
        gonality1_groups cat:= [ rec<recformat<level, genus, index: RngIntElt, gens : SeqEnum, CPname: MonStgElt> | 
                    level:=r`level, genus:=r`genus, index:=r`index, CPname:=r`name, gens:=[]> ];
        continue k; 
    end if;

    GL2:=GL2Ambient(r`level);
    SL2:=SL2Ambient(r`level);
    B:=GL2![-1,0,0,1];    
    gens1:=r`matgens cat [[-1,0,0,-1]];
    gens1:=[SL2!g: g in gens1];
    gens2:=[SL2!(B*g*B^(-1)): g in gens1];

    H1:=sub<SL2|gens1>;
    H2:=sub<SL2|gens2>;
    if IsConjugate(SL2,H1,H2) then
        m[k]:=1;
    else
        m[k]:=2;
    end if;

    groups:=[ rec<recformat<level, genus, index: RngIntElt, gens : SeqEnum, CPname: MonStgElt> | 
                    level:=r`level, genus:=r`genus, index:=r`index, CPname:=r`name, gens:=[[Integers()!a: a in Eltseq(A)]: A in gens1]> ];
    if m[k] eq 2 then                    
        groups cat:=[ rec<recformat<level, genus, index: RngIntElt, gens : SeqEnum, CPname: MonStgElt> | 
                    level:=r`level, genus:=r`genus, index:=r`index, CPname:=r`name, gens:=[[Integers()!a: a in Eltseq(A)]: A in gens2]> ];
    end if;

    if k in gonality1 then
        gonality1_groups cat:=groups;
    elif k in gonality2 then
        gonality2_groups cat:=groups;
    elif k in gonality3 then
        gonality3_groups cat:=groups;
    end if;

    if k in smooth_plane_quintic then
        quintic_groups cat:=groups;
    end if;
    if k in bielliptic then
        bielliptic_groups cat:=groups;
    end if;

end for;


I:=Open("groups/gonality1.dat", "w");
WriteObject(I, gonality1_groups);

I:=Open("groups/gonality2.dat", "w");
WriteObject(I, gonality2_groups);

I:=Open("groups/gonality3.dat", "w");
WriteObject(I, gonality3_groups);

I:=Open("groups/quintic.dat", "w");
WriteObject(I, quintic_groups);

I:=Open("groups/bielliptic.dat", "w");
WriteObject(I, bielliptic_groups);
