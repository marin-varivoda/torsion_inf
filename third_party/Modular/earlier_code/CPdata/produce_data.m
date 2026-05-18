
// Run this code in Magma from the directory this file is in.
cp_data:=[];
filename:="CPdata.dat";  
I:=Open(filename, "r"); 
_,cp_data:=ReadObjectCheck(I); 

// We only keep track of the entries that are useful for us.
CPRec := recformat< level, genus, index : RngIntElt,
                supergroups, matgens, cusps,
                name: MonStgElt>;
cp_data_new:=[];
for i in [1..#cp_data] do
    r:=cp_data[i];
    s:= rec<CPRec | level:=r`level, genus:=r`genus, index:=r`index, supergroups:=r`supergroups, cusps:=r`cusps, name:=r`name>;
    s`matgens:=[[Integers()!a:a in A] : A in r`matgens];
    cp_data_new:=cp_data_new cat [s];
end for;

// The following will output the data as text that Magma can run.  This is hardcoded into the file "CPdata.m".
PrintFileMagma("temp.m",cp_data_new);