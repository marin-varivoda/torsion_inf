This folder contains data from Cummins and Pauli (https://mathstats.uncg.edu/sites/pauli/congruence/) who classify low genus congruence subgroups.
The data file "CPdata.dat" is obtained from their code by running the following in Magma:
    load "pre.m";
    load "csg.m";
    load "csg24.dat";
    filename:="CPdata.dat";
    I:=Open(filename,"w");
    WriteObject(I, L);

Starting with "CPdata.dat", we then make a more useful version of the data by running "produce_data.m".
This outputs text which we have hardcoded into the file "CPdata.m";

At the end of "CPdata.m" are functions that we will use to access the Cummins-Pauli data.