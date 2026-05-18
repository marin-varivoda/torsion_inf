# Groups

This folder contains the groups in the classifications from the paper "Classification of modular curves with low gonality".
It classifies the congruences subgroups `Gamma` of `SL(2,Z)`, up to conjugacy, for which the modular curve `X_Gamma` has gonality 1, 
gonality 2, gonality 3, or is bielliptic (also isomorphic to a smooth plane quintic).

To load the gonality 3 groups to a sequence `groups` use the code:

    I:=Open("gonality3.dat", "r");
    _,groups:=ReadObjectCheck(I);
You can load the other groups in the same manner.

The file `count.m` verifies the numerical values in the main theorem and Table 1 of the paper.

The file `CPnames.m` contains the classification in terms of the Cummins-Pauli labels.