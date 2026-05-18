# Classification

This folder contains the code for the paper "Classification of modular curves with low gonality".
It classifies the congruences subgroups `Gamma` of `SL(2,Z)`, up to conjugacy, for which the modular curve `X_Gamma` has gonality 1, gonality 2, gonality 3, or is bielliptic.

Loading the file `all_checks.m` in Magma from this directory will verify the classification!
(   The verification is broken up into pieces given by the files  `gonality2.m`, `gonality3.m`, `plane_qunitic.m`, `bielliptic.m`. 
    The file `classification_functions.m` contains functions needed for these computations. )
The verifications took a little under 2 hours on a machine at Cornell's math department; your times may vary.

The file `save_groups.m` saves the groups in the classification to files in the subfolder `groups`.
