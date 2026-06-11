X1(37), F_2 degree-9 translate test: compact first-run package
================================================================

Files
-----

X1_37_prop49_F2_short_complete.m
    Compact Magma script.  It keeps the completed-run safeguards but removes
    alternate slow paths and unused chunk modes.  The computation is always the
    successful type-chunked, streamed, lazy fast-linear W0 computation.

X1_37_X0plus_map_data.m
    Hard-coded S37_hard_coeffs and T37_hard_coeffs.

run_first.sh
    First-run driver that was used to verify our results.  It verifies the hard-coded map and the cached generator
    pair <1,4> once, then launches all chunks.  Workers skip only checks already
    proved in the verification log.

logs_first_runs contains the log of a run that verified the result, and was run on 16 cores.

Usage
-----
    ./run_first.sh 16 16

The first argument is the number of concurrent Magma processes.  The second is
the number of degree-type chunks; it defaults to 16.

Safety markers
--------------

The driver refuses to launch workers unless verify.log contains:

    MAP_VERIFICATION: SUCCESS
    cyclic_delta_ok           = true
    Using candidate fiber pair <1,4> as a verified order-5 generator over F_2.

The final successful run prints:

    GLOBAL_RESULT: SUCCESS
