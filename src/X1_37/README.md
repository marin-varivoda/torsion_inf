Minimal ST-relation code for X1(37), F_2, degree 9.

Files:
  X1_37_ST_minimal.m  main mathematical driver
  ff_search.m         generic effective-divisor search routines
  verify.m            one-time map/order verification
  worker.m            one chunk of the search
  run_first.sh        verify once, then run chunks in parallel

Run:
  ./run_first.sh 16 16

The code uses the quotient image relation

  q(S,T) = S^38 + S^32*T^5 + S^24*T^13 + S^14*T^20
           + S^8*T^29 + S^8*T^25 + S^2*T^34 + T^37 + T^33.

It pulls back two rational places on the normalization of q(S,T)=0:

  infinity:       1/S, T^36/S^37
  S=0, T/S=1:    S, T/S+1

Their pullbacks both have degree 36, and D0 is the difference
  D0 = infinity fiber - (S=0, T/S=1) fiber.
