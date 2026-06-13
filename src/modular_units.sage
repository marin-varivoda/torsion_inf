from mdsage import *
from mdsage.maartens_sage_functions import *
from sage.all import euler_phi


T7 =  [(1, n) for n in list(range(1, 25)) + [26, 27, 28, 30]] \
    + [(2, 2*n) for n in range(1, 11)]

T8 =  [(1, n) for n in list(range(1, 29)) + [30, 32, 36]] \
    + [(2, 2*n) for n in range(1, 13)] \
    + [(3, 3*n) for n in range(1, 5)] \
    + [(4, 4*n) for n in range(1, 4)] \
    + [(5, 5), (6, 6)]

T9 =  [(1, n) for n in list(range(1, 29)) + [30, 36]] \
    + [(2, 2*n) for n in range(1, 13)]

has_func = None
modular_unit = None

for d, T in [(7, T7), (8, T8), (9, T9)]:
    print(f"\n===== d={d} =====")

    verified = 0
    failures = []

    for M, N in T:
        G = Gamma11(M, N // M)
        genus = G.genus()
        target_degree = d // euler_phi(M)

        if (genus == 0):
            has_func = True
            modular_unit = None

            print(f"(M,N)=({M:>1},{N:>2})  "
                  f"Gamma11({M},{N//M:>2})  "
                  f"g={genus:>2}  "
                  f"deg={target_degree}  "
                  f"has_func={has_func}  ")

        else:
            # We only construct explicit modular unit when genus > 0
            has_func, modular_unit = has_modular_unit_of_degree(G, target_degree, l2_step=2)
            print(f"(M,N)=({M:>1},{N:>2})  "
                  f"Gamma11({M},{N//M:>2})  "
                  f"g={genus:>2}  "
                  f"deg={target_degree}  "
                  f"has_func={has_func}  "
                  f"modular_unit={modular_unit}")

        if has_func:
            verified += 1
        else:
            failures.append((M, N))

    print(f"\nSummary for d={d}: verified {verified}/{len(T)} pairs.")

    if failures:
        print(f"Failures for d={d}: {failures}")
    else:
        print(f"All pairs verified for d={d}.")