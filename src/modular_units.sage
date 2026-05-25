from mdsage import *
from mdsage.maartens_sage_functions import *

X1_N_with_deg7_function = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30]
X1_N_with_deg8_function = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 32, 36]
X1_N_with_deg9_function = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 36]

for N in X1_N_with_deg7_function:
    G = Gamma11(1, N)
    g = G.genus() # maybe compare this with list?
    print(has_modular_unit_of_degree(G,7, l2_step=2))

for N in X1_N_with_deg8_function:
    G = Gamma11(1, N)
    g = G.genus() # maybe compare this with list?
    print(has_modular_unit_of_degree(G, 8, l2_step=2))

for N in X1_N_with_deg9_function:
    G = Gamma11(1, N)
    g = G.genus() # maybe compare this with list?
    print(has_modular_unit_of_degree(G,9, l2_step=2))





