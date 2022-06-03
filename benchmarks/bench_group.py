from timeit import timeit
import sys
sys.path.append('../')
from gmpy2 import powmod
import group
from group_params import *

import secrets

n_bench = 1000
k_range = range(1, 18)

def time_this(s, iterations=10):
    # returns running time in milliseconds
    t = timeit(s, globals=globals(), number=iterations) * 1000 / iterations
    return t

e_table = [randq() for _ in range(n_bench)]
t_gmpy2 = time_this('[powmod(g, e, p) for e in e_table]')
print("Exponentiation using gmpy2: ", t_gmpy2)

print("{:>8} {:>12} {:>12}".format("k", "pre", str(n_bench) + " exp"))

for k in k_range:
    # Timing the precomputation
    t_pre = time_this('group.PowRadix(g, k)')
    # Timing the exponentiation
    gpow = group.PowRadix(g, k, False).pow
    e_table = [randq() for _ in range(n_bench)]
    t_exp = time_this('[gpow(e) for e in e_table]')
    print("{:>8} {:>12} {:>12}".format(k, round(t_pre, 1), round(t_exp, 1)))
