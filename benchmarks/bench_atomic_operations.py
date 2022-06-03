from timeit import timeit
from gmpy2 import mpz
import sys
sys.path.append('../')
import multiEG
from group_params import *
import group

def time_this(s, iterations=1000000):
    t = timeit(s, globals=globals(), number=iterations) * 1000 / iterations
    return t

myEG = multiEG.multiEG(m = 1, k = 8)
myEG.precompute()

t_rand = time_this("randq()")
t_mult = time_this("(randq()*randq())%q") - 2*t_rand
t_add = time_this("(randq()+randq())%q") - 2*t_rand
t_exp = time_this("myEG.gpow(randq())") - t_rand

print("t_exp:", t_exp)
print("t_mult:", t_mult)
print("t_add:", t_add)
print("t_rand:", t_rand)

def f(N, t_mult, t_exp, t_add, t_rand):
    n = 0
    while 2**n < N:
        n += 1
    return (4 * n + 4)* t_exp + (N * ((n**2 + 5*n+40)/4)+2*n+4) * t_mult + (N * (6 * n + 6) + 8) * t_add + (7 + 3 * n) * t_rand


n_list = [2**i for i in range(14)]
print("Theoretical time:")
for n in n_list:
    print(n, f(n, t_mult, t_exp, t_add, t_rand))

print("Real time:")
for n in n_list:
    v = [mpz(i % 2) for i in range(n)]
    c_table = myEG.long_multi_enc(v)
    r = [c[2] for c in c_table]
    t = time_this('myEG.long_log_proof(v, r, n)', 10)
    print(n, t)
    
