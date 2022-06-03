from timeit import timeit
from gmpy2 import mpz
import sys
sys.path.append('../')
import multiEG

n_bench = 16
n_list = [2**(4*i) for i in range(1, 5)]
#k_range = list(range(1, 14))
#m_range = list(range(1, 7))
k_range = list(range(1, 14))
m_range = list(range(1, 7))


def time_this(s, iterations=10):
    t = timeit(s, globals=globals(), number=iterations) * 1000 / iterations
    return t


def print_table(result):
    print("{:>8} {:>8} {:>12} {:>12} {:>12}".format(
            "n", "k", "pre.", "online", "total"))
    for row in result:
        print("{:>8} {:>8} {:>12} {:>12} {:>12}".format(
            row['n'], row['k'], row['pre'], row['online'], row['total']))


result_pre = [[0 for k in k_range] for m in m_range]
result_online = [[[0 for n in n_list] for k in k_range] for m in m_range]

# Exploring all options
for k_idx, k in enumerate(k_range):
    myEG = multiEG.multiEG(1, k)
    myEG.precompute()
    t_pre = time_this('myEG.precompute()')
    print('\nExploring for k = ', k)
    print('m =')
    for m_idx, m in enumerate(m_range):
        print(m)
        # Precomputation
        myEG = multiEG.multiEG(m, k)
        myEG.precompute()
        result_pre[m_idx][k_idx] =  (t_pre * (m + 1))/2
        # Computation
        v = [mpz(i % 2) for i in range(n_bench)]
        t = time_this('myEG.long_multi_enc(v)')
        for n_idx, n in enumerate(n_list):
            result_online[m_idx][k_idx][n_idx] = n*t/n_bench
    print('\nPrecomputation for each k:')
    print([result_pre[m_idx][k_idx] for m_idx in range(len(m_range))])
    print('Computation for each m and k:')
    print([result_online[m_idx][k_idx] for m_idx in range(len(m_range))])

# print(result_pre)
# print(result_online)

# Finding best choice of k for each (m, n) pair

import matplotlib.pyplot as plt
plt.axes().set_aspect(2)
for n_idx, n in enumerate(n_list):
    result = []
    res = []
    for m_idx, m in enumerate(m_range):
        best_time = 2**100
        for k_idx, k in enumerate(k_range):
            total_time = result_pre[m_idx][k_idx] + result_online[m_idx][k_idx][n_idx]
            if total_time < best_time:
                best_time = total_time
                best_k = k
                best_k_idx = k_idx
        result.append({'n': n,
                       'k': best_k,
                       'pre': result_pre[m_idx][best_k_idx],
                       'online': result_online[m_idx][best_k_idx][n_idx],
                       'total': best_time})
        res.append((result[-1]["pre"] + result[-1]["online"])/n)
    plt.plot(m_range, res, "-o", label = "{} votes".format(n))
plt.legend(loc="best")
plt.ylim(ymin=0)
plt.xlabel("Number of bases")
plt.ylabel("Time (ms.)/bit")
plt.savefig("exp_m.png")
plt.show()
"""
for m_idx, m in enumerate(m_range):
    result = []
    for n_idx, n in enumerate(n_list):
        best_time = 2**100
        for k_idx, k in enumerate(k_range):
            total_time = result_pre[m_idx][k_idx] + result_online[m_idx][k_idx][n_idx]
            if total_time < best_time:
                best_time = total_time
                best_k = k
                best_k_idx = k_idx
        result.append({'n': n,
                       'k': best_k,
                       'pre': result_pre[m_idx][best_k_idx],
                       'online': result_online[m_idx][best_k_idx][n_idx],
                       'total': best_time})
    print("m =", m)
    print_table(result)
""" 
