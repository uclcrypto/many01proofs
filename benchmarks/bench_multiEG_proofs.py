from timeit import timeit
from gmpy2 import mpz
import sys
sys.path.append('../')
import multiEG

n_bench = 1000
n_list = [10**i for i in range(1, 7)]
k_range = range(1, 17)
m_range = [1]


def time_this(s, iterations=10):
    t = timeit(s, globals=globals(), number=iterations) * 1000 / iterations
    return t


def print_table(result):
    print("{:>8} {:>8} {:>12} {:>12} {:>12} {:>12}".format(
            "n", "k", "pre.", "enc", "proof", "total"))
    for row in result:
        print("{:>8} {:>8} {:>12} {:>12} {:>12} {:>12}".format(
            row['n'], row['k'], row['pre'], row['t_enc'], row['t_proof'], row['total']))


result_pre = [[0 for k in k_range] for m in m_range]
result_online_enc = [[[0 for n in n_list] for k in k_range] for m in m_range]
result_online_proof = [[[0 for n in n_list] for k in k_range] for m in m_range]

# Exploring all options
for m_idx, m in enumerate(m_range):
    print('\nExploring with m = ', m)
    print('k =', end=' ')
    for k_idx, k in enumerate(k_range):
        print(k, end=' ')
        # Precomputation
        myEG = multiEG.multiEG(m, k)
        myEG.precompute()
        t = time_this('myEG.precompute()')
        result_pre[m_idx][k_idx] = t
        # Computation of encryption
        v = [mpz(i % 2) for i in range(n_bench)]
        t = time_this('myEG.long_multi_enc(v)')
        for n_idx, n in enumerate(n_list):
            result_online_enc[m_idx][k_idx][n_idx] = n * t / n_bench
        # Computation of proof
        c_table = myEG.long_multi_enc(v)
        t = time_this('myEG.long_multi_proof(v, c_table)')
        for n_idx, n in enumerate(n_list):
            result_online_proof[m_idx][k_idx][n_idx] = n * t / n_bench
    print('\nPrecomputation for each k:')
    print(result_pre[m_idx])
    print('Encryption for each k:')
    print(result_online_enc[m_idx])
    print('Proof for each k:')
    print(result_online_proof[m_idx])

# Finding best choice of k for each (m, n) pair
for m_idx, m in enumerate(m_range):
    result = []
    for n_idx, n in enumerate(n_list):
        best_time = 2**100
        for k_idx, k in enumerate(k_range):
            total_time = result_pre[m_idx][k_idx] + result_online_enc[m_idx][k_idx][n_idx] + result_online_proof[m_idx][k_idx][n_idx]
            if total_time < best_time:
                best_time = total_time
                best_k = k
                best_k_idx = k_idx
        result.append({'n': n,
                       'k': best_k,
                       'pre': round(result_pre[m_idx][best_k_idx]),
                       't_enc': round(result_online_enc[m_idx][best_k_idx][n_idx]),
                       't_proof': round(result_online_proof[m_idx][best_k_idx][n_idx]),
                       'total': round(best_time)})
    print("\nm =", m)
    print_table(result)
