import matplotlib.pyplot as plt
from math import sqrt
from timeit import timeit, repeat
from gmpy2 import mpz
import sys
sys.path.append('../')
import multiEG

#def time_this(s, setup = "",iterations=10):
#    t = min(repeat(s, setup = setup, globals=globals(), number=iterations)) * 1000 / iterations
#    return t
def time_this(s, iterations=10):
    t = timeit(s, globals=globals(), number=iterations) * 1000 / iterations
    return t

m = 4
n_bench = m * 4

n_list = [2**(4*i) for i in range(1, 5)]
k_list = list(range(1, 15))

result_pre = [0 for k in k_list] 
result_enc = [[0 for k in k_list] for n in n_list]
result_multi_proof = [[0 for k in k_list] for n in n_list]
result_batch_proof = [[0 for k in k_list] for n in n_list]
result_compressed_batch_proof = [[0 for k in k_list] for n in n_list]


best_k = [2, 5, 7, 11]
# Computed before with:
"""
for k_i, k in enumerate(k_list):
    print(k_i)
    myEG = multiEG.multiEG(m = m, k = k)
    t = time_this('myEG.precompute()', 1)
    myEG.precompute()
    result_pre[k_i] = t
    for n_i, n in enumerate(n_list):
        v = [mpz(i % 2) for i in range(n_bench)]
        t = time_this('myEG.long_multi_enc(v)', 1)
        result_enc[n_i][k_i] = (t*n/n_bench)

        # Computation of proof
        c_table = myEG.long_multi_enc(v)
        r = [c[2] for c in c_table]
        t = time_this('myEG.long_multi_proof(v, r)', 1)
        result_multiEG[n_i][k_i] = (t*n/n_bench)

        #%t = time_this('myEG.batch_proof(v, r)')
        #result_batch_proof[n_i][k_i] = t
        
        #t = time_this('myEG.compressed_batch_proof(v, r)')
        #result_compressed_batch_proof[n_i][k_i] = t

result = [result_pre, result_enc, result_multiEG, result_batch_proof, result_compressed_batch_proof]

best_k = [-1 for i in n_list]
for n_i, n in enumerate(n_list):
    for k_i, k in enumerate(k_list):
        if best_k[n_i] == -1 or result_pre[k_i] + result_enc[n_i][k_i] + result_batch_proof[n_i][k_i] < result_pre[best_k[n_i]] + result_enc[n_i][best_k[n_i]] + result_batch_proof[n_i][best_k[n_i]]:
            best_k[n_i] = k_i
print(best_k)
"""
for n_i, n in enumerate(n_list):
    print(n_i)
    k = best_k[n_i]
    myEG = multiEG.multiEG(m = m, k = k)
    t = time_this('myEG.precompute()', 10)
    myEG.precompute()
    result_pre[k] = round(t)
    v = [mpz(i % 2) for i in range(n)]
    t = time_this('myEG.long_multi_enc(v)', 10)
    result_enc[n_i][k] = round(t)

    # Computation of proof
    c_table = myEG.long_multi_enc(v)
    r = [c[2] for c in c_table]
    t = time_this('myEG.long_multi_proof(v, r)', 10)
    result_multi_proof[n_i][k] = round(t)

    t = time_this('myEG.batch_proof(v, r)')
    result_batch_proof[n_i][k] = round(t)
    
    t = time_this('myEG.compressed_batch_proof(v, r)')
    result_compressed_batch_proof[n_i][k] = round(t)

table = [[n, best_k[n_i], result_pre[best_k[n_i]], result_multiEG[n_i][best_k[n_i]], result_batch_proof[n_i][best_k[n_i]], result_compressed_batch_proof[n_i][best_k[n_i]]] for n_i, n in enumerate(n_list)]
print("\\begin{center}")
print("\\begin{tabular}{|c|c|c|c|c|c|c}")
labels = ["Number of votes", "Optimal k", "Precomputation", "Encryption", "Multi EG", "Batch Proof", "Compressed Batch Proof"]
print("\hline\t" + "&".join(labels))
for n_i, n in enumerate(n_list):
    print("\hline\t "  + "& ".join(map(lambda x : str(x), table[n_i])) + "\\\\")
print("\\end{tabular}")
print("\\end{center}")




