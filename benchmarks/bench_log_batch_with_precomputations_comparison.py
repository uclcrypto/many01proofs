from timeit import timeit
from gmpy2 import mpz
import sys
sys.path.append('../')
import multiEG
import matplotlib.pyplot as plt
#plt.axes().set_aspect(2)

def time_this(s, iterations=10):
    t = timeit(s, globals=globals(), number=iterations) * 1000 / iterations
    return t


n_log_list = range(0, 10)
n_list = [2**i for i in n_log_list]
k_list = list(range(1, 13))
votes = [[mpz(i % 2) for i in range(n)] for n in n_list]
m = 4

n_bench = 64
v_bench = [mpz(i % 2) for i in range(n_bench)]

result_pre_log = [-1 for k in k_list]
result_enc_log = [[-1 for k in k_list] for n in n_list]
result_proof_log = [[-1 for k in k_list] for n in n_list]
result_pre_batch = [-1 for k in k_list]
result_enc_batch = [[-1 for k in k_list] for n in n_list]
result_proof_batch = [[-1 for k in k_list] for n in n_list]

result_log = [-1] * len(n_list)
result_batch = [-1] * len(n_list)

# Log proof computations
myEG = multiEG.multiEG(m = m, k = 13)
myEG.precompute()

ciphertexts = [myEG.long_basic_enc(v) for v in votes]
for k_i, k in enumerate(k_list):
    print(k_i)
    myEG.recompute(k)
    result_pre_log[k_i] = time_this('myEG.recompute(k)')*2/(m+1)
    t_enc_bench = time_this('myEG.long_basic_enc(v_bench)')
    for n_i, n in enumerate(n_list):
        if n >= n_bench:
            result_enc_log[n_i][k_i] = t_enc_bench*n/n_bench
        else:
            result_enc_log[n_i][k_i] =  time_this('myEG.long_basic_enc(votes[n_i])')
        r = [c[2] for c in ciphertexts[n_i]]
        result_proof_log[n_i][k_i] = time_this('myEG.long_log_proof(votes[n_i], r, len(votes[n_i]))')

for n_i, n in enumerate(n_list):
    result_log[n_i] = sum(min(zip(result_pre_log, result_enc_log[n_i],
        result_proof_log[n_i]), key = lambda x: sum(x)))
    
# Batch proof computations

ciphertexts = [myEG.long_multi_enc(v) for v in votes]
for k_i, k in enumerate(k_list):
    print(k_i)
    myEG.recompute(k)
    result_pre_batch[k_i] = time_this('myEG.recompute(k)')
    t_enc_bench = time_this('myEG.long_multi_enc(v_bench)')
    for n_i, n in enumerate(n_list):
        if n >= n_bench:
            result_enc_batch[n_i][k_i] = t_enc_bench*n/n_bench
        else:
            result_enc_batch[n_i][k_i] = time_this('myEG.long_multi_enc(votes[n_i])')

        r = [c[2] for c in ciphertexts[n_i]]
        result_proof_batch[n_i][k_i] = time_this('myEG.compressed_batch_proof(votes[n_i], r)')

for n_i, n in enumerate(n_list):
    result_batch[n_i] = sum(min(zip(result_pre_batch, result_enc_batch[n_i],
        result_proof_batch[n_i]), key = lambda x: sum(x)))
    
print(result_log)
print(result_batch)

plt.plot(n_log_list, result_log, label = "Log proof")
plt.plot(n_log_list, result_batch, label = "Batch proof")

plt.xlabel("log_2(n) for n bits")
plt.ylabel("Time (ms.)")
plt.legend()
plt.savefig("exp_log_comparison.png")
#plt.show()
