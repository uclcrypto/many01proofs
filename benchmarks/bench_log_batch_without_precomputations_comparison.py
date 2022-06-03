from timeit import timeit
from gmpy2 import mpz
import sys
sys.path.append('../')
import multiEG
import matplotlib.pyplot as plt
#plt.axes().set_aspect(2)

def time_this(s, iterations=20):
    t = timeit(s, globals=globals(), number=iterations) * 1000 / iterations
    return t


n_log_list = range(0, 10)
n_list = [2**i for i in n_log_list]
#n_list = list(range(1, 17))
k_list = [8, 13, 16]
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

result_log = [[-1 for n in n_list] for k in k_list]
result_batch = [[-1 for n in n_list] for k in k_list]

# Log proof computations
myEG = multiEG.multiEG(m = m, k = 13)
myEG.precompute()

ciphertexts_log = [myEG.long_basic_enc(v) for v in votes]
ciphertexts_batch = [myEG.long_multi_enc(v) for v in votes]
for k_i, k in enumerate(k_list):
    print(k)
    myEG.recompute(k)
    t_pre = time_this('myEG.recompute(k)') 
    result_pre_log[k_i] = t_pre*2/(m+1)
    result_pre_batch[k_i] = t_pre
    t_enc_bench = time_this('myEG.long_multi_enc(v_bench)')
    for n_i, n in enumerate(n_list):
        result_enc_log[n_i][k_i] =  time_this('myEG.long_basic_enc(votes[n_i])')
        if n >= n_bench:
            result_enc_batch[n_i][k_i] = t_enc_bench*n/n_bench
        else:
            result_enc_batch[n_i][k_i] = time_this('myEG.long_multi_enc(votes[n_i])')

        r_log = [c[2] for c in ciphertexts_log[n_i]]
        r_batch = [c[2] for c in ciphertexts_batch[n_i]]
        result_proof_log[n_i][k_i] = time_this('myEG.long_log_proof(votes[n_i], r_log, len(votes[n_i]))')
        result_proof_batch[n_i][k_i] = time_this('myEG.compressed_batch_proof(votes[n_i], r_batch)')
        result_log[k_i][n_i] = result_enc_log[n_i][k_i] + result_proof_log[n_i][k_i]
        result_batch[k_i][n_i] = result_enc_batch[n_i][k_i] + result_proof_batch[n_i][k_i]

    plt.plot(n_log_list, result_log[k_i], label = "Log proof")
    plt.plot(n_log_list, result_batch[k_i], label = "Batch proof")

    plt.xlabel("log_2(n) for n bits")
    #plt.xlabel("number of bits")
    plt.ylabel("Time (ms.)")
    plt.legend()
    plt.savefig("exp_log_nopre_comparison_" + str(k_list[k_i]) + ".png")
    plt.cla()
print(result_log)
print(result_batch)
