"""
3. Pour la preuve en log, j'ai l'impression que l'histoire la plus intéressante
est l'évolution du speedup de la preuve log par rapport à une preuve de base
(la baseline? ou ce sera peut-être plus marquant en donnant la batch
        compressé?) en fonction de n (on prend n= 2^i pour i allant de 1 à 16
                par exemple, en prenant le k optimal pour chaque valeur de n).
        Ici, on compare juste le coût des preuves, en laissant tomber le coût
        du chiffrement. Message: (i) même pour de petits n genre n=8, la preuve
        en log est sympa; (ii) quand n grandit, la preuve en log écrase les
        autres (speedup énorme) -- c'est montrer cela qu'on ne met pas le coût
        du chiffrement dans le graphe. (iii) même si la preuve en log est
        asymptotiquement plus coûteuse, cela ne se manifeste pas pour des
        valeurs réalistes. (iv) Par contre, on observe un pic bien intéressant
        de l'efficactité autour de quelques milliers de preuves à calculer.
        Leçon: si on beaucoup de preuves à calculer, mieux vaut les calculer
        par paquets de quelques milliers: ce sera plus rapide, et la taille de
        la preuve reste très petite par rapport à celle des chiffrés de toute
        manière.
"""
import matplotlib.pyplot as plt
from math import sqrt
from timeit import timeit
from gmpy2 import mpz
import sys
sys.path.append('../')
import multiEG

def time_this(s, iterations=10):
    t = timeit(s, globals=globals(), number=iterations) * 1000 / iterations
    return t

result_compressed_batch_proof = []
result_log_proof = []

n_bench = 64

n_log_list = range(17)
n_list = [2**x for x in n_log_list]
k_list = list(range(1, 14))


result_baseline_proof = [[0 for k in k_list] for n in n_list]
result_log_proof = [[0 for k in k_list] for n in n_list]

speedup = []
"""
for n_i, n in enumerate(n_list):
    k_v = -1
    k_best = -1
    #print(n_i)
    for k_i, k in enumerate(k_list):
        myEG = multiEG.multiEG(1, k)
        if k != 0: 
            myEG.precompute()
        t_pre = time_this('myEG.precompute()', 10)

        v = [mpz(i % 2) for i in range(n_bench)]
        t_enc = (time_this('myEG.long_multi_enc(v)', 10)* n / n_bench)
        c_table = myEG.long_multi_enc(v)
        r = [c[2] for c in c_table]
        t_proof = (time_this('myEG.long_baseline_proof(v, r)', 10) * n / n_bench)

        #print(k, t_pre, t_enc, t_proof, t_pre + t_enc+ t_proof)
        if k_best == -1 or t_enc + t_pre + t_proof < k_v:
            k_v = t_enc + t_pre + t_proof
            k_best = k
        else:
            break
    print(k_best)
"""

k_opt = [1, 3, 3, 4, 5, 5, 6, 7, 8, 8, 9, 10, 11, 11, 11, 12, 12]
#k_opt = [13]*17

# Exploring all options
myEG = multiEG.multiEG(1, 1)
myEG.precompute()
for n_i, n in enumerate(n_list):
    myEG.recompute(k_opt[n_i])
    print(n)
    # Computation
    v = [mpz(i % 2) for i in range(n)]
    c_table = myEG.long_multi_enc(v)
    r = [c[2] for c in c_table]
    t = time_this('myEG.long_log_proof(v, r, len(v))', 50)
    result_log_proof[n_i] = t
    
    v = [mpz(i % 2) for i in range(n_bench)]
    c_table = myEG.long_multi_enc(v)
    r = [c[2] for c in c_table]
    t = time_this('myEG.long_baseline_proof(v, r)', 25)
    result_baseline_proof[n_i] = (t*n/n_bench)

    speedup.append(result_baseline_proof[n_i]/result_log_proof[n_i])
    print(n, k_opt[n_i], speedup[-1])
#for n_i, n in enumerate(n_list):
#    k_best = -1
#    for k_i, k in enumerate(k_list):
#        if k_best == -1 or result_baseline_proof[k_best] > result_baseline_proof[k_i]:
#            k_best = k_i
#    print(k_best)
#    speedup.append(round(result_baseline_proof[n_i][k_best])/round(result_log_proof[n_i][k_best]))

plt.axes().set_aspect(0.1)
plt.plot(n_log_list, speedup, "-")
plt.ylim(ymin=0)
plt.xlabel("log_2(n) for n bits")
plt.ylabel("Speed-up factor")
plt.savefig("exp_log.png")
plt.show()

