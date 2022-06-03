from timeit import timeit
import sys
sys.path.append('../')
import multiEG

n_bench = 1000
n_list = [10**i for i in range(1, 6)]
k_range = range(1, 15)

def time_this(s, iterations=10):
    t = timeit(s, globals=globals(), number=iterations) * 1000 / iterations
    return t


def print_table(result):
    print("{:>8} {:>8} {:>12} {:>12} {:>12} {:>12}".format(
            "n", "k", "pre", "enc", "proof.", "total"))
    for row in result:
        print("{:>8} {:>8} {:>12} {:>12} {:>12} {:>12}".format(
            row['n'], row['k'], row['t_pre'], row['t_enc'], row['t_proof'], row['t_pre'] + row['t_enc'] + row['t_proof']))


myEG = multiEG.multiEG(m=1, k=8)

# Encryption
v = [i % 2 for i in range(n_bench)]
t = time_this('myEG.long_basic_enc(v, precomputation=False)')

result = []
for n in n_list:
    result.append({'n': n, 'k': 'no', 't_pre': 0, 't_enc': round(n*t/n_bench)})

# CDS proof
myEG.precompute()
c = myEG.long_basic_enc(v)
t = time_this('myEG.long_CDS_proof(v, c, precomputation=False)')

for r in result:
    r['t_proof'] = round(r['n']*t/n_bench)

print("Without fixed_base exponentiation")
print_table(result)

# Testing with precomputation
result_k = {}

for k in k_range:
    myEG = multiEG.multiEG(m=1, k=k)
    myEG.precompute()
    t_pre = time_this('myEG.precompute()')

    # Encryption
    v = [i % 2 for i in range(n_bench)]
    t_enc = time_this('myEG.long_basic_enc(v, precomputation=True)')
    # CDS proof
    c = myEG.long_basic_enc(v, precomputation=True)
    t_proof = time_this('myEG.long_CDS_proof(v, c, precomputation=True)')

    result_k[k] = {}
    for n in n_list:
        t_enc_n = n*t_enc/n_bench
        t_proof_n = n*t_proof/n_bench
        result_k[k][n] = {'t_pre': t_pre, 't_enc': t_enc_n, 't_proof': t_proof_n}

# Looking for best k:
result = []
for n in n_list:
    best_time = 2**100
    for k in k_range:
        total_time = result_k[k][n]['t_pre'] + result_k[k][n]['t_enc'] + result_k[k][n]['t_proof']
        if total_time < best_time:
            best_time = total_time
            best_k = k
    result.append({'n': n, 'k': best_k, 't_pre': round(result_k[best_k][n]['t_pre']), 't_enc': round(result_k[best_k][n]['t_enc']), 't_proof': round(result_k[best_k][n]['t_proof'])})

print("With fixed_base exponentiation")
print_table(result)
print(result_k)



