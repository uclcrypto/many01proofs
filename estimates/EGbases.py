from math import ceil

n = 100000
range_m = 20 # How many ElGamal bases?
range_k = 15 # How much precomputation?

def EGcost(n):
    # Computing the total cost of encrypting n values (precomputation included)
    # We test various values of k (precomputation) and m (EG width)
    # We print the resulting table

    #first index is k, second is m
    total = [[0 for i in range(range_m)] for j in range(range_k)]

    for m in range(1, range_m+1):
        l = ceil(n/m)
        for k in range(1, range_k+1):
            t = ceil(256/k)
            pre = (m+1) * t * (2**k-1)
            comp = (n+l) * (t-1)
            total[k-1][m-1] = pre + comp

    print("Cost for", n, "exponentiations, precomputation included.")
    print("Each line corresponds to a value of k")
    i = 0
    for row in total:
        i += 1
        print(i, row)

EGcost(n)
