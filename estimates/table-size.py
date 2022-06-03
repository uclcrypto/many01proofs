from math import ceil

modulus_length = 4096
e_size = 256
max_k = 18

for k in range(3, max_k + 1):
    t = ceil(e_size/k)
    cells = t * (2**k -1)
    table_size_bits = modulus_length * cells
    print("k = ", k, "storage (MB) = ", round(table_size_bits/8000000, 2), "mult = ", t-1)
