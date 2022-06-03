from math import ceil

e_length = 256
group_elem_length = 4096

n_list = [10 ** i for i in [1, 2, 3, 4, 5, 6]]
kg_list = range(2, 20)
ky_list = range(2, 20)
m_list = range(1, 6)


def pre_length(k):
    # Returns the length of the precomputation table in MB
    t = ceil(e_length / k)
    nbr_group_elem = t * (2 ** k - 1)
    bit_length = nbr_group_elem * group_elem_length
    return bit_length / 8000000


def c_proof_length(n, m):
    # Returns the length of the ciphertext + proof in MB
    l = ceil(n / m)
    c_group_elem = l + n  # mod p
    proof_group_elem = 2 * l + 2 * n  # mod q
    bit_length = c_group_elem * group_elem_length + proof_group_elem * e_length
    return bit_length / 8000000


def find_best(n, m):
    best_cost = 2 ** 100
    l = ceil(n / m)
    for kg in kg_list:
        for ky in ky_list:
            pre_g = ceil(e_length / kg) * (2 ** kg - 1)
            pre_y = m * ceil(e_length / ky) * (2 ** ky - 1)
            online_g = (ceil(e_length / kg) - 1) * (2 * l)
            online_y = (ceil(e_length / ky) - 1) * (3 * n)
            total = pre_g + pre_y + online_g + online_y
            if best_cost > total:
                best_cost = total
                best_m = m
                best_kg = kg
                best_ky = ky
                best_pre = pre_g + pre_y
                best_online = online_g + online_y
    # Computing the precomputation memory requirements
    best_pre_size = round(pre_length(best_kg) + m * pre_length(best_ky), 1)
    best_c_proof_size = round(c_proof_length(n, m), 1)
    return {'n': n,
            'cost': best_cost,
            'pre': best_pre,
            'm': best_m,
            'kg': best_kg,
            'ky': best_ky,
            'pre_size': best_pre_size,
            'c_p_size': best_c_proof_size}

def print_table(result):
    print("{:>8} {:>12} {:>12} {:>12} {:>12} {:>5} {:>5} {:>5}".format(
            "n", "cost", "pre.", "pre. (MB)", "c+proof (MB)", "m", "kg", "ky"))
    for row in result:
        print("{:>8} {:>12} {:>12} {:>12} {:>12} {:>5} {:>5} {:>5}".format(
            row['n'], row['cost'], row['pre'], row['pre_size'], row['c_p_size'], row['m'], row['kg'], row['ky']))


# Searching for best choice for each n
result = []
for n in n_list:
    best_cost = 2 ** 100
    for m in m_list:
        res_m = find_best(n, m)
        if res_m['cost'] < best_cost:
            best_cost = res_m['cost']
            best_sol = res_m
    result.append(best_sol)

print("Overall best choice for each n")
print_table(result)

# Exploring what various values of m have to offer
print("\n \n Exploring values for each m")
for m in m_list:
    result = []
    print("\n m = ", m)
    for n in n_list:
        result.append(find_best(n, m))
    print_table(result)

