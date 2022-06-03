from group_params import *
from montgomery import Montgomery
# Class supporting exponentiation with precomputation (radix method)
class PowRadix:
    def __init__(self, base, k=1, montgomery = False):

        self.montgomery = montgomery

        if self.montgomery:
            self.M = Montgomery(p, 4096)
            base = self.M.repr(base)
            self.r_one = self.M.repr(mpz(1))

        # Precomputation in base
        self.table_length = -(-q.bit_length() // k)  # Double negative to take the ceiling
        self.k = k
        table = []
        row_base = base
        running_base = row_base
        for _ in range(self.table_length):
            if self.montgomery:
                row = [self.r_one]
            else:
                row = [1]
            for j in range(1, 2**k):
                row.append(running_base)
                if self.montgomery:
                    running_base = self.M.mult(running_base, row_base)
                else:
                    running_base = running_base * row_base % p
            table.append(row)
            row_base = running_base
        self.table = table
        self.max_exponent = 2**256

    def pow(self, e):
        # Computing pow(base, e, p)
        assert 0 <= e < self.max_exponent, 'exponent out of range'
        y = mpz(1)
        if self.montgomery:
            y = self.r_one

        for i in range(self.table_length):
            e_slice = e[i * self.k: (i+1) * self.k]
            if self.montgomery:
                y = self.M.mult(y, self.table[i][e_slice]) 
            else:
                y = y * self.table[i][e_slice] % p
        if self.montgomery:
            y = self.M.reduce(y)
        return y
