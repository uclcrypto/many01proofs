from gmpy2 import powmod, invert
from group_params import *
import group

class multiEG:
    def __init__(self, m = 1, k = 8):
        self.m = m
        self.k = k
        self.x = [randq() for i in range(m)]
        self.y = [powmod(g, self.x[i], p) for i in range(m)]
        self.gpow = lambda x: powmod(g, x, p)
        self.ypow = [lambda x: powmod(self.y[i], x, p) for i in range(m)]

    def precompute(self):
        g_radix = group.PowRadix(g, self.k)
        self.gpow = g_radix.pow
        y_radix = [group.PowRadix(self.y[i], self.k) for i in range(self.m)]
        self.ypow = [y_radix[i].pow for i in range(self.m)]

    def recompute(self, k):
        # re-precompute the fixed-based exponentiations with another precomputation factor k
        self.k = k
        g_radix = group.PowRadix(g, self.k)
        self.gpow = g_radix.pow
        y_radix = [group.PowRadix(self.y[i], self.k) for i in range(self.m)]
        self.ypow = [y_radix[i].pow for i in range(self.m)]
    
    def baseline_enc(self, v):
        # textbook ElGamal encryption: single base, precomputation
        r = randq()
        c = self.gpow(r)
        d = self.ypow[0](r) * self.gpow(v) % p
        return c, d, r

    def long_baseline_enc(self, votes):
        # textbook ElGamal encryption: single base, produces len(v) ciphertexts
        ciphertexts = []
        for v in votes:
            c, d, r = self.baseline_enc(v)
            ciphertexts.append([c, d, r])
        return ciphertexts 

    def basic_enc(self, v):
        # textbook ElGamal encryption: single base, precomputation
        r = randq()
        c = self.gpow(r)
        d = self.ypow[0](r + v) % p
        return c, d, r

    def long_basic_enc(self, votes):
        # textbook ElGamal encryption: single base, produces len(v) ciphertexts
        ciphertexts = []
        for v in votes:
            c, d, r = self.basic_enc(v)
            ciphertexts.append([c, d, r])
        return ciphertexts 
    
    def multi_enc(self, votes):
        # returns ciphertexts and the associated randomness, for a vector v of length at most m
        assert len(votes) <= self.m
        r = randq()
        c = self.gpow(r)
        d = [self.ypow[i](r + votes[i]) for i in range(len(votes))]
        return c, d, r

    def long_multi_enc(self, votes):
        # returns an array of ciphertexts and associated randomness, for a vector v of any length
        ciphertexts = []
        for i in range(0, len(votes), self.m):
            c, d, r = self.multi_enc(votes[i:i+self.m])
            ciphertexts.append([c, d, r])
        return ciphertexts 

    def baseline_proof(self, v, r, full_proof = False):
        # Computes a proof that (c, d) is an encryption of v using randomness r

        s = randq()
        t = randq()
        w = randq()

        # Computing commitments
        temp = self.gpow(w)
        a0 = self.gpow(s)
        a1 = self.gpow(t)
        if v == 0:
            b0 = self.ypow[0](s)
            b1 = self.ypow[0](t) * temp % p
        else:
            b0 = self.ypow[0](s) * temp % p
            b1 = self.ypow[0](t)

        # Picking challenge, benchmarking only
        # CAUTION: In reality, this will be computed via Fiat-Shamir!!!
        e = randq()

        # Computing response
        e0 = ((1-v)*e - w)  % q
        e1 = (v*e + w) % q

        f0 = (s + e0 * r) % q
        f1 = (t + e1 * r) % q
        
        if not full_proof:
            return e0, e1, f0, f1
        else:
            return a0, b0, a1, b1, e0, e1, f0, f1
    
    def verify_baseline_proof(self, c, d, a0, b0, a1, b1, e0, e1, f0, f1):
        result = True
        if self.gpow(f0) != powmod(c, e0, p) * a0 % p:
            result = False
        if self.gpow(f1) != powmod(c, e1, p) * a1 % p:
            result = False
        if self.ypow[0](f0) != powmod(d, e0, p) * b0 % p:
            result = False
        if self.gpow(e1) * self.ypow[0](f1) % p != powmod(d, e1, p) * b1 % p:
            result = False
        return result
    
    def long_baseline_proof(self, votes, ciphertexts, full_proof = False):
        proof = []
        for v, c in zip(votes, ciphertexts):
            rolling_proof = self.baseline_proof(v, c[2], full_proof)
            proof.append(rolling_proof)
        return proof 

    def long_verify_baseline_proof(self, ciphertexts, proofs):
        for proof, ciphertext in zip(ciphertexts, proofs):
            if not self.verify_baseline_proof(*ciphertext, *proof):
                return False
        return True

    def multi_proof(self, v, r, full_proof=False):
        # Computes a proof that (c, d) is an encryption of a 0-1 vector v using randomness r
        # Computing commitment
        s = randq()
        ar = self.gpow(s)
        w = [randq() for _ in range(len(v))]
        a = [self.ypow[i]((w[i] + s) % q) for i in range(len(v))]
        t = [randq() for _ in range(len(v))]
        b = [self.ypow[i]((w[i]*v[i] + t[i]) % q) for i in range(len(v))]
        u = [self.gpow(t[i]) for i in range(len(v))]

        # Picking challenge, benchmarking only
        # CAUTION: In reality, this will be computed via Fiat-Shamir!!!
        e = randq()

        # Computing response
        fr = (s + e * r) % q
        fa = [(w[i] + e * v[i]) % q for i in range(len(v))]
        fb = [(t[i] + r * (e - fa[i])) % q for i in range(len(v))]

        if not full_proof:
            return e, fr, fa, fb
        else:
            return ar, a, b, u, e, fr, fa, fb

    def verify_multi_proof(self, c, d, ar, a, b, u, e, fr, fa, fb):
        result = True
        if self.gpow(fr) != ar * pow(c, e, p) % p:
            result = False
        for i in range(len(d)):
            if self.ypow[i]((fa[i] + fr) % q) != a[i] * pow(d[i], e, p) % p:
                result = False
            if self.ypow[i](fb[i]) != b[i] * pow(d[i], (e - fa[i]) % q, p) % p:
                result = False
            if self.gpow(fb[i]) != u[i] * pow(c, (e - fa[i]) % q, p) % p:
                result = False
        return result

    def long_multi_proof(self, v, ciphertext, full_proof=False):
        proof = []
        for i, c in enumerate(ciphertext):
            rolling_v = v[i*self.m: (i + 1) * self.m]
            if full_proof:
                ar, a, b, u, e, fr, fa, fb = self.multi_proof(rolling_v, c[2], full_proof=True)
                proof.append([ar, a, b, u, e, fr, fa, fb])
            else:
                e, fr, fa, fb = self.multi_proof(rolling_v, c[2], full_proof=False)
                proof.append([e, fr, fa, fb])
        return proof

    def CDS_proof(self, v, r, full_proof=False):
        # Computes a proof that (c, d) is an encryption of a 0-1 vector v using randomness r
        # Uses precompution for the exponentiations
        # Caution: this only works for basic ElGamal ciphertexts (that is, one public key element)
        # Computing simulated proof
        sim_e = randq()
        sim_f = randq()
        sim_a = self.gpow((sim_f - r * sim_e) % q)
        sim_b = self.ypow[0]((sim_f - (r - 1 + 2 * v) * sim_e) % q)
        # Computing real proof
        s = randq()
        true_a = self.gpow(s)
        true_b = self.ypow[0](s)
        # Picking challenge, benchmarking only
        # CAUTION: In reality, this will be computed via Fiat-Shamir and sim_e!!!
        true_e = randq()

        # Computing response
        true_f = (s + true_e * r) % q

        if not full_proof:
            if v == 0:
                return true_e, sim_e, true_f, sim_f
            else:
                return sim_e, true_e, sim_f, true_f
        else:
            if v == 0:
                return true_a, true_b, sim_a, sim_b, true_e, sim_e, true_f, sim_f
            else:
                return sim_a, sim_b, true_a, true_b, sim_e, true_e, sim_f, true_f

    def verify_CDS_proof(self, c, d, a0, b0, a1, b1, e0, e1, f0, f1):
        result = True
        # Verifying 0 proof
        if self.gpow(f0) != a0 * pow(c, e0, p) % p:
            result = False
        if self.ypow[0](f0) != b0 * pow(d, e0, p) % p:
            result = False
        # Verifying 1 proof
        if self.gpow(f1) != a1 * pow(c, e1, p) % p:
            result = False
        d_over_y = d * invert(self.y[0], p) % p
        if self.ypow[0](f1) != b1 * pow(d_over_y, e1, p) % p:
            result = False
        return result

    def long_CDS_proof(self, votes, ciphertext, full_proof=False):
        proof = []
        for v, c in zip(votes, ciphertext):
            partial_proof = self.CDS_proof(v, c[2], full_proof)
            proof.append(partial_proof)
        return proof
  
    def batch_proof(self, votes, ciphertext, full_proof = False):
        # Picking challenge, benchmarking only
        # CAUTION: In reality, this will be computed via Fiat-Shamir!!!
        gamma = randq()
   
        n = len(votes)
        l = len(ciphertext)

        s = [randq() for i in range(l)]
        w = [randq() for i in range(n)]
        t = [randq() for i in range(n)]
       
        # Computing commitments
        A = [self.gpow(s[i]) for i in range(l)] 
        B = [self.ypow[i % self.m]((w[i] + s[i // self.m]) % q) for i in range(n)]
        E = [self.ypow[i % self.m]((t[i] + votes[i]*w[i]) % q) for i in range(n)]
        
        T = mpz(0)
        gamma_curr = 1
        for i in range(n):
            T = (T + t[i] * gamma_curr) % q 
            gamma_curr = (gamma_curr * gamma) % q
        E0 = self.gpow(T)

        # Picking challenge, benchmarking only
        # CAUTION: In reality, this will be computed via Fiat-Shamir!!!
        e = randq()

        # Computing responses
        zr = [(e * ciphertext[j][2] + s[j]) % q for j in range(l)] 
        zv = [(e*votes[i] + w[i]) % q for i in range(n)]
        zt = [(ciphertext[i // self.m][2]*(e - zv[i]) + t[i]) % q for i in range(n)]
        
        if not full_proof:
            return e, zr, zv, zt
        else:
            return gamma, A, B, E, E0, e, zr, zv, zt
    

    def verify_batch_proof(self, c, d, gamma, A, B, E, E0, e, zr, zv, zt):
        result = True

        n = len(B)
        l = len(c)
        for i in range(l):
            if self.gpow(zr[i]) != (pow(c[i], e, p)*A[i]) % p:
                result = False

        for i in range(n):
            j = i // self.m
            k = i % self.m
            if self.ypow[i % self.m]((zv[i] + zr[i // self.m]) % q) != (pow(d[i // self.m][i % self.m], e, p) * B[i]) % p:
                result = False
            if self.ypow[i % self.m](zt[i]) != (pow(d[i // self.m][i % self.m], e - zv[i], p) * E[i]) % p:
                result = False

        A_gamma = mpz(1)
        z_gamma = mpz(0)
        gamma_curr = mpz(1)
        for i in range(n):
            A_gamma = A_gamma * pow(c[i // self.m], (e - zv[i]) * gamma_curr % q, p) % p
            z_gamma = (z_gamma + zt[i] * gamma_curr)% q
            gamma_curr = (gamma_curr * gamma) % q
        
        if (A_gamma * E0) % p != self.gpow(z_gamma):
            result = False
        return result
    
    def compressed_batch_proof(self, votes, ciphertext, full_proof = False):
        # Picking challenge, benchmarking only
        # CAUTION: In reality, this will be computed via Fiat-Shamir!!!
        gamma = randq()
   
        m = self.m
        n = len(votes)
        l = len(ciphertext)

        s = [randq() for i in range(l)]
        w = [randq() for i in range(n)]
        t = [randq() for i in range(n)]
        
        # Computing commitments
        A = [self.gpow(s[i]) for i in range(l)] 
        B = [self.ypow[i % self.m]((w[i] + s[i // self.m]) % q) for i in range(n)]
        2^8 + 2^6 + 2^2 + 1
        T = mpz(0)
        gamma_curr = 1
        for i in range(n):
            T = (T + t[i] * gamma_curr) % q 
            gamma_curr = (gamma_curr * gamma) % q
        A0 = self.gpow(T)

        B0 = [mpz(0) for i in range(m)]
        for k in range(m):
            T = mpz(0)
            gamma_curr = 1
            for i in range(l):
                if i*m + k < n:
                    T = (T + (t[i*m + k] + votes[i*m + k]*w[i*m + k]) * gamma_curr) % q 
                    gamma_curr = (gamma_curr * gamma) % q

            B0[k] = self.ypow[k](T)

        # Picking challenge, benchmarking only
        # CAUTION: In reality, this will be computed via Fiat-Shamir!!!
        e = randq()

        # Computing responses
        zr = [(e * ciphertext[j][2] + s[j]) % q for j in range(l)] 
        zv = [(e*votes[i] + w[i]) % q for i in range(n)]
        zt = [(ciphertext[i // self.m][2]*(e - zv[i]) + t[i]) % q for i in range(n)]
    
        if not full_proof:
            return e, zr, zv, zt
        else:
            return gamma, A, B, A0, B0, e, zr, zv, zt
    
    def verify_compressed_batch_proof(self, c, d, gamma, A, B, A0, B0, e, zr, zv, zt):
        result = True
        n = len(B)
        l = len(c)
        for i in range(l):
            if self.gpow(zr[i]) != (pow(c[i], e, p)*A[i]) % p:
                result = False

        for i in range(n):
            if self.ypow[i % self.m]((zv[i] + zr[i // self.m]) % q) != (pow(d[i // self.m][i % self.m], e, p) * B[i]) % p:
                result = False

        A_gamma = mpz(1)
        z_gamma = mpz(0)
        gamma_curr = mpz(1)
        for i in range(n):
            A_gamma = A_gamma * pow(c[i // self.m], (e - zv[i]) * gamma_curr % q, p) % p
            z_gamma = (z_gamma + zt[i] * gamma_curr)% q
            gamma_curr = (gamma_curr * gamma) % q
        
        if (A_gamma * A0) % p != self.gpow(z_gamma):
            result = False

        for k in range(self.m):
            B_gamma = mpz(1)
            z_gamma = mpz(0)
            gamma_curr = mpz(1)
            for i in range(l):
                if i*self.m + k < n:
                    B_gamma = B_gamma * pow(d[i][k], (e - zv[i*self.m + k]) * gamma_curr % q, p) % p
                    z_gamma = (z_gamma + zt[i*self.m + k] * gamma_curr)% q
                    gamma_curr = (gamma_curr * gamma) % q
            
            if (B_gamma * B0[k]) % p != self.ypow[k](z_gamma):
                result = False
        
        return result


    def log_proof(self, v, r, full_proof=False):
        # Computes a proof that (c, d) is an encryption of a 0-1 vector v using randomness r
        N = len(v)

        # Initial round
        # Get the smallest n such that 2**n >= len(v)
        n = 0
        while 2**n < len(v):
            n += 1
        
        # Picking challenge, benchmarking only
        # CAUTION: In reality, this will be computed via Fiat-Shamir!!!
        gamma = randq()

        y = [mpz(0)]*len(v)
        gamma_i = [mpz(1)]*len(v)
        for i in range(N):
            if i > 0:
                gamma_i[i] = gamma_i[i-1] * gamma % q
            y[i] = gamma_i[i]*(1 - v[i])

        cu  = [mpz(0)]*n
        cw  = [mpz(0)]*n
        rho = [mpz(0)]*(n+1)
        mu = [randq() for _ in range(n)]
        nu = [randq() for _ in range(n)]
        U = [mpz(0)]*n
        W = [mpz(0)]*n

        alpha = [mpz(0)] * n 
        alpha_inv = [mpz(0)] * n
        alpha_I = [mpz(1)]*N
        alpha_I_inv = [mpz(1)]*N

        # Iterative rounds
        for k in range(n):
            for i_p in range(0, N, 2**(k+1)):
                Y0, Y1 = mpz(0), mpz(0)
                X0, X1 = mpz(0), mpz(0)

                for i_m in range(2**k):
                    if i_p + i_m >= N:
                        break

                    if v[i_p + i_m] == 0:
                        Y0 = (Y0 + y[i_p + i_m] * alpha_I_inv[i_m]) % q
                    else:
                        X0 = (X0 + alpha_I[i_m]) % q
                
                for i_m in range(2**k, 2**(k+1)):
                    if i_p + i_m >= N:
                        break

                    if v[i_p + i_m] == 0:
                        Y1 = (Y1 + y[i_p + i_m] * alpha_I_inv[i_m-2**k]) % q
                    else:
                        X1 = (X1 + alpha_I[i_m - 2**k]) % q


                U[k] = (U[k] + X0 * Y1) % q
                W[k] = (W[k] + X1 * Y0) % q
                
            cu[k] = [self.gpow(mu[k]), self.ypow[0]((U[k] + mu[k])%q)]
            cw[k] = [self.gpow(nu[k]), self.ypow[0]((W[k] + nu[k])%q)]

            # Picking challenge, benchmarking only
            # CAUTION: In reality, this will be computed via Fiat-Shamir!!!
            alpha[k] = randq()
            alpha_inv[k] = invert(alpha[k], q)
            
            alpha_I[2**(k)] = alpha[k]
            alpha_I_inv[2**k] = alpha_inv[k]
            for i in range(2**k):
                if 2**k + i >= N:
                    break
                alpha_I[2**k + i] = alpha_I[i]*alpha_I[2**k] % q
                alpha_I_inv[2**k + i] = alpha_I_inv[i]*alpha_I_inv[2**k] % q

            rho[k+1] = (mu[k]*alpha_inv[k] + rho[k] + nu[k]*alpha[k]) % q
      
        # Penultimate round
        X = mpz(0)
        R = mpz(0)
        S = mpz(0)
        Y = mpz(0)

        gamma_curr = mpz(1)
        for i in range(N):
            X = (X + v[i]*alpha_I[i]) % q 
            Y = (Y + y[i]*alpha_I_inv[i]) % q
            R = (R + r[i] * alpha_I[i]) % q
            S = (S - (r[i] * gamma_i[i] % q)* alpha_I_inv[i]) % q

        T = (rho[n] - S*X) % q

        Xp = randq()
        Yp = randq()
        Rp = randq()
        Sp = randq()
        Tp = randq()

        Cp = [self.gpow(Rp), self.ypow[0]((Xp + Rp) % q)]
        Ep = [self.gpow((S * Xp + Tp) % q), self.ypow[0]((Y*Xp + S*Xp + Tp) % q)]

        # Picking challenge, benchmarking only
        # CAUTION: In reality, this will be computed via Fiat-Shamir!!!
        beta = randq()

        # Final round
        zx = (beta*X + Xp) % q
        zr = (beta*R + Rp) % q
        zt = (beta*T + Tp) % q

        return cu, cw, alpha, beta, gamma, Cp,  Ep, zx, zr, zt
   
    def verify_log_proof(self, c, cu, cw, alpha, beta, gamma, Cp, Ep, zx, zr, zt):
        result = True
     
        # Get the smallest n such that 2**n >= len(c)
        n = 0
        while 2**n < len(c):
            n += 1
        assert n == len(alpha)

        cv = [[mpz(1), mpz(1)]]
        alpha_inv = [mpz(0)] * n

        for i in range(n):
            alpha_inv[i] = invert(alpha[i], q)

            new_cv = [mpz(0), mpz(0)]
            for j in range(2):
                new_cv[j] = (powmod(cu[i][j], alpha_inv[i], p) * cv[i][j] * powmod(cw[i][j], alpha[i], p)) % p
            cv.append(new_cv)
       
        C, D = [mpz(1), mpz(1)], [mpz(1), mpz(1)]
        gamma_curr = mpz(1)
        for i in range(len(c)):
            alpha_curr = mpz(1)
            alpha_inv_curr = mpz(1)
            for b in range(n):
                if i & (1 << b):
                    alpha_curr = (alpha_curr * alpha[b]) % q
                    alpha_inv_curr = (alpha_inv_curr * alpha_inv[b]) % q
                
            C[0] = (C[0] * powmod(c[i][0], alpha_curr, p)) % p
            C[1] = (C[1] * powmod(c[i][1], alpha_curr, p)) % p
            d = [powmod(invert(c[i][0], p), gamma_curr, p), powmod((self.y[0] * invert(c[i][1], p)) % p, gamma_curr, p)]
            D[0] = (D[0] * powmod(d[0], alpha_inv_curr, p)) % p
            D[1] = (D[1] * powmod(d[1], alpha_inv_curr, p)) % p

            gamma_curr = gamma_curr * gamma % q

        E = cv[n]
        if ((powmod(C[0], beta, p) * Cp[0]) % p != self.gpow(zr) or
            (powmod(C[1], beta, p) * Cp[1]) % p != self.ypow[0]((zx + zr) % q) or
            (powmod(E[0], beta, p) * Ep[0]) % p != (powmod(D[0], zx, p)*self.gpow(zt)) % p or
            (powmod(E[1], beta, p) * Ep[1]) % p != (powmod(D[1], zx, p)*self.ypow[0](zt)) % p):
            result = False 

        return result

    def long_log_proof(self, v, ciphertext, block_size, full_proof=False):
        table = []
        n = len(v)
        for i in range((n+block_size-1)//block_size):
            rolling_v = v[i*block_size: (i + 1) * block_size]
            rolling_r = [c[2] for c in ciphertext[i*block_size: min((i + 1) * block_size, n)]]
            proof = self.log_proof(rolling_v, rolling_r, full_proof)
            table.append(proof)
        return table
    
    def long_verify_log_proof(self, ciphertext, proof, block_size, full_proof=False):
        table = []
        for i in range(len(proof)):
            rolling_c = [(c[0], c[1]) for c in ciphertext[i*block_size: min((i + 1) * block_size, len(ciphertext))]]
            if not self.verify_log_proof(rolling_c, *proof[i]):
                return False
        return True
