import sys
sys.path.append('../')

import unittest
import multiEG
from gmpy2 import mpz
# from group_params import *


class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.k = 8
        self.m = 4
        self.n = 23
        self.myEG = multiEG.multiEG(m = self.m,k = self.k)
        self.myEG.precompute()

        self.l = -(-self.n//self.m)
        self.votes = [mpz(i % 2) for i in range(self.n)]

    def test_baseline_proof(self):
        v = mpz(0)
        c, d, r = self.myEG.baseline_enc(v)
        proof = self.myEG.baseline_proof(v, r, full_proof=True)
        self.assertTrue(self.myEG.verify_baseline_proof(c, d, *proof))
        
        v = mpz(1)
        c, d, r = self.myEG.baseline_enc(v)
        proof = self.myEG.baseline_proof(v, r, full_proof=True)
        self.assertTrue(self.myEG.verify_baseline_proof(c, d, *proof))
       
        v = mpz(2)
        c, d, r = self.myEG.baseline_enc(v)
        proof = self.myEG.baseline_proof(v, r, full_proof=True)
        self.assertFalse(self.myEG.verify_baseline_proof(c, d, *proof))
    
    def test_long_multi_proof(self):

        ciphertexts = self.myEG.long_baseline_enc(self.votes)
        proof = self.myEG.long_baseline_proof(self.votes, ciphertexts, full_proof=True)
        for i in range(self.l):
            self.assertEqual(len(ciphertexts[i]), 3)
            self.assertEqual(len(proof[i]), 8)
            c, d = ciphertexts[i][0:2]
            partial_proof = proof[i]
            self.assertTrue(self.myEG.verify_baseline_proof(c, d, *partial_proof))

        # Checking if short proofs have the right format
        proof = self.myEG.long_baseline_proof(self.votes, ciphertexts, full_proof=False)
        for i in range(self.l):
            self.assertEqual(len(proof[i]), 4)

    def test_multi_proof(self):
        c, d, r = self.myEG.multi_enc(self.votes[:4])
        proof = self.myEG.multi_proof(self.votes[:4], r, full_proof=True)
        self.assertTrue(self.myEG.verify_multi_proof(c, d, *proof))
        
        v = [mpz(i) for i in range(self.m)]
        c, d, r = self.myEG.multi_enc(v)
        proof = self.myEG.multi_proof(v, r, full_proof=True)

    def test_long_multi_proof(self):

        ciphertext = self.myEG.long_multi_enc(self.votes)
        proof = self.myEG.long_multi_proof(self.votes, ciphertext, full_proof=True)

        for i in range(self.l):
            self.assertEqual(len(ciphertext[i]), 3)
            self.assertEqual(len(proof[i]), 8)
            c, d = ciphertext[i][0:2]
            partial_proof = proof[i]
            slice_length = len(d)
            self.assertTrue(slice_length <= self.m)
            self.assertEqual(len(d), slice_length)
            self.assertEqual(len(partial_proof[1]), slice_length)
            self.assertEqual(len(partial_proof[2]), slice_length)
            self.assertEqual(len(partial_proof[6]), slice_length)
            self.assertEqual(len(partial_proof[7]), slice_length)
            self.assertTrue(self.myEG.verify_multi_proof(c, d, *partial_proof))

        # Checking if short proofs have the right format
        proof = self.myEG.long_multi_proof(self.votes, ciphertext, full_proof=False)
        for i in range(self.l):
            self.assertEqual(len(ciphertext[i]), 3)
            self.assertEqual(len(proof[i]), 4)
            c, d = ciphertext[i][0:2]
            partial_proof = proof[i]
            slice_length = len(d)
            self.assertTrue(slice_length <= self.m)
            self.assertEqual(len(d), slice_length)
            self.assertEqual(len(partial_proof[2]), slice_length)
            self.assertEqual(len(partial_proof[3]), slice_length)
    
    def test_CDS_proof(self):
        v = mpz(0)
        c, d, r = self.myEG.basic_enc(v)
        proof = self.myEG.CDS_proof(v, r, full_proof=True)
        self.assertTrue(self.myEG.verify_CDS_proof(c, d, *proof))

        v = mpz(1)
        c, d, r = self.myEG.basic_enc(v)
        proof = self.myEG.CDS_proof(v, r, full_proof=True)
        self.assertTrue(self.myEG.verify_CDS_proof(c, d, *proof))
        
        v = mpz(2)
        c, d, r = self.myEG.basic_enc(v)
        proof = self.myEG.CDS_proof(v, r, full_proof=True)
        self.assertFalse(self.myEG.verify_CDS_proof(c, d, *proof))

    def test_long_CDS_proof(self):
        ciphertext = self.myEG.long_basic_enc(self.votes)

        # Checking if full proofs have the right format and verify
        proof = self.myEG.long_CDS_proof(self.votes, ciphertext, full_proof=True)
        for i in range(self.n):
            self.assertEqual(len(ciphertext[i]), 3)
            self.assertEqual(len(proof[i]), 8)
            c, d = ciphertext[i][0:2]
            partial_proof = proof[i]
            self.assertTrue(self.myEG.verify_CDS_proof(c, d, *partial_proof))

        # Checking if short proofs have the right format
        proof = self.myEG.long_multi_proof(self.votes, ciphertext, full_proof=False)
        for i in range(self.n):
            self.assertEqual(len(ciphertext[i]), 3)
            self.assertEqual(len(proof[i]), 4)

    
    def test_batch_proof(self):
        ciphertext = self.myEG.long_multi_enc(self.votes)
        
        # Checking if full proofs have the right format and verify
        proof = self.myEG.batch_proof(self.votes, ciphertext, full_proof = True)
        self.assertEqual(len(proof), 9)

        self.assertEqual(len(proof[1]), self.l)
        self.assertEqual(len(proof[2]), self.n)
        self.assertEqual(len(proof[3]), self.n)
        self.assertEqual(len(proof[6]), self.l)
        self.assertEqual(len(proof[7]), self.n)
        self.assertEqual(len(proof[8]), self.n)
        
        c = [c[0] for c in ciphertext]
        d = [d[1] for d in ciphertext]
        self.assertTrue(self.myEG.verify_batch_proof(c, d, *proof))
        
        # Checking if full proofs have the right format and verify
        proof = self.myEG.batch_proof(self.votes, ciphertext, full_proof = False)
        self.assertEqual(len(proof), 4)

        self.assertEqual(len(proof[1]), self.l)
        self.assertEqual(len(proof[2]), self.n)
        self.assertEqual(len(proof[3]), self.n)
    
    def test_compressed_batch_proof(self):
        ciphertext = self.myEG.long_multi_enc(self.votes)

        proof = self.myEG.compressed_batch_proof(self.votes, ciphertext, full_proof = True)
        self.assertEqual(len(proof), 9)

        self.assertEqual(len(proof[1]), self.l)
        self.assertEqual(len(proof[2]), self.n)
        self.assertEqual(len(proof[4]), self.m)
        self.assertEqual(len(proof[6]), self.l)
        self.assertEqual(len(proof[7]), self.n)

        c = [c[0] for c in ciphertext]
        d = [c[1] for c in ciphertext]
        self.assertTrue(self.myEG.verify_compressed_batch_proof(c, d, *proof))
        
        # Checking if full proofs have the right format and verify
        proof = self.myEG.compressed_batch_proof(self.votes, ciphertext, full_proof = False)
        self.assertEqual(len(proof), 4)

        self.assertEqual(len(proof[1]), self.l)
        self.assertEqual(len(proof[2]), self.n)
        self.assertEqual(len(proof[3]), self.n)
    
    def test_log_proof(self):
        for vote_size in [1, 2, 5, 10, 16]:
            with self.subTest(i = vote_size):
                votes = [mpz(i % 2) for i in range(vote_size)]
                ciphertext = self.myEG.long_basic_enc(votes)
                r = [c[2] for c in ciphertext]
                com = [[c[0], c[1]] for c in ciphertext]

                proof = self.myEG.log_proof(votes, r, full_proof=True)
                self.assertTrue(self.myEG.verify_log_proof(com, *proof))
    
    def test_long_log_proof(self):
        for block_size in [1, 2, 4, 16]:
            with self.subTest(i = block_size):
                votes = [i % 2 for i in range(16)]
                ciphertext = self.myEG.long_basic_enc(votes)
                proof = self.myEG.long_log_proof(votes, ciphertext, block_size, full_proof=True)
                self.assertTrue(self.myEG.long_verify_log_proof(ciphertext, proof, block_size))


if __name__ == '__main__':
    unittest.main()
