import sys
sys.path.append('../')

import unittest
import multiEG
from group_params import *


class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.m = 4
        self.v = [mpz(i) for i in range(self.m)]
        self.myEG = multiEG.multiEG(self.m)
        self.myEG.precompute()

    def test_key(self):
        self.assertEqual(len(self.myEG.y), self.m)
        self.assertEqual(self.myEG.y[1], pow(g, self.myEG.x[1], p))
        self.assertNotEqual(self.myEG.y[1], pow(g, self.myEG.x[0], p))

    def test_enc(self):
        c, d, r = self.myEG.multi_enc(self.v)
        self.assertEqual(c, pow(g, r, p))
        self.assertEqual(d[1], pow(self.myEG.y[1], r + self.v[1], p))
        self.assertEqual(d[2] * pow(c, q - self.myEG.x[2], p) % p, pow(self.myEG.y[2], self.v[2], p))

    def test_partial_enc(self):
        c, d, r = self.myEG.multi_enc(self.v)
        self.assertEqual(c, pow(g, r, p))
        self.assertEqual(d[1], pow(self.myEG.y[1], r + self.v[1], p))
        self.assertEqual(d[0] * pow(c, q - self.myEG.x[0], p) % p, pow(self.myEG.y[0], self.v[0], p))

    def test_long_multi_enc(self):
        n = 23
        v = [mpz(i) for i in range(n)]
        t = self.myEG.long_multi_enc(v)
        self.assertEqual(len(t), -(-len(v)//self.m))
        self.assertEqual(len(t[-1][1]), len(v) % self.m)
        self.assertEqual(t[0][0], pow(g, t[0][2], p))
        self.assertEqual(t[3][1][2], pow(self.myEG.y[2], t[3][2] + v[self.m * 3 + 2], p))
    
    def test_baseline_enc(self):
        c, d, r = self.myEG.baseline_enc(self.v[0])
        self.assertEqual(c, pow(g, r, p))
        self.assertEqual(d, pow(self.myEG.y[0], r, p) * pow(g, self.v[0], p) % p)

    def test_long_baseline_enc(self):
        t = self.myEG.long_baseline_enc(self.v)

        self.assertEqual(len(t), len(self.v))
        self.assertEqual(t[0][0], pow(g, t[0][2], p))
        self.assertEqual(t[-1][1], pow(self.myEG.y[0], t[-1][2], p) * pow(g, self.v[-1], p) % p)

    def test_basic_enc(self):
        c, d, r = self.myEG.basic_enc(self.v[0])
        self.assertEqual(c, pow(g, r, p))
        self.assertEqual(d, pow(self.myEG.y[0], r + self.v[0], p))

    def test_long_basic_enc(self):
        t = self.myEG.long_basic_enc(self.v)

        self.assertEqual(len(t), len(self.v))
        self.assertEqual(t[0][0], pow(g, t[0][2], p))
        self.assertEqual(t[-1][1], pow(self.myEG.y[0], t[-1][2] + self.v[-1], p))

if __name__ == '__main__':
    unittest.main()
