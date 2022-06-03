import sys
sys.path.append('../')
import unittest
from gmpy2 import xmpz
import group
from group_params import *

class TestPowRadix(unittest.TestCase):
    def setUp(self):
        self.k = 7
        self.g_radix = group.PowRadix(g, self.k)

    def test_ElectionGuardParams(self):
        self.assertEqual(q.bit_length(), 256)
        self.assertEqual(p.bit_length(), 4096)
        self.assertEqual(q, 2**256 - 189)

    def test_precomputation(self):
        table = self.g_radix.table
        self.assertEqual(len(table), -(-q.bit_length() // self.k))
        self.assertEqual(len(table), self.g_radix.table_length)

    def test_pow(self):
        e = xmpz(189723973987397928729872)
        self.assertEqual(self.g_radix.pow(e), pow(g, e, p))


if __name__ == '__main__':
    unittest.main()
