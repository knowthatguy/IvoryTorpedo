import os
import sys
import unittest
cmd_folder = os.path.abspath(os.path.join('..', '..'))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

from callculation.hydrodynamic_coefficients import *

class TestHFrm(unittest.TestCase):

    def test_f1(self):
        val = HFrm(1,1,1).f1(0)
        self.assertEqual(val, 1)

if __name__ == '__main__':
    unittest.main()