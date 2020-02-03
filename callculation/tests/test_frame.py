import os
import sys
import unittest

from callculation.hydrodynamic_coefficients import *

cmd_folder = os.path.abspath(os.path.join("..", ".."))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)


class TestHFrm(unittest.TestCase):
    def setUp(self):
        self.h = HFrm(0.5, 1, 1)

    def test_f1(self):
        self.assertEqual(self.h.f1(0), 1)

    def test_rmb(self):
        self.assertEqual(self.h.rmb(self.h.f1, 0.0, 1), 1.0572508754006549)

    def test_fB(self):
        self.assertEqual(self.h.fB, 0.078553374303142381)

    def test_flmb(self):
        self.assertEqual(self.h.flmb, 0.0018334138273949584)

    def test_fmu(self):
        self.assertEqual(self.h.fmu, 0.95363264851429019)


if __name__ == "__main__":
    unittest.main()
