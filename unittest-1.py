import unittest
from filterdesign import *

class TestFilters(unittest.TestCase):

    def assertListAlmostEqual(self, l0, l1, places):
        self.assertEqual(len(l0), len(l1))
        for k in range(0, len(l0)):
            self.assertAlmostEqual(l0[k], l1[k], places=places)

    def test_0(self):
        # Creating the Butterworth coefficients
        # Taken from ARRL Handbook 2016 pg. 11.12
        expected_g_list = [0.7654, 1.848, 1.848, 0.7654 ]
        g_list = butterworthNormalizedComponents(4)
        self.assertListAlmostEqual(expected_g_list, g_list, 3)

        # Creating Chebyshev components
        # Taken from IRFD pg. 65
        expected_g_list = [ 0.8430, 0.6220 ]
        g_list = chebyshevNormalizedComponents(2, 0.1)
        self.assertListAlmostEqual(expected_g_list, g_list, 3)

        expected_g_list = [ 1.0316, 1.1474, 1.0316 ]
        g_list = chebyshevNormalizedComponents(3, 0.1)
        self.assertListAlmostEqual(expected_g_list, g_list, 3)

    def test_1(self):
        # Taken from Hayward "Designing Narrow-Badwith Ladder Filters" pg. 40
        self.assertAlmostEqual(1.1347, cutoffAdjustmentChebyshev(0.1, 5), places=4)

        # Taken from Hayward "Designing Narrow-Badwith Ladder Filters" pg. 41
        # NOTE: We are only using the first two g parameters here
        g_list = chebyshevNormalizedComponents(5, 0.1)
        self.assertAlmostEqual(0.7028, couplingCoefficientsChebyshev(0.1, g_list)[0], places=4)
 
        # Taken from Hayward "Designing Narrow-Badwith Ladder Filters" pg. 41
        # NOTE: We are only using the first g parameter here
        g_list = chebyshevNormalizedComponents(5, 0.1)
        self.assertAlmostEqual(1.3013, endSectionCoefficientChebyshev(0.1, g_list), places=4)

        # Taken from Hayward "Designing Narrow-Badwith Ladder Filters" pg. 5
        g_list = butterworthNormalizedComponents(4)
        expected_k_list = [ 0.841, 0.541, 0.841 ]
        k_list = couplingCoefficientsButterworth(g_list)
        self.assertListAlmostEqual(expected_k_list, k_list, 3)
        self.assertAlmostEqual(0.7654, endSectionCoefficientButterworth(g_list), places=4)


if __name__ == '__main__':
    unittest.main()
