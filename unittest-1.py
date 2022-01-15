import unittest
from filterdesign import *

class TestFilters(unittest.TestCase):

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

    def test_lpf_1(self):
        """ Design a low pass filter and denormalize. """
        g_list = butterworthNormalizedComponents(5)
        fc = 50000000
        r0 = 50
        c1 = denormalizeC(g_list[0], fc, r0)
        l2 = denormalizeL(g_list[1], fc, r0)
        c3 = denormalizeC(g_list[2], fc, r0)
        l4 = denormalizeL(g_list[3], fc, r0)
        c5 = denormalizeC(g_list[4], fc, r0)
        self.assertAlmostEqual(393.4e-12, c1, places=1)
        self.assertAlmostEqual(2.575e-6, l2, places=3)
        self.assertAlmostEqual(1273e-12, c3, places=0)

    def test_lpf_2(self):
        """ Design a low pass filter and denormalize. 
            This example is from IRFD pg 84.
        """
        g_list = chebyshevNormalizedComponentsAdjusted3db(2, 0.1)
        fc = 100000000
        bw = 400000
        r0 = 50

        # Check the gk values from the example
        self.assertAlmostEqual(1.638, g_list[0], places=3)
        self.assertAlmostEqual(1.2087, g_list[1], places=4)

        # This funtion assumes a load of 1, so we need to re-normalize 
        # in order to adjust to a source of 1.
        rs, rl = chebyshevNormalizedSourceAndLoadResistances(2, 0.1)
        self.assertAlmostEqual(1.355, 1.0 / rs, places=3)

        

    def assertListAlmostEqual(self, l0, l1, places):
        self.assertEqual(len(l0), len(l1))
        for k in range(0, len(l0)):
            self.assertAlmostEqual(l0[k], l1[k], places=places)

if __name__ == '__main__':
    unittest.main()
