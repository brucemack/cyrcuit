import unittest
from filterdesign import *

class TestFilters(unittest.TestCase):

    def test_1(self):
        # Take from Hayward "Designing Narrow-Badwith Ladder Filters" pg. 40
        self.assertAlmostEqual(1.1347, cutoffAdjustmentChebyshev(0.1, 5), places=4)

        # Take from Hayward "Designing Narrow-Badwith Ladder Filters" pg. 41
        # NOTE: We are only using the first two g parameters here
        g_list = [ 1.1468, 1.3712, 1, 1, 1 ]
        self.assertAlmostEqual(0.7028, couplingCoefficientsChebyshev(0.1, g_list)[0], places=4)
 
if __name__ == '__main__':
    unittest.main()
