# Functions related to filter design
import math 

def cutoffAdjustmentChebyshev(ripple_db, n):
    """ Creates the w coefficient that adjusts the filter for a -3dB cutoff (as opposed to the 
        normal "maximal ripple" cutoff that would come from the standard parameters).
    """
    e = math.sqrt(10 ** (ripple_db / 10) - 1) 
    y = (1.0 / n) * math.log(1/e + math.sqrt(1.0 / e**2 - 1.0))
    w = (1.0 / 2.0) * (math.exp(y) + math.exp(-y))
    return w

def couplingCoefficientsChebyshev(ripple_db, g_params):
    """  Creates the knm coupling coefficients  given a set of g parameters from the
         standard low-pass polynomials.  Expect to get n-1 coupling coefficients 
         back.
    """
    n = len(g_params)
    w = cutoffAdjustmentChebyshev(ripple_db, n)
    k_list = []
    for i in range(0, n-1):
        k = (1.0 / w) * (1.0 / ( math.sqrt(g_params[i] * g_params[i+1]) ))
        k_list.append(k)
    return k_list

def endSectionCoefficientChebyshev(ripple_db, g_params, rs_ohms):
    """  Creates the q1/qn end coefficient given a set of g parameters from the
         standard low-pass polynomials.  We are assuming symmetry here.
    """
    n = len(g_params)
    # Enforce oddness
    if n % 2 == 0:
        raise Exception("Enven order Chebyshev not supported")
    w = cutoffAdjustmentChebyshev(ripple_db, n)
    return g_params[0] * w / rs_ohms

