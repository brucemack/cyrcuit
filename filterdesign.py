# Functions related to filter design
import math 

def butterworthNormalizedComponents(n):
    """ Creates the g(k) components of a standard Butterworth low-pass filter. """
    g_list = []
    for k in range(1, n+1):
        g_list.append(2.0 * math.sin((2.0 * k - 1) * math.pi / (2.0 * n)))
    return g_list

def chebyshevNormalizedComponents(n, ripple_db):
    """ Creates the g(k) components of a standard Chebyshev low-pass filter. """
    g_list = []
    d = ripple_db / 8.68589
    B = (1.0 / (2.0 * n)) * math.log((math.exp(d) + 1) / (math.exp(d) - 1))
    B0 = 2 * n * B
    N = (1.0 / 2.0) * (math.exp(B) - math.exp(-B))
    for k in range(1, n+1):
        ak = math.sin(((2.0 * k - 1.0) * math.pi ) / (2.0 * n))
        bk = N ** 2.0 + math.sin(k * math.pi / n) ** 2.0
        if k == 1:
            gk = 2.0 * ak / N
        else:
            gk = (4 * ak_minus_1 * ak) / (bk_minus_1 * gk_minus_1)
        g_list.append(gk)
        # Keep this stuff for the next iteration
        ak_minus_1 = ak
        bk_minus_1 = bk
        gk_minus_1 = gk
    return g_list

def cutoffAdjustmentChebyshev(ripple_db, n):
    """ Creates the w coefficient that adjusts the filter for a -3dB cutoff (as opposed to the 
        normal "maximal ripple" cutoff that would come from the standard parameters).
    """
    e = math.sqrt(10 ** (ripple_db / 10) - 1) 
    y = (1.0 / n) * math.log(1/e + math.sqrt(1.0 / e**2 - 1.0))
    w = (1.0 / 2.0) * (math.exp(y) + math.exp(-y))
    return w

def couplingCoefficientsButterworth(g_params):
    """  Creates the knm coupling coefficients  given a set of g parameters from the
         standard low-pass polynomials.  Expect to get n-1 coupling coefficients 
         back.
    """
    return _couplingCoefficients(1.0, g_params)

def couplingCoefficientsChebyshev(ripple_db, g_params):
    """  Creates the knm coupling coefficients  given a set of g parameters from the
         standard low-pass polynomials.  Expect to get n-1 coupling coefficients 
         back.
    """
    n = len(g_params)
    w = cutoffAdjustmentChebyshev(ripple_db, n)
    return _couplingCoefficients(w, g_params)

def _couplingCoefficients(w, g_params):
    """  Creates the knm coupling coefficients  given a set of g parameters from the
         standard low-pass polynomials.  Expect to get n-1 coupling coefficients 
         back.
    """
    n = len(g_params)
    k_list = []
    for i in range(0, n-1):
        k = (1.0 / w) * (1.0 / ( math.sqrt(g_params[i] * g_params[i+1]) ))
        k_list.append(k)
    return k_list

def endSectionCoefficientChebyshev(ripple_db, g_params):
    """  Creates the q1/qn end coefficient given a set of g parameters from the
         standard low-pass polynomials.  We are assuming symmetry here.
    """
    n = len(g_params)
    # Enforce oddness
    if n % 2 == 0:
        raise Exception("Enven order Chebyshev not supported")
    w = cutoffAdjustmentChebyshev(ripple_db, n)
    return g_params[0] * w

def endSectionCoefficientButterworth(g_params):
    """  Creates the q1/qn end coefficient given a set of g parameters from the
         standard low-pass polynomials.  We are assuming symmetry here.
    """
    return g_params[0]

