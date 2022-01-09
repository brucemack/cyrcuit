# Functions related to filter design
import math 

def cutoffAdjustmentChebyshev(rippleDb, n):
    """ Gets the w coefficient that adjusts the filter for a -3dB cutoff (as opposed to the 
        normal "maximal ripple" cutoff that would come from the standard parameters).
    """
    e = math.sqrt(10 ** (rippleDb / 10) - 1) 
    y = (1.0 / n) * math.log(1/e + math.sqrt(1.0 / e**2 - 1.0))
    w = (1.0 / 2.0) * (math.exp(y) + math.exp(-y))
    return w

def couplingCoefficientsChebyshev(rippleDb, g_params):
    n = len(g_params)
    w = cutoffAdjustmentChebyshev(rippleDb, n)
    k_list = []
    for i in range(0, n-1):
        k = (1.0 / w) * (1.0 / ( math.sqrt(g_params[i] * g_params[i+1]) ))
        k_list.append(k)
    return k_list
