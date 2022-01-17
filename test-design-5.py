# Demonstration of the design of N=4 bandpass filter with 
# series crystals.

# Source: Hayward "Designing Narrow-Bandwidth Ladder Filters"
# section 9 page 22.
#
# The coupling capacitors are shunt in this design.
#
import math
from sympy import symbols, Matrix, simplify, zeros, parse_expr, I, re, im
import matplotlib.pyplot as plt
import numpy as np
from network import Network
from filterdesign import *

# Setup parmaters
fc = 5000000
bw = 1000
Rs = 1000
N = 4
# Measured
Lm = 0.1
Qu_crystal = 100000
# Butterworth LPF parameters
g_list = butterworthNormalizedComponents(N)
# Convert into BPF normalized parameters
k_list = couplingCoefficientsButterworth(g_list)
k12, k23, k34 = k_list[0], k_list[1], k_list[0]
q = endSectionCoefficientButterworth(g_list)
# Angular center frequency
wc = 2 * math.pi * fc 
# Calculate motional capacitance Cm, given Lm
Cm = 1 / (wc**2 * Lm)
# The filter Q is a a function of the desired bandwidth
Qfilter = fc / bw 
# Compute the denormalized end Q (by definition, doesn't matter
# whether this is series or parallel)
Qe = 1.0 / ((1.0 / (q * Qfilter)) - (1.0 / Qu_crystal))
print("Qe", Qe)
# Compute the denormalized coupling capacitors. These 
# are shunt so the formulae are a bit different
C12 = Cm * Qfilter / k12
C23 = Cm * Qfilter / k23
C34 = Cm * Qfilter / k34
# Figure out the end resistance needed to properly load 
# the end resonator to acheive the desired end Q.  This 
# is by definition of "Q" for an inductor.  This is a series
# situation so the formulae are a bit different.
Rse = (1 / Qe) * wc * Lm
# Figure out what *parallel* end capacitor is needed to 
# acheive the necessary end loading, given the Rs that 
# we already have.
Ce = math.sqrt(( Rs - Rse ) / ( Rse * wc**2 * Rs**2))
# Figure out what the equivalent *series* capacitance
# is formed by the Rs and the parallel Ce.  This will be needed
# later during the re-tuning process.
Ce_series = (Ce**2 * wc**2 * Rs**2 + 1) / (Ce * wc**2 * Rs**2)

# Calculate the total capacitance in the end mesh
Cnet_1 = 1 / (1 / Ce_series + 1 / Cm + 1 / C12)
# Compute the resonant freqeuncy in the end mesh
f1 = 1 / (2 * math.pi * math.sqrt(Lm * Cnet_1))

# Calculate the total capacitance in the second mesh                                                                                                  
Cnet_2 = 1 / (1 / Cm + 1 / C12 + 1 / C23)
# Compute the resonant freqeuncy in the second mesh
f2 = 1 / (2 * math.pi * math.sqrt(Lm * Cnet_2))

# NOTICE: f1 and f2 are not the same!!
# Since f1 is higher, we will "tune up" f2 (and f3) so
# that it matches f1 (and f4).  Figure out what capacitance
# needs to be added to compensate:

# Pick the highest frequency and we'll tune all meshes to that
fh = f1
wh = 2 * math.pi * fh
Ch = 1 / (wh**2 * Lm)

# Figure out what extra needs to be added to mesh 2 to re-tune,
# keeping in mind that Cm, C12, and C23 are already fixed.
Cadjust_2 = 1 / (1 / Ch - 1 / Cm - 1 / C12 - 1 / C23)
# Figure out what extra needs to be added to mesh 3 to re-tune
Cadjust_3 = 1 / (1 / Ch - 1 / Cm - 1 / C23 - 1 / C34)

# Another way to get the adjustment directly from the 
# existing mesh capacitiances:
Cadjust = (Cnet_1 * Cnet_2) / (Cnet_2 - Cnet_1)
print("Cadjust", Cadjust)

# Re-tuning: Now we need to adjust the resonator capacitor in order 
# maintain the same fc.  We are "backing out" the capacitance
# that was introduced in the loop as a consequence of the loading and/or
# coupling capacitors

# Display
print("End capacitor Ce        ", Ce)
print("Couping capacitor C12   ", C12)
print("Couping capacitor C23   ", C23)
print("Couping capacitor C34   ", C34)
print("Adjusting capacitor 2   ", Cadjust_2)
print("Adjusting capacitor 3   ", Cadjust_3)
