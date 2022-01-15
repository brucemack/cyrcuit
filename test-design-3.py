# Demonstration of the design of N=2 bandpass filter with 
# series LC tuned circuits.

# Source: Hayward "Designing Narrow-Bandwidth Ladder Filters"
# section 6 page 15.
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
bw = 200000
Rs = 50
N = 2
# This can be chosen
L0 = 7e-6
# As measured on the inductor via lab test
Qu_inductor = 200
# Butterworth LPF parameters
g_list = butterworthNormalizedComponents(N)
# Convert into BPF normalized parameters
k_list = couplingCoefficientsButterworth(g_list)
k12, k23 = k_list[0], k_list[0]
q = endSectionCoefficientButterworth(g_list)
# Angular center frequency
wc = 2 * math.pi * fc 
# Calculate nodal capacitance C0, given the selected L0
C0 = 1 / (wc ** 2 * L0)
# The filter Q is a a function of the desired bandwidth
Qfilter = fc / bw 
# Compute the normalized q0.  We're looking for something
# more than 2.
q0 = Qu_inductor / Qfilter
print("q0 (looking for >2)", q0)
# Compute the denormalized end Q (by definition, doesn't matter
# whether this is series or parallel)
Qe = 1.0 / ((1.0 / (q * Qfilter)) - (1.0 / Qu_inductor))
# Compute the denormalized coupling capacitors. These 
# are shunt so the formulae are a bit different
C12 = C0 * Qfilter / k12
C23 = C0 * Qfilter / k23
# Figure out the end resistance needed to properly load 
# the end resonator to acheive the desired end Q.  This 
# is by definition of "Q" for an inductor.  This is a series
# situation so the formulae are a bit different.
Rse = (1 / Qe) * wc * L0
# Figure out what *parallel* end capacitor is needed to 
# acheive the necessary end loading, given the Rs that 
# we already have.
Ce = math.sqrt(( Rs - Rse ) / ( Rse * wc * wc * Rs * Rs))
# Figure out what the equivalent series capacitance
# is formed by the Rs and the parallel Ce.  This will be needed
# later during the re-tuning process.
Ce_series = (Ce**2 * wc**2 * Rs**2 + 1) / (Ce * wc**2 * Rs**2)
# Re-tuning: Now we need to adjust the resonator capacitor in order 
# maintain the same fc.  We are "backing out" the capacitance
# that was introduced in the loop 
Cretuned_1 = 1 / (1 / C0 - 1 / Ce_series - 1 / C12)
Cretuned_2 = 1 / (1 / C0 - 1 / Ce_series - 1 / C23)

# Display
print("Resonator inductor      ", L0)
print("Resonator resistor      ", Rse)
print("Resonator capacitor 1   ", Cretuned_1)
print("Resonator capacitor 2   ", Cretuned_2)
print("End capacitor Ce        ", Ce)
print("Couping capacitor C12   ", C12)
print("Couping capacitor C23   ", C23)
