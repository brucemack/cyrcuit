# Demonstration of the design of N=4 bandpass filter with 
# parallel LC tuned circuits.

# Source: Hayward "Designing Narrow-Bandwidth Ladder Filters"
# page 11.
#
import math
from sympy import symbols, Matrix, simplify, zeros, parse_expr, I, re, im
import matplotlib.pyplot as plt
import numpy as np
from network import Network
from filterdesign import *

# Twp parallel resonators with capacitor coupler

# Setup parmaters
fc = 5000000
bw = 200000
Rs = 50
N = 4
# This can be chosen
L0 = 6e-6
# As measured on the inductor via lab test
Qu_inductor = 200
# Butterworth LPF parameters
g_list = butterworthNormalizedComponents(N)
# Convert into BPF normalized parameters
k_list = couplingCoefficientsButterworth(g_list)
k12, k23, k34 = k_list[0], k_list[1], k_list[0]
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
# Compute the denormalized end Q
Qe = 1.0 / ((1.0 / (q * Qfilter)) - (1.0 / Qu_inductor))
# Compute the denormalized coupling capacitors
C12 = C0 * k12 / Qfilter
C23 = C0 * k23 / Qfilter
C34 = C0 * k34 / Qfilter
# Figure out the end resistance needed to properly load 
# the end resonator to acheive the desired end Q.  This 
# is by definition of "Q" for an inductor.
Rpe = Qe * wc * L0
# Figure out what *series* end capacitor is needed to 
# acheive the necessary end loading, given the Rs that 
# we already have.
Ce = 1 / (wc * math.sqrt(Rpe * Rs - Rs ** 2))
# Figure out what the equivalent parallel (shunt) capacitance
# is formed by the Rs and the series Ce.  This will be needed
# later during the re-tuning process.
Ce_parallel = Ce / ((Rs * wc * Ce) ** 2 + 1)
# Re-tuning: Now we need to adjust the resonator capacitor in order 
# maintain the same fc.  We are "backing out" the capacitance
# that was introduced in the loop 
Cretuned_1 = C0 - Ce_parallel - C12
Cretuned_2 = C0 - C12 - C23
Cretuned_3 = C0 - C23 - C34
Cretuned_4 = C0 - C34 - Ce_parallel

# Display
print("Resonator inductor      ", L0)
print("Resonator resistor      ", Rpe)
print("Resonator capacitor 1   ", Cretuned_1)
print("Resonator capacitor 2   ", Cretuned_2)
print("Resonator capacitor 3   ", Cretuned_3)
print("Resonator capacitor 4   ", Cretuned_4)
print("End capacitor Ce        ", Ce)
print("Couping capacitor C12   ", C12)
print("Couping capacitor C23   ", C23)
print("Couping capacitor C34   ", C34)
