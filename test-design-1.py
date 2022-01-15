# Demonstration of the design of doubly-tuned-circuit.
# Source: Hayward "Designing Narrow-Bandwidth Ladder Filters"
# page 9.
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
# This can be chosen
L0 = 3e-6
# As measured on the inductor via lab test
Qu_inductor = 200
# Butterworth LPF parameters
g_list = butterworthNormalizedComponents(2)
# Convert into BPF normalized parameters
k_list = couplingCoefficientsButterworth(g_list)
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
# Compute the denormalized coupling capacitor
C12 = C0 * k_list[0] / Qfilter
# Figure out the end resistance needed to properly load 
# the end resonator to acheive the desired end Q.  This 
# is by definition of "Q" for an inductor.
Rpe = Qe * wc * L0
# Figure out what *series* end capacitor is needed to 
# acheive the necessary end loading, given the Rs that 
# we already have.
Ce = 1 / (wc * math.sqrt(Rpe * Rs - Rs ** 2))
# Figure out what the equivalent parallel (shunt) capacitance
# is formed by the Rpe and the series Ce.  This will be needed
# later during the re-tuning process.
Ce_parallel = Ce / ((Rpe * Ce) ** 2 + 1)
# Re-tuning: Now we need to adjust the resonator capacitor in order 
# maintain the same fc.  We are "backing out" the capacitance
# that was introduced in the loop 
Cretuned = C0 - Ce_parallel - C12
# Compute the insertion loss of the filter
il = 20 * math.log10(q0 / (q0 - q))
# Display
print("Resonator inductor    ", L0)
print("Resonator resistor    ", Rpe)
print("Resonator capacitor   ", Cretuned)
print("End capacitor Ce      ", Ce)
print("Couping capacitor C12 ", C12)
print("Insertion loss (db)    ", il)
