# A complete example of creating a four crystal filter network
# and simulating it across a range of frequecies.
#
# The example is taken from EMRFD on page 3.19.
#
import math
from sympy import symbols, Matrix, simplify, zeros, parse_expr, I, re, im, lambdify
import matplotlib.pyplot as plt
import numpy as np
from network import Network

print("Starting ...")
"""
 Vin-Zs-Va-Zc-Vb-Z3-Vc-Zc-Vd-Z5-Ve-Zc-Vf-Zc-Vout
         |     |           |           |    |  | 
         Z1    Z2          Z4          Z6   Z7 Zl
         |     |           |           |    |  |
         =     =           =           =    =  =

Zc represents the Crystal
"""
# Test stuff
network = Network()

network.add_element("vin", "va", "zs")
network.add_element("va", "gnd", "z1")
# Crystal
network.add_element("va", "vb", "zc")
network.add_element("vb", "gnd", "z2")
network.add_element("vb", "vc", "z3")
# Crystal 2
network.add_element("vc", "vd", "zc")
network.add_element("vd", "gnd", "z4")
network.add_element("vd", "ve", "z5")
# Crystal 3
network.add_element("ve", "vf", "zc")
network.add_element("vf", "gnd", "z6")
# Crystal 4
network.add_element("vf", "vout", "zc")
network.add_element("vout", "gnd", "z7")
network.add_element("vout", "gnd", "zl")
# Shut off KCL for this node
network.set_input("vin")

# Get the solution for the node voltages
print("Solving linear system ...")
x, names = network.get_solution()
print(names)

# Setup the complex impedances using RLC values
s, w = symbols("s w")
rs, rl, c1, c2, c3, c4, c5, c6, c7, rc, cc, lc, cpc = symbols("rs rl c1 c2 c3 c4 c5 c6 c7 rc cc lc cpc")

z_values = [
    # Crystals are a parallel combination of (rc, lc, and cc) and (cpc)
    #(symbols("zc"), 1.0 / (1.0 / (rc + s*lc + 1.0 / (s*cc)) + 1.0 / (1.0 / (s*cpc))) ),
    (symbols("zc"), (rc + s*lc + 1.0 / (s*cc)) ),
    # Caps
    (symbols("z1"), 1.0 / (s*c1)),
    (symbols("z2"), 1.0 / (s*c2)),
    (symbols("z3"), 1.0 / (s*c3)),
    (symbols("z4"), 1.0 / (s*c4)),
    (symbols("z5"), 1.0 / (s*c5)),
    (symbols("z6"), 1.0 / (s*c6)),
    (symbols("z7"), 1.0 / (s*c7)),
    # Source/load
    (symbols("zs"), rs),
    (symbols("zl"), rl),
]
x = x.subs(z_values)

# Quality factor (unloaded) for the crystal.  Modeling this 
# using a series resistor
quc = 240000
# Center freq and bandwidth
f0 = 5
bw = 400
# The ESR of the filter
rc = (1.2e8 * f0) / (bw * quc)

# Plug in the actual values.
lcr_values = [
    # Source/Load
    (symbols("rs"), 50 * 9),
    (symbols("rl"), 50 * 9),
    # Crystal
    # ESR of the crystal (computed from Qu)
    (symbols("rc"), rc),
    (symbols("lc"), 0.098),
    (symbols("cc"), 0.010339e-12),
    # Parallel capacitance measured
    (symbols("cpc"), 3e-12),
    # Caps
    (symbols("c1"), 47e-12),
    (symbols("c2"), 154e-12),
    (symbols("c3"), 328e-12),
    (symbols("c4"), 288e-12),
    (symbols("c5"), 328e-12),
    (symbols("c6"), 154e-12),
    (symbols("c7"), 47e-12),
]
x = x.subs(lcr_values)

# Change s->jw and simplify.  Notice the use of I (imaginary component)
x = x.subs(s, w * I)

# Create the transfer function H(jw) = vout(jw) / vin(jw)
# But we are assuming vin(jw) = 1.0
H = x[8]

# Create the sweep of frequencies in rad/sec
bw = 2000
fc = 5000000
input_angles = np.linspace(fc - bw, fc + bw, 100)
h, w = symbols("h w")

# Make a quicker version of the function
h_fast = lambdify(w, H)

# Used to evaluate the transfer function at a specific frequency
def eval_h(f):
    # Convert the frequency in Hz to frequency in Radias/second 
    omega = f * 2.0 * np.pi
    # Evaluate the vout expression with the specified frequency
    return h_fast(omega)

# Evaluate the transfer function across the phase range
print("Evaluating ...")
result_complex = np.vectorize(eval_h)(input_angles)
print("Done")
# This is a power ratio, so use 10*log10() for the dB calculation.
# Also notice the 4.0*(Vout**2) formulation since we are using available
# power for a matched load.
result_mag_db = 10.0 * np.log10(4.0 * (np.absolute(result_complex) ** 2.0))
# Linear result
#result_mag_db = np.absolute(result_complex)
# Simple phase calculation
#result_phase = np.angle(result_complex, True)

# matplotlib stuff
plt.plot(input_angles, result_mag_db)
plt.ylabel('H(f)')
plt.xlabel('Frequency (Hz)')
plt.xscale("log")
plt.grid(color='grey', linestyle='dashed', linewidth=1)
plt.title("EMRFD Crystal Filter - pg. 3.19")
plt.show()
