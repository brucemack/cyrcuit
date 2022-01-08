# A complete example of creating an DTC network (Band-Pass Filter)
# and simulating it across a range of frequecies.
#
# The example is taken from EMRFD on page 3.14, table 3A, frist row.
#
import math
from sympy import symbols, Matrix, simplify, zeros, parse_expr, I, re, im
import matplotlib.pyplot as plt
import numpy as np
from network import Network

print("Starting ...")

# Test stuff
network = Network()

# A DTC
network.add_element("vin", "va", "rs")
network.add_element("va", "vb", "z1")
network.add_element("vb", "gnd", "z2")
network.add_element("vb", "gnd", "z3")
network.add_element("vb", "vc", "z4")
network.add_element("vc", "gnd", "z5")
network.add_element("vc", "gnd", "z6")
network.add_element("vc", "vout", "z7")
network.add_element("vout", "gnd", "rl")
# Shut off KCL for this node
network.set_input("vin")

# Get the solution for the node voltages
print("Solving linear system ...")
x, names = network.get_solution()
print(names)

# Setup the complex impedances using RLC values
s, w, l_qu, rs, c1, l2, c3, c4, l5, c6, c7, rl = symbols("s w l_qu rs c1 l2 c3 c4 l5 c6 c7 rl")

# Quality factor (unloaded) for the inductor.  Modeling this 
# using a series resistor
l_qu = 200

z_values = [
    (symbols("z1"), 1.0 / (c1 * s)),
    # Notice that the Q-driven resistor depends on w (not s!)
    (symbols("z2"), (l2 * s) + (l2 * w / l_qu)),
    (symbols("z3"), 1.0 / (c3 * s)),
    (symbols("z4"), 1.0 / (c4 * s)),
    (symbols("z5"), (l5 * s) + (l5 * w / l_qu)),
    (symbols("z6"), 1.0 / (c6 * s)),
    (symbols("z7"), 1.0 / (c7 * s))
]
x = x.subs(z_values)

# Plug in the actual values.
lcr_values = [
    # Source Z
    (symbols("rs"), 50),
    # Ce
    (symbols("c1"), 250e-12),
    # First LC
    (symbols("l2"), 6.98e-6),
    (symbols("c3"), 775e-12),
    # Coupling12
    (symbols("c4"), 41e-12),
    # Second LC
    (symbols("l5"), 6.98e-6),
    (symbols("c6"), 775e-12),
    # Ce
    (symbols("c7"), 250e-12),
    # Load Z
    (symbols("rl"), 50)
]
x = x.subs(lcr_values)

# Change s->jw and simplify.  Notice the use of I (imaginary component)
x = simplify(x.subs(s, w * I))

# Create the transfer function H(jw) = vout(jw) / vin(jw)
# But we are assuming vin(jw) = 1.0
H = x[5]

# Create the sweep of frequencies in rad/sec
input_angles = np.linspace(1300000, 2100000, 50)
h, w = symbols("h w")

# Used to evaluate the transfer function at a specific frequency
def eval_h(f):
    # Convert the frequency in Hz to frequency in Radias/second 
    phi = f * 2.0 * np.pi
    # Evaluate the vout expression with the specified frequency
    h = H.subs(w, phi)
    # Return as a native Python complex value
    return complex(re(h), im(h))

# Evaluate the transfer function across the phase range
print("Evaluating ...")
result_complex = np.vectorize(eval_h)(input_angles)
print("Done")
# This is a power ratio, so use 10*log10() for the dB calculation.
# Also notice the 4.0*(Vout**2) formulation since we are using available
# power for a matched load.
result_mag_db = 10.0 * np.log10(4.0 * (np.absolute(result_complex) ** 2.0))
# Simple phase calculation
result_phase = np.angle(result_complex, True)

# matplotlib stuff
plt.plot(input_angles, result_mag_db)
plt.ylabel('H(f)')
plt.xlabel('Frequency (Hz)')
plt.xscale("log")
plt.grid(color='grey', linestyle='dashed', linewidth=1)
plt.title("EMRFD DTC - pg. 3.14")
plt.show()
