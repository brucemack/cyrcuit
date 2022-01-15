# A complete example of creating an RLC network (Low-Pass Filter)
# and simulating it across a range of frequecies.
#
# The example is taken from EMRFD on page 3.4
#
import math
from sympy import symbols, Matrix, simplify, zeros, parse_expr, I, re, im
import matplotlib.pyplot as plt
import numpy as np
from network import Network

# Test stuff
network = Network()

# A very simple circuit EMRFD page 3.4
network.add_element("vin", "va", "rs")
network.add_element("va", "vb", "z1")
network.add_element("vb", "gnd", "z2")
network.add_element("vb", "vout", "z3")
network.add_element("vout", "gnd", "z4")
network.add_element("vout", "gnd", "rl")
# Shut off KCL for this node
network.set_input("vin")

# Get the solution for the node voltages
x, names = network.get_solution()
print(names)

# Setup the complex impedances using RLC values
s, w, l1, c2, l3, c4, r = symbols("s w l1 c2 l3 c4 r")

z_values = [
    (symbols("rs"), r),
    (symbols("z1"), l1 * s),
    (symbols("z2"), 1.0 / (c2 * s)),
    (symbols("z3"), l3 * s),
    (symbols("z4"), 1.0 / (c4 * s)),
    (symbols("rl"), r)
]
x = x.subs(z_values)

# Plug in the actual values.
lcr_values = [
    (symbols("r"), 50),
    (symbols("l1"), 0.609e-6),
    (symbols("c2"), 580e-12),
    (symbols("l3"), 1.472e-6),
    (symbols("c4"), 244e-12)
]
x = x.subs(lcr_values)

# Change s->jw and simplify.  Notice the use of I (imaginary component)
x = simplify(x.subs(s, w * I))

# Create the transfer function H(jw) = vout(jw) / vin(jw)
# But we are assuming vin(jw) = 1.0
H = x[4]

# Create the sweep of frequencies in rad/sec
input_angles = np.linspace(0, 50000000, 100)
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
result_complex = np.vectorize(eval_h)(input_angles)
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
plt.title("EMRFD Low Pass Filter Example - pg. 3.4")
plt.show()
