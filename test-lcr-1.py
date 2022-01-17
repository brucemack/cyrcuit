# A complete example of creating an RLC network (Low-Pass Filter)
# and simulating it across a range of frequecies.
#
# The example is taken from the Radio Frequency Design book 
# on page 48.
#
import math
from sympy import symbols, Matrix, simplify, zeros, parse_expr, I, re, im
import matplotlib.pyplot as plt
import numpy as np
from network import Network

# Test stuff
network = Network()

# A very simple circuit (From Radio Frequency Design (Hayward) pg. 48)
network.add_element("vin", "va", "rs + other")
network.add_element("va", "vout", "z1")
network.add_element("vout", "gnd", "z2")
network.add_element("vout", "gnd", "rl")
# Shut off KCL for this node
network.set_input("vin")

# Get the solution for the node voltages
x, names = network.get_solution()
print(names)

# Setup the complex impedances using RLC values
s, w, l, c, r = symbols("s w l c r")

z_values = [
    (symbols("rs"), r),
    (symbols("z1"), l * s),
    (symbols("z2"), 1.0 / (c * s)),
    (symbols("rl"), r),
    (symbols("other"), 0)
]
x = x.subs(z_values)

# Plug in the actual values.  These are all normlized to a 1.0 impedance
# (symertric) and a 1.0 rad/sec cut-off frequency.
z_values = [
    (symbols("r"), 1.0),
    (symbols("l"), math.sqrt(2.0)),
    (symbols("c"), math.sqrt(2.0)),
]
x = x.subs(z_values)

# Change s->jw and simplify.  Notice the use of I (imaginary component)
x = simplify(x.subs(s, w * I))

# Create the transfer function H(jw) = vout(jw) / vin(jw)
H = x[3]

# Create the sweep of frequencies in rad/sec
input_angles = np.linspace(0.01, 5, 100)
h, w = symbols("h w")

# Used to evaluate the transfer function at a specific frequency
def eval_h(x):
    # Evaluate the vout expression with the specified frequency
    h = H.subs(w, x)
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
plt.xlabel('Frequency rad/sec')
plt.grid(color='grey', linestyle='dashed', linewidth=1)
plt.title("RFD Low Pass Filter Example - pg. 48")
plt.show()
