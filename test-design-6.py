# Demonstration of the design of a generalized filter 
# using an arbitrary number of crystals >2.
#
# The coupling capacitors are shunt in this design.
#
import math
from sympy import symbols, Matrix, simplify, zeros, parse_expr, I, re, im, lambdify
import matplotlib.pyplot as plt
import numpy as np
from network import Network
from filterdesign import *

def simulate_filter(N, Lm, Cm, Cp, Rse, Cs_list, Ck_list, use_cp):

    # Test stuff
    network = Network()

    # Source impedance
    network.add_element("vin", "v1", "rs")
    # Shut off KCL for the input node
    network.set_input("vin")
    # Input Mesh
    network.add_element("v1", "v2", "zs1 + zx")
    network.add_element("v2", "gnd", "zk1")
    # Output Mesh with load impedance
    network.add_element("v" + str(N), "vout", "zs" + str(N) + " + zx")
    network.add_element("vout", "gnd", "rl")
    # Go through the inner meshes (stops before output mesh)
    for mesh in range(2, N):
        network.add_element("v" + str(mesh), "v" + str(mesh+1), "zs" + str(mesh) + " + zx")
        network.add_element("v" + str(mesh + 1), "gnd", "zk" + str(mesh))

    # Get the solution for the node voltages
    x, names = network.get_solution()
    print(names)

    s, w, rs, rl = symbols("s w rs rl")

    def parallel(z1, z2):
        return 1 / ((1 / z1) + (1 / z2))

    z_values = []
    # Setup the source/load impedances
    z_values.append((symbols("rs"), Rse))
    z_values.append((symbols("rl"), Rse))
    # Setup the crystal values
    if not use_cp:
        z_values.append((symbols("zx"), s * Lm + 1.0 / (s * Cm) + Rx))
    else:
        z_values.append((symbols("zx"), parallel(s * Lm + 1.0 / (s * Cm) + Rx, (1.0 / (s * Cp)))))
    # Setup the capactors for each mesh, except for the end
    for mesh in range(1, N):
        cs = Cs_list[mesh-1]
        if cs == 0:
            # The no capacitor case looks like a short-circuit
            z_values.append((symbols("zs" + str(mesh)), 0))
        else:
            z_values.append((symbols("zs" + str(mesh)), 1.0 / (cs * s)))
        ck = Ck_list[mesh-1]
        z_values.append((symbols("zk" + str(mesh)), 1.0 / (ck * s)))
    # Output mesh
    cs = Cs_list[N-1]
    z_values.append((symbols("zs" + str(N)), 1.0 / (cs * s)))

    x = x.subs(z_values)

    # Change s->jw.  Notice the use of I (imaginary component)
    x = x.subs(s, w * I)

    return x

# Setup parameters
fc = 5000000
bw = 3000
N = 4
# Measured
Lm = 0.098
Cm = 0.010339e-12
Cp = 4e-12
Qu_crystal = 240000
# The ESR of the filter
Rx = (1.2e8 * (fc / 1000000)) / (bw * Qu_crystal)

# Angular center frequency
wc = 2 * math.pi * fc 
# Calculate motional capacitance Cm, given Lm
#Cm = 1 / (wc**2 * Lm)
# The filter Q is a a function of the desired bandwidth
Qfilter = fc / bw 

# Butterworth LPF parameters
g_list = butterworthNormalizedComponents(N)
# Convert into BPF normalized parameters
k_list = couplingCoefficientsButterworth(g_list)
q = endSectionCoefficientButterworth(g_list)

# Compute the denormalized end Q (by definition, doesn't matter
# whether this is series or parallel)
Qe = 1.0 / ((1.0 / (q * Qfilter)) - (1.0 / Qu_crystal))

# Figure out the end resistance needed to properly load 
# the end resonator to acheive the desired end Q.  This 
# is by definition of "Q" for an inductor.
#
# THIS HAS NOTING TO DO WITH THE DESIRED SYSTEM IMPEDANCE R0!
# THIS IS WHAT IS REQUIRED TO HIT THE END Qe.
Rse = (1 / Qe) * wc * Lm

# Produce the denormalized coupling capacitors
Ck_list = [Cm * Qfilter / x for x in k_list]

# Compute the total capacitance of the second mesh since this will
# be used to establish the frequency of resonance across all of the
# meshes.  This works becase the second mesh has the highest 
# frequency because the crystal is connected in series with the 
# smallest coupling capacitances on both sides.

# The second mesh has the first and second shunt coupling capactors
# and the motional capacitance of the crystal all connected in 
# series.
Cmesh2 = 1 / ( 1 / Cm + 1 / Ck_list[0] + 1 / Ck_list[1])
# Convert the second mesh capacitance to the network mesh frequency
fmesh2 = 1 / (2 * math.pi * math.sqrt(Lm * Cmesh2))

# Walk through the meshes and determine the appropriate tuning 
# capacitor that is needed to make the mesh resonate at the 
# same frequency as mesh #2.
Cs_list = []

for i in range(0, N):
    # Special case for the first mesh (an end)
    if i == 0:
        # Adjustment is whatever is neededd to match the loop
        # capacitance of the second mesh.  We pull out the series 
        # equivalent of the load match, the motional capacitance 
        # of the crystal, and the coupling capacitor. 
        Cs = 1 / (1 / Cmesh2 - 1 / Cm - 1 / Ck_list[i])
    # Special case for the last mesh (an end)
    elif i == N-1:
        # See above case
        Cs = 1 / (1 / Cmesh2 - 1 / Cm - 1 / Ck_list[i-1])
    # Special case for the second mesh (and the mirror) since
    # they don't require an adjustment (alread the higest freq)
    elif i == 1 or i == N-2:
        Cs = 0
    else:
        # Here we pull out the motional capacitance and the two coupling 
        # capacitances
        Cs = 1 / (1 / Cmesh2 - 1 / Cm - 1 / Ck_list[i-1] - 1 / Ck_list[i])
    
    Cs_list.append(Cs)    

# Re-tuning: Now we need to adjust the resonator capacitor in order 
# maintain the same fc.  We are "backing out" the capacitance
# that was introduced in the loop as a consequence of the loading and/or
# coupling capacitors

# Display
print("Computed Cm                  ", Cm)
print("Required load resistance     ", Rse)
print("Center frequency             ", fmesh2)
print("Couping capacitors           ", Ck_list)
print("Adjusting capacitors         ", Cs_list)

# Create the sweep of frequencies in rad/sec
input_angles = np.linspace(fc - bw * 2, fc + bw * 2, 200)

# ----- First Curve -----

# Create the system of equations for the filter
x_a = simulate_filter(N, Lm, Cm, Cp, Rse, Cs_list, Ck_list, False)

# Create the transfer function H(jw) = vout(jw) / vin(jw)
# But we are assuming vin(jw) = 1.0
H_a = x_a[5]

# Make a quicker version of the function
h_fast_a = lambdify(symbols("w"), H_a)

# Used to evaluate the transfer function at a specific frequency
def eval_h_a(f):
    # Convert the frequency in Hz to frequency in Radias/second 
    omega = f * 2.0 * np.pi
    # Evaluate the vout expression with the specified frequency
    return h_fast_a(omega)

# Evaluate the transfer function across the phase range
print("Evaluating ...")
result_complex_a = np.vectorize(eval_h_a)(input_angles)
print("Done")
# This is a power ratio, so use 10*log10() for the dB calculation.
# Also notice the 4.0*(Vout**2) formulation since we are using available
# power for a matched load.
result_mag_db_a = 10.0 * np.log10(4.0 * (np.absolute(result_complex_a) ** 2.0))

# ----- Second Curve -----

# Create the system of equations for the filter
x_b = simulate_filter(N, Lm, Cm, Cp, Rse, Cs_list, Ck_list, True)

# Create the transfer function H(jw) = vout(jw) / vin(jw)
# But we are assuming vin(jw) = 1.0
H_b = x_b[5]

# Make a quicker version of the function
h_fast_b = lambdify(symbols("w"), H_b)

# Used to evaluate the transfer function at a specific frequency
def eval_h_b(f):
    # Convert the frequency in Hz to frequency in Radias/second 
    omega = f * 2.0 * np.pi
    # Evaluate the vout expression with the specified frequency
    return h_fast_b(omega)

# Evaluate the transfer function across the phase range
result_complex_b = np.vectorize(eval_h_b)(input_angles)
# This is a power ratio, so use 10*log10() for the dB calculation.
# Also notice the 4.0*(Vout**2) formulation since we are using available
# power for a matched load.
result_mag_db_b = 10.0 * np.log10(4.0 * (np.absolute(result_complex_b) ** 2.0))

# Figure out the max
peak = np.amax(result_mag_db_a)
print("Peak",peak)

# Look for the places where we cross -3dB
in_passband = False
min_found = False
max_found = False

for i in range(0, len(result_complex_a)):
    f = input_angles[i]
    mag = result_mag_db_a[i]
    if not in_passband:
        if mag - peak > -3:
            in_passband = True
            min_found = True
            min_mag = mag
            min_freq = f
    else:
        if mag - peak < -3:
            in_passband = False
            max_found = True
            max_mag = mag
            max_freq = f
            break

if min_found and max_found:
    width_3db = max_freq - min_freq
    center_3db = (max_freq + min_freq) / 2
    print("Min", min_freq)
    print("Max", max_freq)
    print("Center", center_3db)
    print("3dB Bandwidth", width_3db)

# matplotlib stuff
plt.plot(input_angles, result_mag_db_a)
plt.plot(input_angles, result_mag_db_b, 'r')
if min_found:
    plt.plot(min_freq, min_mag, 'bo')
if max_found:
    plt.plot(max_freq, max_mag, 'bo')
plt.ylabel('H(f)')
plt.xlabel('Frequency (Hz)')
plt.xscale("log")
plt.grid(color='grey', linestyle='dashed', linewidth=1)
if min_found and max_found:
    plt.title("Filter center: " + str(int(center_3db)) + " width: " + str(int(width_3db)))
plt.show()
