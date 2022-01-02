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

class Node:
    """ Observable node in the circuit """
    def __init__(self, name: str, input: bool, ground: bool):
        self.name = name 
        self.input = input
        self.ground = ground
        self.ordinal = 0

class Edge:
    """ The connection between two notes in the circuit """
    def __init__(self, start: Node, end: Node, imp: str):
        self.start = start
        self.end = end
        self.imp = imp

class Network:
    def __init__(self):
        self.edges = []
        self.nodes = {}
        # Automatically create ground
        n = self.get_or_create_node("gnd")
        n.ground = True

    def get_or_create_node(self, name: str):
        if name not in self.nodes:
            c = len(self.nodes)
            n = self.nodes[name] = Node(name, False, False)
            n.ordinal = c
        return self.nodes[name]

    def add_element(self, node0: str, node1: str, imp: str):
        n0 = self.get_or_create_node(node0)
        n1 = self.get_or_create_node(node1)
        self.edges.append(Edge(n0, n1, imp))

    def set_input(self, name: str):
        self.nodes[name].input = True

    def get_linear_system(self):

        # Create the system of equations based on the KCL for each node.
        # Start off with zeros
        A = Matrix(len(self.nodes), len(self.nodes), lambda x, y: 0)
        b = Matrix(len(self.nodes), 1, lambda x, y: 0)
        names = [None] * len(self.nodes)

        # Consider each node individually
        for node in self.nodes.values():
            names[node.ordinal] = node.name
            # Don't write an KCL expression for the input node
            if node.input == True:
                A[node.ordinal, node.ordinal] += 1
                b[node.ordinal] = 1
            # Don't write an KCL expression for the ground node
            elif node.ground == True:
                A[node.ordinal, node.ordinal] += 1
                b[node.ordinal] = 0
            # All other nodes need KCL
            else:
                # Sum of the currents is always zero
                b[node.ordinal] = 0
                # Look at all of the edges that touch the node.  Edges that 
                # touch add/subtract current.
                for edge in self.edges:
                    if edge.start == node or edge.end == node:
                        contrib = 1.0 / symbols(edge.imp)
                        # The direction of contribution depends on whether this 
                        # is an out-flowing branch or an in-flowing one.
                        if edge.end == node:
                            A[node.ordinal, edge.start.ordinal] -= contrib
                            A[node.ordinal, edge.end.ordinal] += contrib
                        else:
                            A[node.ordinal, edge.start.ordinal] += contrib
                            A[node.ordinal, edge.end.ordinal] -= contrib
        return A, b, names

# Test stuff
network = Network()

# A very simple circuit (From Radio Frequency Design (Hayward) pg. 48)
network.add_element("vin", "va", "rs")
network.add_element("va", "vout", "z1")
network.add_element("vout", "gnd", "z2")
network.add_element("vout", "gnd", "rl")
# Shut off KCL for this node
network.set_input("vin")

# Create the matrix representation of the circuit
A, b, names = network.get_linear_system()
# Solve for the node voltages
x = simplify(A.LUsolve(b))
print(names)

# Setup the complex impedances using RLC values
s, w, l, c, r = symbols("s w l c r")

z_values = [
    (symbols("rs"), r),
    (symbols("z1"), l  * s),
    (symbols("z2"), 1.0 / (c * s)),
    (symbols("rl"), r),
]
x = x.subs(z_values)

# Plug in the actual values.  These are all normlized to a 1.0 impedance
# (symertric) and a 1.0 rad/sec cut-off frequency.
z_values = [
    (symbols("r"), 1),
    (symbols("l"), math.sqrt(2.0)),
    (symbols("c"), math.sqrt(2.0)),
]
x = x.subs(z_values)

# Change s->jw and simplify.  Notice the use of I (imaginary component)
x = simplify(x.subs(s, w * I))

# Create the transfer function H(jw) = vout(jw) / vin(jw)
H = x[3]

# Create the sweep of frequencies
input_angles = np.linspace(0, 5, 100)
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
plt.grid(color='grey', linestyle='dashed', linewidth=1)
plt.show()
