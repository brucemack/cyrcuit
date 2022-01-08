from sympy import symbols, Matrix, simplify, zeros, parse_expr, I, re, im
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
                A[node.ordinal, node.ordinal] += 1.0
                b[node.ordinal] = 1.0
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

    def get_solution(self):
        a, b, names = self.get_linear_system()
        # Solve for the node voltages
        x = a.LUsolve(b)
        return x, names
