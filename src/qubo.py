import itertools
from collections import defaultdict
import math

from qiskit.quantum_info import SparsePauliOp

def build_qubo(R, S, mapping, adjlist, seqH, A = 10.0, B = 10.0, C = 10.0):
    linear = defaultdict(float)
    quadrs = defaultdict(float)
    const_ = 0.0


    def add_linear(i, coeff): linear[i] += coeff
    def add_quadratic(i, j, coeff): quadrs[tuple(sorted((i, j)))] += coeff

    def site_value(r):
        bits = mapping[r]
        weights = [2**k for k in range(len(bits))]
        return weights, bits
    
    def add_adj_penalty(r1, r2, coeff):
        w1, b1 = site_value(r1)
        w2, b2 = site_value(r2)

        for i, wi in enumerate(w1):
            for j, wj in enumerate(w2):
                add_linear(b1[i], coeff * wi**2)
                add_linear(b2[j], coeff * wj**2)
                add_quadratic(b1[i], b2[j], -2 * coeff * wi * wj)

        for i, wi in enumerate(w1):
            add_linear(b1[i], -2 * coeff * wi)
        for j, wj in enumerate(w2):
            add_linear(b2[j], +2 * coeff * wj)

        nonlocal const_
        const_ += coeff

    def add_self_collision_penalty(r1, r2, coeff):
        w1, b1 = site_value(r1)
        w2, b2 = site_value(r2)

        for i, wi in enumerate(w1):
            for j, wj in enumerate(w2):
                add_linear(b1[i], coeff * wi ** 2)
                add_linear(b2[j], coeff * wj ** 2)
                add_quadratic(b1[i], b2[j], -2 * coeff * wi * wj)

    def add_contact_reward(r1, r2, coeff):
        nonlocal const_
        add_adj_penalty(r1, r2, -coeff)
        const_ += coeff

    for r in range(R - 1):
        add_adj_penalty(r, r + 1, C)
    
    for i in range(R):
        for j in range(i + 1, R):
            add_self_collision_penalty(i, j, B)

    for i in range(R):
        for j in range(i + 1, R):
            if seqH[i] and seqH[j]:
                add_contact_reward(i, j, A)

    return linear, quadrs, const_

def qubo_to_pauli(linear, quadratic, constant, num_qubits):
    pauli_map = {}

    def add_pauli(pauli_label, coeff):
        pauli_map[pauli_label] = pauli_map.get(pauli_label, 0.0) + coeff

    for i, c in linear.items():
        # Global Term
        add_pauli('I' * num_qubits, c * 0.5) 

        # Z-Term
        label = ['I'] * num_qubits
        label[num_qubits - 1 - i] = 'Z'
        add_pauli(''.join(label), -0.5 * c)

    for (i, j), c in quadratic.items():
        # Global Term
        add_pauli('I'* num_qubits, 0.25 * c)

        # - 1/4 Z_i
        lab_i = ['I'] * num_qubits
        lab_i[num_qubits - 1 - i] = 'Z'
        add_pauli(''.join(lab_i), -0.25 * c)

        # - 1/4 Z_j
        lab_j = ['I'] * num_qubits
        lab_j[num_qubits - 1 - j] = 'Z'
        add_pauli(''.join(lab_j), -0.25 * c)

        # + 1/4 Z_i Z_j
        lab_ij = ['I'] * num_qubits
        lab_ij[num_qubits - 1 - i] = 'Z'
        lab_ij[num_qubits - 1 - j] = 'Z'
        add_pauli(''.join(lab_ij), 0.25 * c)


    pauli_list = [(k, v) for k, v in pauli_map.items() if abs(v) > 1e-12]
    
    max_coeff = max([abs(v) for _, v in pauli_list])
    if max_coeff > 1e-8:
        scale = 1.0/max_coeff
        pauli_list = [(p, v * scale) for p, v in pauli_list]
        constant *= scale

    op = SparsePauliOp.from_list(pauli_list)
    return op, constant