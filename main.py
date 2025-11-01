import itertools

from qiskit.quantum_info import SparsePauliOp
from qiskit_aer import Aer
from qiskit_algorithms.minimum_eigensolvers import VQE
from qiskit_algorithms.optimizers import COBYLA
from qiskit.circuit.library import TwoLocal
from qiskit_aer.primitives import EstimatorV2
def build_lattice(nx: int, ny: int, nz: int):
    # will add a type param later, so that it can build FCC/BCC/etc.
    coords = []
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                coords.append((x, y, z))

    return coords

def build_adj_list(coords):
    # arhamaneeq/sw-pathfinders.git wasnt useless
    S = len(coords)
    adj = {
        i: set()
        for i in range(S)
    }

    for i, a in enumerate(coords):
        for j, b in enumerate(coords):
            if i == j:
                continue
            if sum(abs(a[k] - b[k]) for k in range(3)) == 1:
                adj[i].add(j)

    return adj

def map_idx_to_qubit(R, S):
    # so R is the residue number, i.e., the Rth (0-indexed) amino acid in the chain
    # row-maj flattening from (r,p) -> idx
    # grouped by residue
    
    mapping = {}

    for r in range(R):
        for p in range(S):
            mapping[(r, p)] = r * S + p

    return mapping

# QUBOOOOOOOOOOOOOOO
def build_qubo(R, S, mapping, adjlist, seqH, A = 10.0, B = 10.0, C = 10.0):
    linear = {}
    quadrs = {}
    constant = 0.0

    def add_linear(term_dict: dict, i, coeff):
        term_dict.setdefault((i,), 0.0)
        term_dict[(i,)] += coeff
    def add_quadratic(term_dict: dict, i, j, coeff):
        key = tuple(sorted((i, j)))
        term_dict.setdefault(key, 0.0)
        term_dict[key] += coeff


    for i in range(R):
        for j in range(i + 1, R):
            if abs(i - j) == 1 or not (seqH[i] and seqH[j]):
                continue
            else:
                Eij = -1.0
            for p in range(S):
                for q in adjlist[p]:
                    if q <= p:
                        continue
                    idx_i = mapping[(i, p)]
                    idx_j = mapping[(j, q)]
                    add_quadratic(quadrs, idx_i, idx_j, Eij)

    for r in range(R):
        ss = [mapping[(r, p)] for p in range(S)]

        for i in ss:
            add_linear(linear, i, A * 1.0)
        for (i, j) in itertools.combinations(ss, 2):
            add_quadratic(quadrs, i, j, 2.0 * A)

        constant += A * 1.0

        for i in ss:
            add_linear(linear, i, -2.0 * A)

    for p in range(S):
        rs = [mapping[(r, p)] for r in range(R)]
        for (i, j) in itertools.combinations(rs, 2):
            add_quadratic(quadrs, i, j, B)


    for r in range(R - 1):
        for p in range(S):
            idx = mapping[(r, p)]
            add_linear(linear, idx, C * 1.0)

        for p in range(S):
            idx_rp = mapping[(r, p)]
            for q in adjlist[p]:
                idx_r1q = mapping[(r + 1, q)]
                add_quadratic(quadrs, idx_rp, idx_r1q, -C * 1.0)

    return linear, quadrs, constant

def qubo_to_pauli(linear, quadratic, constant, num_qubits):
    pauli_map = {}

    def add_pauli(pauli_label, coeff):
        pauli_map[pauli_label] = pauli_map.get(pauli_label, 0.0) + coeff

    for (i, ), c in linear.items():
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
    labels = [p for p, _ in pauli_list]
    coeffs = [v for _, v in pauli_list]
    op = SparsePauliOp.from_list(list(zip(labels, coeffs)))

    return op, constant

# QAOAOAOAOAOAO
if __name__ == "__main__":
    nx, ny, nz = 2, 2, 1
    coords = build_lattice(nx,ny,nz)
    adjlist = build_adj_list(coords)
    S = len(coords)

    seqH = [True, False, True]
    R = len(seqH)

    mapping = map_idx_to_qubit(R,S)
    A, B, C = 5.0, 5.0, 5.0

    linear, quadratic, constant = build_qubo(R, S, mapping, adjlist, seqH, A, B, C)
    num_qubits = R * S

    pauli_op, const = qubo_to_pauli(linear, quadratic, constant, num_qubits)
    
    est = EstimatorV2()
    
    ansatz = TwoLocal(
        num_qubits, 
        ['ry', 'rz'], 'cz', 
        reps=2
    ).decompose().decompose()

    vqe = VQE(
        ansatz=ansatz, 
        optimizer=COBYLA(), 
        estimator=est
    )
    res = vqe.compute_minimum_eigenvalue(
        operator=pauli_op
    )
    print(res)