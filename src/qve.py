from qiskit_aer import Aer
from qiskit_algorithms.minimum_eigensolvers import VQE
from qiskit_algorithms.optimizers import COBYLA
from qiskit.circuit.library import TwoLocal
from qiskit_aer.primitives import EstimatorV2

from .lattice import build_adj_list, build_lattice, map_idx_to_qubit, choose_lattice_dims
from .qubo import qubo_to_pauli, build_qubo

def Estimate(seqH):
    R = len(seqH)

    nx, ny, nz = choose_lattice_dims(R)
    coords = build_lattice(nx, ny, nz, "fcc")
    adjlist = build_adj_list(coords)
    S = len(coords)

    mapping, num_bits = map_idx_to_qubit(R,S)
    num_qubits = R * num_bits
    A, B, C = 5.0, 5.0, 5.0

    linear, quadratic, constant = build_qubo(R, S, mapping, adjlist, seqH, A, B, C)
    print(f"R:{R} | S:{S} | N_q: {num_qubits}\nQubit Density: {R * S / num_bits}")

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