from qiskit_aer import Aer
from qiskit_algorithms.minimum_eigensolvers import VQE
from qiskit_algorithms.optimizers import COBYLA
from qiskit.circuit.library import TwoLocal
from qiskit_aer.primitives import EstimatorV2
from qiskit.quantum_info import Statevector

from .lattice import build_adj_list, build_lattice, map_idx_to_qubit, choose_lattice_dims, decode_positions, positions_to_coords, decode_bitstring
from .qubo import qubo_to_pauli, build_qubo
from .interpreter import classical_ground_state, visualize_folds
from .amino import is_hydrophobic

def Estimate(amino_seq):
    seqH = [is_hydrophobic(amine) for amine in list(amino_seq)]
    R = len(seqH)

    nx, ny, nz = choose_lattice_dims(R)
    coords = build_lattice(nx, ny, nz, "fcc")
    adjlist = build_adj_list(coords)
    S = len(coords)

    mapping, num_bits = map_idx_to_qubit(R,S)
    num_qubits = R * num_bits
    A, B, C, D = 5.0, 500.0, 5.0, 1.0

    linear, quadratic, constant = build_qubo(R, S, mapping, adjlist, seqH, A, B, C, D)
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

    opt_qc = res.optimal_circuit
    print(type(opt_qc))
    opt_params = res.optimal_parameters
    bound_qc = opt_qc.assign_parameters(opt_params)
    statevec = Statevector(bound_qc)
    probs = statevec.probabilities_dict()
    probs = sorted(probs.items(), key = lambda x: -x[1])
    
    for b, p in probs[:10]:
        print(f"{b}: {p:.6f}")

    decoded_confs = []
    for b, p in probs[:3]:
        if p > 0.1:
            decoded = decode_bitstring(b, mapping, coords)
            decoded_confs.append((b, p, decoded))

    print("\nDecoded lattice positions of top states:")
    for bitstring, prob, decoded in decoded_confs:
        print(f"State {bitstring} (p={prob:.4f}): {decoded}")

    visualize_folds(decoded_confs, seqH, amino_seq)


    # Classical Visualisation via Brute Force

    # bits, E = classical_ground_state(linear, quadratic)
    # positions = decode_positions(bits, mapping, S)
    # decoded_coords = positions_to_coords(positions, coords)
    # print(decoded_coords)
    # visualize_fold(decoded_coords, seqH)

