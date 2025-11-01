> [!NOTE]
> This repository was built as a submission for the *IISC Quantum Fall Fest 2025 | Cleveland Clinic QFF Hackathon Prompt*

> [!IMPORTANT]
> This README.md contains MathJax. Best viewed on GitHub on a Desktop Browser

# EigenFold

## The Plan (for now)

### QUBO

We use QUBO, i.e., for variables $b_i \in \{0, 1\}$ we have

$$
E(b_1, b_2, ..., b_N) = \sum_i a_i b_i + \sum_{i \lt j} b_i Q_{ij} b_j + c
$$

#### One-Hot Encoding in QUBO

We can't have hard constraints in QUBO, but we can penalties if the constraint is violated. Suppose we have a one-hot encoded variable, then

$$
\sum_p b_{r,p} - 1 =0
$$

which we can turn into a squared penalty

$$
A \times \Big(\sum_p b_{r,p} - 1\Big)^2
$$

where $A$ is a large penalty. Now, keeping in mind that $b^2_{r,p} = b_{r, p}$, since $b \in \{0, 1\}$, we can express the penalty as a sum of a linear, quadratic, and constant term.

$$
A \times \Big(\sum_p b_{r,p} - 1\Big)^2 = A \times \Big( -\sum_p b_{r,p} + 2 \sum_{p \lt q} b_{r,p}b_{r, q} + 1\Big)
$$

The $2A\sum_{p\lt q} b_{r,p}b_{r,q}$ term discourages residue duplication, i.e., no residue can occupy multiple sites at once, while the $-A\sum_p b_{r,p}$ term prevents residue disappearance, i.e., each residue must appear at least once.

We must also encode the fact that sequential residues must be placed next to each other in the lattice, i.e., enforce a rule backbone adjacency rule. This is achieved by the term $-C\sum_p \sum_{q\in N(p)} b_{r, p} b_{r+1, q}$.

#### Mapping QUBO Terms to Pauli Words

The QUBO energy function, where $c$ is a constant offset, $a_i$ is the linear coefficient for $b_i$, and $Q_{ij}$ is the quadratic coupling betweens bits $i$ and $j$, is of the form 

$$
E(\bold{b}) = c + \sum_i a_i b_i + \sum_{i \lt j} Q_{ij} b_i b_j
$$

Each binary variable $b_i$ is replaced by a qubit operator through the mapping $b_i \mapsto \frac{1- Z_i}{2}$



The terms can thus be transformed into pauli words via

$$
a_i b_i \mapsto \frac{a_i}{2}I - \frac{a_i}{2}Z_i
$$

$$
Q_{ij} b_i b_j \mapsto \frac{Q_{ij}}{4}(I - Z_i - Z_j + Z_i Z_j)
$$
### Minimal Prototype

- Simple Cubic Lattice
- Using the HP Model (Hydrophobic-polar protein folding model)[ \[ wiki \] ](https://en.wikipedia.org/wiki/Hydrophobic-polar_protein_folding_model)
- QUBOize to ready for VQE
- VQE

### Unminimal Prototype
- Replace HP Model with Miyazawa-Jernigan
- Turn encoding
- Better mixers

## References
- *Variational Quantum Eigensolver for Protein Folding using Neutral Atom Platforms*, Gefen Barnes, 2024
- *A perspective on protein structure prediction using quantum computers*, Doga et al, 2023
- *Resource-efficient quantum algorithm for protein folding*, Robert et al., 2021
- *Estimation of Effective Interresidue Contact Energies from Protein Crystal Structures: Quasi-Chemical Approximation*, Miyazawa & Jernigan, 1985
