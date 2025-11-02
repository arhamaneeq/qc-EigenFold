> [!NOTE]
> This repository was built as a submission for the *IISC Quantum Fall Fest 2025 | Cleveland Clinic QFF Hackathon Prompt*

> [!IMPORTANT]
> This README.md contains MathJax. Best viewed on GitHub on a Desktop Browser

## Setting up The Repository
After cloning, make sure you're using `Python 3.11.6`. 
On windows, set up the python environment using
```powershell
python -m venv .venv
./.venv/Scripts/Activate
pip install -r requirements.txt
```

You can then run the whole pipeline
```powershell
python main.py
python analysis.py # once main.py finishes running
```
    
# EigenFold


## QUBO

We use QUBO, i.e., for variables $b_i \in \{0, 1\}$ we have

$$
E(b_1, b_2, ..., b_N) = \sum_i a_i b_i + \sum_{i \lt j} b_i Q_{ij} b_j + c
$$


### One-Hot Encoding in QUBO

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


### Mapping QUBO Terms to Pauli Words

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

> [!NOTE]
> Algorithm ran successfully for `R = [True, False, True]` and `(nx, ny, nz) = (2, 2, 1)`, after running for 87.01s.

### Binary Encoding
One hot encoding is exopentially inefficient in the number of sites, $S$. For a 10 residue system in 3×3×3 lattice, $N = R \times S = 270$, which is completely infeasible for VQE, even on simulators. Thus it is completely imperative to switch to a more spatially efficient encoding scheme.

Instead of storing a seperate bit for every lattice site, we store the binary representation of the site index. Thus, for a lattice with $S$ sites, we have $n_{\text{bits per residue}} = \lceil {\log_2(S)} \rceil$. We can then assign each residue's position as:

$$
p_r = \sum_{k=0}^{n-1} {2^k x_{r,k}}
$$
$$
p_r = \sum_i w_i x_i
$$

where, $x_{r,k} \in \{0,1\}$.

Since $E = E_{backbone} + E_{collision} + E_{contact} + E_{const}$, we have

$$
E_{backbone}(r, r + 1) = C \times (|p_{r+1}-p_{r}| - 1)^2 
$$
which reduces to 
$$
C\times \big[(p_{r+1} - p{r})^2 - 2 p_{r + 1} + 2p_r + 1\big]
$$

$$
E_{collision}(i, j) = B \times (p_i - p_j)^2
$$
which reduces to 
$$
B \times \Big[\sum_i w_i^2 x_i + \sum_j w_j^2 x_j - 2 \sum_{x,j} w_i w_j x_i x_j \Big]
$$

$$
E_{contact}(i,j) = -A (1 - (|p_i - p_j| - 1)^2) 
$$

> [!NOTE]
> For `R = 3` and `S = 12`, we have `num_qubits = 12`, an improvement by a factor of 9, and takes 29.68s to execute.

The mapping to pauli words remains relatively unchanged.

## References
- *Variational Quantum Eigensolver for Protein Folding using Neutral Atom Platforms*, Gefen Barnes, 2024
- *A perspective on protein structure prediction using quantum computers*, Doga et al, 2023
- *Resource-efficient quantum algorithm for protein folding*, Robert et al., 2021
- *Estimation of Effective Interresidue Contact Energies from Protein Crystal Structures: Quasi-Chemical Approximation*, Miyazawa & Jernigan, 1985
