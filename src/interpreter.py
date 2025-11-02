import numpy as np
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def classical_ground_state(linear, quadratic):
    max_linear = max((k if isinstance(k, int) else k[0]) for k in linear.keys())
    max_quad = max(max(pair) for pair in quadratic.keys()) if quadratic else -1
    n = max(max_linear, max_quad) + 1

    best_bits, best_E = None, float('inf')

    for bits in itertools.product([0, 1], repeat=n):
        E = 0.0
        # Linear terms
        for k, coeff in linear.items():
            idx = k if isinstance(k, int) else k[0]
            E += coeff * bits[idx]
        # Quadratic terms
        for (i, j), coeff in quadratic.items():
            E += coeff * bits[i] * bits[j]
        if E < best_E:
            best_E, best_bits = E, bits

    return np.array(best_bits), best_E

def visualize_folds(decoded_confs, seqH, name_prefix):
    n = len(decoded_confs)
    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols

    fig = plt.figure(figsize=(6*ncols, 6*nrows))

    for i, (bitstring, prob, decoded_coords) in enumerate(decoded_confs):
        ax = fig.add_subplot(nrows, ncols, i + 1, projection='3d')

        xs, ys, zs = zip(*decoded_coords)
        ax.plot(xs, ys, zs, '-o', color='black', lw=2, alpha=0.6)

        colors = ['orange' if h else 'cyan' for h in seqH]
        for (x, y, z), col in zip(decoded_coords, colors):
            ax.scatter(x, y, z, s=200, c=col, edgecolor='k')

        ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
        ax.set_title(f"{bitstring}\np={prob:.3f}")
        ax.set_box_aspect([1,1,1])

    plt.suptitle(f"Top {n} Stable Conformations for {name_prefix}", fontsize=14)
    plt.tight_layout()
    plt.savefig(f"images/{name_prefix}.png", dpi=300)
    # plt.show()

