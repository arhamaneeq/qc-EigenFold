import math
import itertools

def build_lattice(nx: int, ny: int, nz: int, type: str):
    # will add a type param later, so that it can build FCC/BCC/etc.
    bases = {
         "sc": [(0.0,0.0,0.0)],
        "bcc": [(0.0,0.0,0.0), 
                (0.5,0.5,0.5)],
        "fcc": [(0.0,0.0,0.0),
                (0.0,0.5,0.5),
                (0.5,0.0,0.5),
                (0.5,0.5,0.0)],
        "hcp": [(0.0,0.0,0.0),
                (2/3,1/3,0.5)]
    }

    coords = []
    offsets = bases[type]

    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for ox, oy, oz in offsets:
                        coords.append((
                            x + ox, 
                            y + oy, 
                            z + oz
                        ))

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

    num_bits = math.ceil(math.log2(S))
    mapping = {
        r: [r * num_bits + k for k in range(num_bits)] 
        for r in range(R)
    }

    return mapping, num_bits

def choose_lattice_dims(R: int, ratio: float = 1.0):
    S_target = int(math.ceil(R * ratio))
    nx = ny = int(math.floor(math.sqrt(S_target)))
    nz = math.ceil(S_target / (nx * ny))
    return nx, ny, nz

def decode_positions(bits, mapping, S):
    positions = []
    for r, qubits in mapping.items():
        val = sum(bits[q] * (2 ** i) for i, q in enumerate(qubits))
        val %= S
        positions.append(val)
    return positions

def positions_to_coords(positions, coords):
    return [coords[p] for p in positions]

def decode_bitstring(bitstring, mapping, coords):
    bits = list(map(int, bitstring[::-1]))  # reverse for little-endian
    positions = []
    for r, qubits in mapping.items():
        val = sum(bits[q] * (2 ** i) for i, q in enumerate(qubits))
        val %= len(coords)
        positions.append(coords[val])
    return positions