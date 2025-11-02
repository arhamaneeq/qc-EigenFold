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
    # so R is the residue number, i.e., the Rth (0-indexed) amino acid in the chain
    # row-maj flattening from (r,p) -> idx
    # grouped by residue
    
    mapping = {}

    for r in range(R):
        for p in range(S):
            mapping[(r, p)] = r * S + p

    return mapping

def choose_lattice_dims(R, ratio=1.3):
    # Target number of sites slightly above R
    S_target = int(math.ceil(R * ratio))
    side = math.ceil(S_target ** (1/3))
    nx = ny = nz = side
    return nx, ny, nz