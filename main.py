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