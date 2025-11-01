import itertools

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
def add_linear(term_dict: dict, i, coeff):
    term_dict.setdefault((i,), 0.0)
    term_dict[(i,)] += coeff

def add_quadratic(term_dict: dict, i, j, coeff):
    key = tuple(sorted((i, j)))
    term_dict.setdefault(key, 0.0)
    term_dict[key] += coeff

def build_qubo(R, S, mapping, adjlist, seqH, A = 10.0, B = 10.0, C = 10.0):
    linear = {}
    quadrs = {}
    constant = 0.0

    for i in range(R):
        for j in range(i + 1, R):
            if seqH[i] and seqH[j]:
                Eij = -1.0
            else:
                Eij =  0.0
                continue

            for p in range(S):
                for q in adjlist[p]:
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
            add_quadratic(quadrs, i, j, 2.0 * B)

        constant += B * 1.0

    for r in range(R):
        rs = [mapping[(r, p)] for p in range(S)]
        rp = r + 1 if not r == R - 1 else 0

        for i in range(S):
            add_linear(linear, i, C * 1.0)
        for p in range(S):
            idx_rp = mapping[(r, p)]
            for q in adjlist[p]:
                idx_r1p = mapping[(rp, q)]
                add_quadratic(quadrs, idx_rp, idx_r1p, C * -1.0)

    return linear, quadrs, constant


