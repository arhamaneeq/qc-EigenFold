import itertools, functools
import pandas as pd

from src.amino import is_hydrophobic
from src.qve import Estimate

if __name__ == "__main__":
    lattices = ["fcc", "bcc", "sc"]

    peptides = [
    # 16 dipeptides
    "AV", "RF", "DG", "HQ", "MN", "CT", "YW", "KE",
    "IL", "PS", "FA", "RE", "GN", "VK", "DY", "WL",

    # 16 tripeptides
    "AVF", "RDE", "HQL", "MNG", "CTY", "YWA", "KEP", "ILR",
    "PSV", "FAR", "REG", "GNK", "VLD", "DYH", "WLC", "TMA",

    # 3 tetrapeptides
    "AVFY", "RDEG", "MNGH"
]

    results = []

    for lattice in lattices:
        for peptide in peptides:
            try:
                res = Estimate(amino_seq=peptide, lattice=lattice)
                results.append(res)
            except Exception as e:
                print(f"Failed for {peptide} on {lattice}: {e}")
    
    df = pd.DataFrame(results)
    df.to_csv("data/eigenfold_summary.csv")

