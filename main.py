import itertools, functools
import pandas as pd

from src.amino import is_hydrophobic
from src.qve import Estimate

if __name__ == "__main__":
    lattices = ["fcc", "bcc", "sc"]
    peptides = ["ABC", "XYZ", "VAR", "ARH", "EVA", "RON", "FEY"]

    results = []

    for lattice in lattices:
        for peptide in peptides:
            try:
                res = Estimate(amino_seq=peptide, lattice=lattice)
                results.append(res)
            except Exception as e:
                print(f"Failed for {peptide} on {lattice}: {e}")
    
    df = pd.DataFrame(results)
    df.to_csv("eigenfold_summary.csv")

