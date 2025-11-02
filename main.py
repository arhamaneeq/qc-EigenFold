import itertools, functools

from src.amino import is_hydrophobic
from src.qve import Estimate

if __name__ == "__main__":
    while True:
        poly = input("ENTER POLYPEPTIDE SEQUENCE\n")
        if len(poly) <= 5:
            break
        else: 
            continue

    Estimate(amino_seq=poly)