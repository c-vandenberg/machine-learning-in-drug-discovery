import numpy as np
from numpy import ndarray
from pandas import Series
from tqdm import tqdm
from rdkit.Avalon import pyAvalonTools
from typing import List


def generate_avalon_fingerprints(data: Series) -> ndarray:
    avalon_fingerprints: List = []

    for mol in tqdm(data):

        # Each Avalon molecular fingerprint is a binary vector of 4096 bits/4096 binary positions
        avalon_fingerprint = pyAvalonTools.GetAvalonFP(mol, nBits=4096)
        avalon_fingerprints.append(avalon_fingerprint)

    return np.array(avalon_fingerprints)
