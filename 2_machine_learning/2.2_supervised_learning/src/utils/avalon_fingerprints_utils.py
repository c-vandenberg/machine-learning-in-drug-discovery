import numpy as np
from numpy import ndarray
from pandas import Series
from tqdm import tqdm
from rdkit import Chem
from rdkit.Avalon import pyAvalonTools
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from typing import List, Union


def generate_avalon_fingerprints_from_mol(mol_structure_data: Series, nbits: int) -> ndarray:
    return _avalon_fingerprints_from_mol(mol_structure_data, nbits)


def generate_avalon_fingerprints_from_smiles(smiles_data: Series, nbits: int) -> ndarray:
    mols: List = [Chem.MolFromSmiles(smiles) for smiles in smiles_data if smiles is not None]

    return _avalon_fingerprints_from_mol(mols, nbits)


def _avalon_fingerprints_from_mol(mol_structure_data: Union[Series, List], nbits: int) -> ndarray:
    avalon_fingerprints: List = []

    for mol in tqdm(mol_structure_data):

        # Each Avalon molecular fingerprint is a binary vector of `nbits` bits (`nbits` binary positions)
        avalon_fingerprint: ExplicitBitVect = pyAvalonTools.GetAvalonFP(mol, nBits=nbits)
        avalon_fingerprints.append(avalon_fingerprint)

    return np.array(avalon_fingerprints)
