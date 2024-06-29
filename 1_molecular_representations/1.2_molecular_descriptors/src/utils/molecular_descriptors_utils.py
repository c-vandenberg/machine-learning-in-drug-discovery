import logging
import numpy
from typing import Union, List, Tuple
from numpy import ndarray
from pandas import Series
from pandas import DataFrame
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Mol
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
from mordred import Calculator, descriptors

logger = logging.getLogger(__name__)


def generate_canonical_smiles(smiles_dataset: Series) -> List:
    # Convert SMILES string into RDKit Mol object, which is a molecular graph representation
    mols: List = [Chem.MolFromSmiles(smiles_string) for smiles_string in smiles_dataset]

    # Convert each RDKit Mol object back into SMILES string, but in the canonical form
    canonical_smiles: List = [Chem.MolToSmiles(mol) for mol in mols]

    return canonical_smiles


def calculate_maccs_keys_fingerprints(smiles_dataset: Series) -> ndarray:
    maccs_keys: List = []
    # Convert each SMILES string in the dataset into an RDKit Mol object
    for smiles_string in smiles_dataset:
        mol: Mol = Chem.MolFromSmiles(smiles_string)

        # Calculate MACCS-keys fingerprint for Mol object
        bit_vect_maccs_keys_fingerprint: ExplicitBitVect = MACCSkeys.GenMACCSKeys(mol)
        maccs_keys_fingerprint: ndarray = numpy.array(bit_vect_maccs_keys_fingerprint)
        maccs_keys.append(maccs_keys_fingerprint)

    return numpy.array(maccs_keys)


def calculate_avalon_fingerprints(smiles_dataset: Series) -> ndarray:
    avalon_fingerprints: List = []
    for smiles_string in smiles_dataset:
        # Convert SMILES string into an RDKit Mol object
        mol: Mol = Chem.MolFromSmiles(smiles_string)

        # Calculate Avalon fingerprint for Mol object
        bit_vect_avalon_fingerprint: ExplicitBitVect = pyAvalonTools.GetAvalonFP(mol, nBits=512)
        avalon_fingerprint: ndarray = numpy.array(bit_vect_avalon_fingerprint)
        avalon_fingerprints.append(avalon_fingerprint)

    return numpy.array(avalon_fingerprints)


def calculate_atom_pairs_fingerprints(smiles_dataset: Series) -> ndarray:
    atom_pairs_fingerprints: List = []
    for smiles_string in smiles_dataset:
        # Convert SMILES string into an RDKit Mol object
        mol: Mol = Chem.MolFromSmiles(smiles_string)

        # Calculate atoms-pairs fingerprint for Mol object
        bit_vect_atom_pairs_fingerprint: ExplicitBitVect = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
            mol,
            nBits=512
        )
        atom_pairs_fingerprint: ndarray = numpy.array(bit_vect_atom_pairs_fingerprint)
        atom_pairs_fingerprints.append(atom_pairs_fingerprint)

    return numpy.array(atom_pairs_fingerprints)


def calculate_tt_fingerprints(smiles_dataset: Series) -> ndarray:
    tt_fingerprints: List = []
    for smiles_string in smiles_dataset:
        # Convert SMILES string into an RDKit Mol object
        mol: Mol = Chem.MolFromSmiles(smiles_string)

        # Calculate topological-torsions fingerprint for Mol object
        bit_vect_tt_fingerprint: ExplicitBitVect = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
            mol,
            nBits=512
        )
        tt_fingerprint: ndarray = numpy.array(bit_vect_tt_fingerprint)
        tt_fingerprints.append(tt_fingerprint)

    return numpy.array(tt_fingerprints)


def calculate_morgan_fingerprints(smiles_dataset: Series) -> ndarray:
    morgan_fingerprints: List = []
    for smiles_string in smiles_dataset:
        # Convert SMILES string into an RDKit Mol object
        mol: Mol = Chem.MolFromSmiles(smiles_string)

        # Calculate Morgan Fingerprint for Mol object. Morgan Fingerprint radius is set to 2 (i.e. the Morgan algorithm
        # considers the atom itself and the atoms up to 2 bonds away/within 2 "hops"), and the size of the binary vector
        # is set to 2048 bits
        bit_vect_morg_fingerprint: ExplicitBitVect = rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol,
            2,
            2048
        )
        morgan_fingerprint: ndarray = numpy.array(bit_vect_morg_fingerprint)
        morgan_fingerprints.append(morgan_fingerprint)

    return numpy.array(morgan_fingerprints)


def calculate_rdkit_descriptors(smiles_dataset: Series) -> Union[List, Tuple]:
    # Convert each SMILES string in the dataset into an RDKit Mol object
    mols: List = [Chem.MolFromSmiles(smiles_string) for smiles_string in smiles_dataset]

    # Descriptors.descList is a list of tuples where each tuple contains a descriptor name and its corresponding
    # function. Here we create a list of descriptor names by extracting the first element (name) from each tuple, and
    # initialise a descriptor calculator with the list of descriptor names
    descriptor_calc: MolecularDescriptorCalculator = MolecularDescriptorCalculator(
        [rdkit_descriptor[0] for rdkit_descriptor in Descriptors.descList]
    )
    descriptor_names: tuple = descriptor_calc.GetDescriptorNames()

    mol_descriptors: List = []
    for mol in mols:
        try:
            # Some descriptors require explicit hydrogen atoms for accurate calculation
            mol: Mol = Chem.AddHs(mol)

            # Calculate all 200 RDKit descriptors for each molecule (e.g. MolWt, NumAtoms, NumBonds,
            # NumRotatableBonds etc.)
            mol_rdkit_descriptors: tuple = descriptor_calc.CalcDescriptors(mol)
            mol_descriptors.append(mol_rdkit_descriptors)
        except Exception as e:
            logger.exception(e)

    return mol_descriptors, descriptor_names


def calculate_mordred_descriptors(smiles_dataset: Series) -> DataFrame:
    """
    Mordred is a Python package that can calculate more than 1800 0D, 1D, 2D and 3D molecular descriptors

    """
    mordred_calculator: Calculator = Calculator(descriptors, ignore_3D=False)
    mols: List = [Chem.MolFromSmiles(smiles_strings) for smiles_strings in smiles_dataset]

    # Calculate all >1800 Mordred descriptors for each Mol in a MordredDataFrame
    return mordred_calculator.pandas(mols)
