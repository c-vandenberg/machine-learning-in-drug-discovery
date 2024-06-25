import logging
import numpy
from typing import Iterable, Union, List, Tuple
from numpy import ndarray
from pandas import DataFrame
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from mordred import Calculator, descriptors

logger = logging.getLogger(__name__)


def generate_canonical_smiles(smiles_dataset: Iterable) -> List:
    # Convert SMILES string into RDKit Mol object, which is a molecular graph representation
    mols = [Chem.MolFromSmiles(smiles_string) for smiles_string in smiles_dataset]

    # Convert each RDKit Mol object back into SMILES string, but in the canonical form
    canonical_smiles = [Chem.MolToSmiles(mol) for mol in mols]
    return canonical_smiles


def calculate_rdkit_descriptors(smiles_dataset: Iterable) -> Union[List, Tuple]:
    # Convert each SMILES string in the dataset into an RDKit Mol object
    mols = [Chem.MolFromSmiles(smiles_string) for smiles_string in smiles_dataset]

    # Descriptors.descList is a list of tuples where each tuple contains a descriptor name and its corresponding
    # function. Here we create a list of descriptor names by extracting the first element (name) from each tuple, and
    # initialise a descriptor calculator with the list of descriptor names
    descriptor_calc = MoleculeDescriptors.MolecularDescriptorCalculator(
        [rdkit_descriptor[0] for rdkit_descriptor in Descriptors.descList]
    )
    descriptor_names = descriptor_calc.GetDescriptorNames()

    mol_descriptors = []
    for mol in mols:
        try:
            # Some descriptors require explicit hydrogen atoms for accurate calculation
            mol = Chem.AddHs(mol)

            # Calculate all 200 RDKit descriptors for each molecule (e.g. MolWt, NumAtoms, NumBonds,
            # NumRotatableBonds etc.)
            mol_rdkit_descriptors = descriptor_calc.CalcDescriptors(mol)
            mol_descriptors.append(mol_rdkit_descriptors)
        except Exception as e:
            logger.exception(e)
    return mol_descriptors, descriptor_names


def calculate_morgan_fingerprints(smiles_dataset: Iterable) -> ndarray:
    morgan_fingerprints = []
    for smiles_string in smiles_dataset:
        # Convert SMILES string into an RDKit Mol object
        mol = Chem.MolFromSmiles(smiles_string)

        # Calculate Morgan Fingerprint for Mol object. Morgan Fingerprint radius is set to 2 (i.e. the Morgan algorithm
        # considers the atom itself and the atoms up to 2 bonds away/within 2 "hops"), and the size of the binary vector
        # is set to 2048 bits
        fingerprints = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, 2048)
        morgan_fingerprint = numpy.array(fingerprints)
        morgan_fingerprints.append(morgan_fingerprint)
    return numpy.array(morgan_fingerprints)


def calculate_mordred_descriptors(smiles_dataset: Iterable) -> DataFrame:
    """
    Mordred is a Python package that can calculate more than 1800 0D, 1D, 2D and 3D molecular descriptors

    """
    mordred_calculator = Calculator(descriptors, ignore_3D=False)
    mols = [Chem.MolFromSmiles(smiles_strings) for smiles_strings in smiles_dataset]

    # Calculate all >1800 Mordred descriptors for each Mol in a MordredDataFrame
    return mordred_calculator.pandas(mols)
