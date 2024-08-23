import logging
import numpy
import pandas
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
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.ML.Cluster import Butina
from mordred import Calculator, descriptors

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)


def generate_canonical_smiles(smiles_dataset: Series) -> List:
    """
    Convert SMILES strings to their canonical form.

    Parameters
    ----------
    smiles_dataset : Series
        A series of SMILES strings.

    Returns
    -------
    List[str]
        A list of canonical SMILES strings.
    """
    # Convert SMILES string into RDKit Mol object, which is a molecular graph representation
    mols: List = [Chem.MolFromSmiles(smiles_string) for smiles_string in smiles_dataset]

    # Convert each RDKit Mol object back into SMILES string, but in the canonical form
    canonical_smiles: List = [Chem.MolToSmiles(mol) for mol in mols]

    return canonical_smiles


def calculate_maccs_keys_fingerprints(smiles_dataset: Series) -> ndarray:
    """
    Calculate MACCS-keys fingerprints for a dataset of SMILES strings.

    Parameters
    ----------
    smiles_dataset : Series
        A series of SMILES strings.

    Returns
    -------
    ndarray
        A numpy array of MACCS-keys fingerprints.
    """
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
    """
    Calculate Avalon fingerprints for a dataset of SMILES strings.

    Parameters
    ----------
    smiles_dataset : Series
        A series of SMILES strings.

    Returns
    -------
    ndarray
        A numpy array of Avalon fingerprints.
    """
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
    """
    Calculate Atom-Pairs fingerprints for a dataset of SMILES strings.

    Parameters
    ----------
    smiles_dataset : Series
        A series of SMILES strings.

    Returns
    -------
    ndarray
        A numpy array of Atom Pairs fingerprints.
    """
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
    """
    Calculate Topological-Torsions fingerprints for a dataset of SMILES strings.

    Parameters
    ----------
    smiles_dataset : Series
        A series of SMILES strings.

    Returns
    -------
    ndarray
        A numpy array of Topological-Torsions fingerprints.
    """
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
    """
    Calculate Morgan fingerprints for a dataset of SMILES strings.

    Parameters
    ----------
    smiles_dataset : Series
        A series of SMILES strings.

    Returns
    -------
    ndarray
        A numpy array of Morgan fingerprints.
    """
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
    """
    Calculate RDKit descriptors for a dataset of SMILES strings.

    Parameters
    ----------
    smiles_dataset : Series
        A series of SMILES strings.

    Returns
    -------
    Union[List, Tuple]
        A list of RDKit descriptors and a tuple of descriptor names.
    """
    # Convert each SMILES string in the dataset into an RDKit Mol object if SMILES string is valid
    mols: List = []
    for smiles_string in smiles_dataset:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            logger.warning(f"Invalid SMILES string: {smiles_string}")
        else:
            mols.append(mol)

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
            Chem.SanitizeMol(mol)
        except Exception as e:
            logger.exception(f"Sanitization failed for molecule: {Chem.MolToSmiles(mol)} with error: {e}")
            continue

        try:
            # Some descriptors require explicit hydrogen atoms for accurate calculation
            mol: Mol = Chem.AddHs(mol)

            # Calculate all 200 RDKit descriptors for each molecule (e.g. MolWt, NumAtoms, NumBonds,
            # NumRotatableBonds etc.)
            mol_rdkit_descriptors: tuple = descriptor_calc.CalcDescriptors(mol)
            mol_descriptors.append(mol_rdkit_descriptors)
        except Exception as e:
            logger.exception(f"Descriptor calculation failed for molecule: {Chem.MolToSmiles(mol)} with error: {e}")

    return mol_descriptors, descriptor_names


def calculate_mordred_descriptors(smiles_dataset: Series) -> DataFrame:
    """
    Calculate Mordred descriptors for a dataset of SMILES strings.

    Mordred is a Python package that can calculate more than 1800 0D, 1D, 2D and 3D molecular descriptors

    Parameters
    ----------
    smiles_dataset : Series
        A series of SMILES strings.

    Returns
    -------
    DataFrame
        A dataframe of Mordred descriptors.
    """
    mordred_calculator: Calculator = Calculator(descriptors, ignore_3D=False)
    mols: List = [Chem.MolFromSmiles(smiles_strings) for smiles_strings in smiles_dataset]

    # Calculate all >1800 Mordred descriptors for each Mol in a MordredDataFrame
    return mordred_calculator.pandas(mols)


def tanimoto_similarity_smiles_search(query_compound_smiles: str, smiles_dataset: Series,
                                      similarity_threshold=0.5) -> DataFrame:
    """
    Perform a Tanimoto similarity search on a dataset of SMILES strings.

    Parameters
    ----------
    query_compound_smiles : str
        The SMILES string of the query compound.
    smiles_dataset : Series
        A series of SMILES strings.
    similarity_threshold : float (Optional)
        The similarity threshold for filtering results (default = 0.5).

    Returns
    -------
    DataFrame
        A dataframe of Tanimoto similarities that meet the threshold.
    """
    query_compound_mol: Mol = AllChem.MolFromSmiles(query_compound_smiles)

    query_compound_morgan_fingerprint: ExplicitBitVect = AllChem.GetMorganFingerprintAsBitVect(
        query_compound_mol,
        2,
        nBits=4096
    )

    # Convert each SMILES string in the dataset into an RDKit Mol object
    dataset_mols: List = [Chem.MolFromSmiles(smiles_string) for smiles_string in smiles_dataset]

    # Calulate dataset Morgan fingerprints
    dataset_morgan_fingerprints: List = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=4096)
                                         for mol in dataset_mols]

    # Calculate Tanimoto coefficient of query compound against each compound in the data set and insert into list
    tanimoto_coeffs: List = [
        DataStructs.FingerprintSimilarity(
            query_compound_morgan_fingerprint,
            dataset_morgan_fingerprint,
            metric=DataStructs.TanimotoSimilarity
        ) for dataset_morgan_fingerprint in dataset_morgan_fingerprints
    ]

    tanimoto_coeffs_dataframe: DataFrame = pandas.DataFrame(tanimoto_coeffs, columns=['Tanimoto Coefficient'])

    # Sort Tanimoto coefficient values in descending order and only return those compounds whose values are less than
    # or equal to the predefined similarity threshold (default value of 0.5)
    sorted_tanimoto_coeffs_dataframe: DataFrame = tanimoto_coeffs_dataframe.sort_values(
        'Tanimoto Coefficient',
        ascending=False)

    return sorted_tanimoto_coeffs_dataframe[
        sorted_tanimoto_coeffs_dataframe['Tanimoto Coefficient'] >= similarity_threshold
    ]


def taylor_butina_tanimoto_clustering(smiles_dataset: Series) -> List:
    """
    Perform a Taylor-Butina cluster analysis on a dataset of SMILES strings based on Tanimoto similarity.

    Parameters
    ----------
    smiles_dataset : Series
        A series of SMILES strings.

    Returns
    -------
    List[int]
        A list of cluster IDs for each molecule.
    """
    # Convert each SMILES string in the dataset into an RDKit Mol object
    dataset_mols: List = [Chem.MolFromSmiles(smiles_string) for smiles_string in smiles_dataset]

    # Calulate dataset Morgan fingerprints
    dataset_morgan_fingerprints: List = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=4096)
                                         for mol in dataset_mols]
    num_morgan_fingerprints: int = len(dataset_morgan_fingerprints)

    dist_matrix: List = []
    for index in range(1, num_morgan_fingerprints):
        # Calculate Tanimoto coefficients between current index fingerprint and all preceding fingerprints in data set
        tanimoto_coefficients = DataStructs.BulkTanimotoSimilarity(
            dataset_morgan_fingerprints[index],
            dataset_morgan_fingerprints[:index]
        )

        # Convert Tanimoto coeffieicnts to distances (1 - Tanimoto coefficients) and extend dist_matrix
        dist_matrix.extend([1 - tanimoto_coeff for tanimoto_coeff in tanimoto_coefficients])

    # Cluster the Tanimoto similarity scores/coefficients using Taylot-Butina algorithm
    taylor_butina_clusters = Butina.ClusterData(dist_matrix, num_morgan_fingerprints, 0.3, isDistData=True)
    cluster_id_list = [0] * len(smiles_dataset)

    for index, cluster in enumerate(taylor_butina_clusters, 1):
        for member in cluster:
            cluster_id_list[member] = index

    return cluster_id_list
