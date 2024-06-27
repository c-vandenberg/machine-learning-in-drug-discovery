#!/usr/bin/env python3
import os
import pandas
from pandas import DataFrame
from rdkit import Chem
from utils import molecular_descriptors_utils
from constants.smiles_constants import SmilesConstants


def main():
    """
    1. Load SMILES data set
        * SMILES data set extracted from https://zivgitlab.uni-muenster.de/m_kueh11/fp-dm-tool
        * Data set gives SMILES strings and HOMO-LUMO energy gap (in meV).
    """
    mol_descriptors_data_dir: str = os.path.join(os.path.dirname(__file__), '../data')
    dataset: DataFrame = pandas.read_csv(os.path.join(mol_descriptors_data_dir, 'Orbital_Energies_input_data.csv'))

    """
    2. Generate Canonical SMILES and Remove Duplicate Canonical SMILES
        * SMILES strings are not unique identifiers. There are many equivalent options to generate a SMILES string for 
          a given structure.
        * We therefore have to canonicalise SMILES strings to generate a unique SMILES representation.
    """
    # Canonicalize SMILES and insert into dataframe
    dataset[SmilesConstants.SMILES] = molecular_descriptors_utils.generate_canonical_smiles(dataset.SMILES)

    # Remove any duplicate canonical SMILES from dataset
    sanitised_dataset: DataFrame = dataset.drop_duplicates(subset=[SmilesConstants.SMILES])

    """
    3. Calculate RDKit Descriptors for All Canonical SMILES in Data Set
    """
    # Calculate RDKit descriptors
    rdkit_mol_descriptors, desc_names = molecular_descriptors_utils.calculate_rdkit_descriptors(
        sanitised_dataset[SmilesConstants.SMILES]
    )
    rdk_descriptors_dataframe = pandas.DataFrame(rdkit_mol_descriptors)

    """
    4. Calculate Morgan Fingerprints for All Canonical SMILES in Data Set
    """
    morgan_fingerprints = molecular_descriptors_utils.calculate_morgan_fingerprints(
        sanitised_dataset[SmilesConstants.SMILES]
    )
    morgan_fingerprints_dataframe = pandas.DataFrame(
        morgan_fingerprints,
        columns=['Col_{}'.format(i) for i in range(morgan_fingerprints.shape[1])]
    )

    """
    5. Calculate Mordred Descriptors for All Canonical SMILES in Data Set
    """
    mordred_descriptors = molecular_descriptors_utils.calculate_mordred_descriptors(
        sanitised_dataset[SmilesConstants.SMILES]
    )


if __name__ == '__main__':
    main()
