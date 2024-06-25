#!/usr/bin/env python3
import os
import pandas
from pandas import DataFrame
from rdkit import Chem
from utils import molecular_descriptors_utils
from constants.smiles_constants import SmilesConstants


def main():
    mol_descriptors_data_dir: str = os.path.join(os.path.dirname(__file__), '../data')

    dataset: DataFrame = pandas.read_csv(os.path.join(mol_descriptors_data_dir, 'Orbital_Energies_input_data.csv'))

    # Canonicalize SMILES and insert into dataframe
    dataset[SmilesConstants.SMILES] = molecular_descriptors_utils.generate_canonical_smiles(dataset.SMILES)

    # Remove any duplicate canonical SMILES from dataset
    sanitised_dataset: DataFrame = dataset.drop_duplicates(subset=[SmilesConstants.SMILES])

    # Calculate RDKit descriptors
    rdkit_mol_descriptors, desc_names = molecular_descriptors_utils.calculate_rdkit_descriptors(
        sanitised_dataset[SmilesConstants.SMILES]
    )
    rdk_descriptors_dataframe = pandas.DataFrame(rdkit_mol_descriptors)

    # Calculate Morgan Fingerprints
    morgan_fingerprints = molecular_descriptors_utils.calculate_morgan_fingerprints(
        sanitised_dataset[SmilesConstants.SMILES]
    )
    morgan_fingerprints_dataframe = pandas.DataFrame(
        morgan_fingerprints,
        columns=['Col_{}'.format(i) for i in range(morgan_fingerprints.shape[1])]
    )

    # Calculate Mordred descriptors
    mordred_descriptors = molecular_descriptors_utils.calculate_mordred_descriptors(
        sanitised_dataset[SmilesConstants.SMILES]
    )


if __name__ == '__main__':
    main()
