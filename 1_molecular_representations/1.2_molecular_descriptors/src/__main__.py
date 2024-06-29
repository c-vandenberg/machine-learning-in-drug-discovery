#!/usr/bin/env python3
import os
import pandas
from pandas import DataFrame
from numpy import ndarray
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
    dataset: DataFrame = pandas.read_csv(os.path.join(mol_descriptors_data_dir, 'orbital_energies_input_data.csv'))

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
    3. Molecular Descriptors in the RDKit Library
       In this section, we will generate the following molecular descriptors via the RDKit library:
        * Molecular ACCess System keys or MACCS-keys
        * Avalon fingerprint
        * Atom-pair fingerprint
        * Topological-Torsions fingerprint.
        * Morgan fingerprint (Circular Fingerprint)
        * RDKit Fingerprint
    """

    """
    3.1  Calculate MACCS-Keys Fingerprint (166 Public Keys) for All Canonical SMILES in Data Set
         The MACCS-keys are a set of structural keys encoding for a set predefined substructures/fragments, with each 
         bit indicating the absence or presence of a particular substructure/fragment.
    """
    maccs_keys_fingerprints: ndarray = molecular_descriptors_utils.calculate_maccs_keys_fingerprints(
        sanitised_dataset[SmilesConstants.SMILES]
    )

    """
    3.2 Calculate Avalon Fingerprints for All Canonical SMILES in Data Set
        Similar to Daylight fingerprints, Avalon uses a fingerprint generator that enumerates certain paths and feature 
        classes of the molecular graph.
    """
    avalon_fingerprints: ndarray = molecular_descriptors_utils.calculate_avalon_fingerprints(
        sanitised_dataset[SmilesConstants.SMILES]
    )

    """
    3.3 Calculate Atom-Pairs Fingerprints for All Canonical SMILES in Data Set
        Atom-pairs is a simple type of substructure that is defined in terms of the atomic environments of, and 
        shortest path separation between all pairs of atoms in the topological representation of a chemical structure.

        The atom-pairs fingerprint is an algorithm for computing atom pairs from this aforementioned representation.

    """
    atom_pairs_fingerprints: ndarray = molecular_descriptors_utils.calculate_atom_pairs_fingerprints(
        sanitised_dataset[SmilesConstants.SMILES]
    )

    """
    3.4 Calculate Topological-Torsions (TT) Fingerprints for All Canonical SMILES in Data Set
        The topological-torsion (TT) fingerprint is a molecular descriptor that consists of four consecutively bonded 
        non-hydrogen atoms along with the number of non-hydrogen branches. The descriptor is basically the topological 
        analogue of the basic conformational element, the torsion angle.

        The atom-pair and TT fingerprints capture and magnify distinct aspects of molecular topology.
    """
    tt_fingerprints: ndarray = molecular_descriptors_utils.calculate_tt_fingerprints(
        sanitised_dataset[SmilesConstants.SMILES]
    )

    """
    3.5 Calculate Morgan Fingerprints for All Canonical SMILES in Data Set
    """
    morgan_fingerprints: ndarray = molecular_descriptors_utils.calculate_morgan_fingerprints(
        sanitised_dataset[SmilesConstants.SMILES]
    )
    morgan_fingerprints_dataframe: DataFrame = pandas.DataFrame(
        morgan_fingerprints,
        columns=['Col_{}'.format(i) for i in range(morgan_fingerprints.shape[1])]
    )

    """
    3.6 Calculate RDKit Descriptors for All Canonical SMILES in Data Set
    """
    # Calculate RDKit descriptors
    rdkit_mol_descriptors, desc_names = molecular_descriptors_utils.calculate_rdkit_descriptors(
        sanitised_dataset[SmilesConstants.SMILES]
    )
    rdk_descriptors_dataframe: DataFrame = pandas.DataFrame(rdkit_mol_descriptors)

    """
    3.7 Calculate Mordred Descriptors for All Canonical SMILES in Data Set
       Mordred is a Python package that can calculate more than 1800 0D, 1D, 2D and 3D molecular descriptors
    """
    mordred_descriptors: DataFrame = molecular_descriptors_utils.calculate_mordred_descriptors(
        sanitised_dataset[SmilesConstants.SMILES]
    )


if __name__ == '__main__':
    main()
