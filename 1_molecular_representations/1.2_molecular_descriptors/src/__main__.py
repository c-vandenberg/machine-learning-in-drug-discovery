#!/usr/bin/env python3
import os
import pandas
from pandas import DataFrame
from numpy import ndarray
from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem import DataStructs
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
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
    4. Calculate Molecular Descriptors (RDKit & Mordred)
    """

    """
    4.1 Calculate RDKit Descriptors for All Canonical SMILES in Data Set
    """


    """
    4.2 Calculate Mordred Descriptors for All Canonical SMILES in Data Set
        Mordred is a Python package that can calculate more than 1800 0D, 1D, 2D and 3D molecular descriptors
    """


    """
    5. Calculate Tanimoto Coefficient
        The Tanimoto similarity algorithm provides a measure of similarity between the molecular fingerprints of two 
        molecules. Usually, these two molecular fingerprints are represented as two sets of fingerprint 'bits', denoted 
        as A and B.

        The Tanimoto coefficient, T(A,B), is calculated as the ratio of the intersection of A and B to the union of A 
        and B (Fig 2). This coefficient ranges from 0, indicating no common bits between the fingerprints, to 1, 
        representing identical fingerprints.
        
        In this section we will calculate the Tanimoto coefficient for the Morgan fingerprints of the first two
         molecules in the data set
    """
    # Calcualte Morgan fingerprints
    molecule_2: Mol = AllChem.MolFromSmiles(sanitised_dataset['SMILES'][2])
    molecule_10: Mol = AllChem.MolFromSmiles(sanitised_dataset['SMILES'][10])

    molecule_2_morgan_fingerprint: ExplicitBitVect = AllChem.GetMorganFingerprintAsBitVect(
        molecule_2,
        2,
        nBits=2048,
        bitInfo={}
    )
    molecule_10_morgan_fingerprint: ExplicitBitVect = AllChem.GetMorganFingerprintAsBitVect(
        molecule_10,
        2,
        nBits=2048,
        bitInfo={}
    )

    # Calculate Tanimoto coefficient betwee Morgan fingerprints using RDKit
    molecule_2_10_tanimoto_coeff = DataStructs.FingerprintSimilarity(
        molecule_2_morgan_fingerprint,
        molecule_10_morgan_fingerprint,
        metric=DataStructs.TanimotoSimilarity
    )

    """
    6. Applications of Molecular Fingerprints
        The two applications of molecular fingerprints we will look at in this section are:
            1. Searching for Compounds (e.g. using a Tanimoto similarity search)
            2. Compound Clustering
    """

    """
    6.1 Searching for Compounds Similar to a Query Compound
        In this section we will query the sanitised dataset for all compounds similar to the compound acetylsalicylic 
        acid (aspirin).

        This will be achieved by performing a Tanimoto similarity search of the Morgan fingerprint for acetylsalicylic 
        acid against the Morgan fingerprints of all compounds in the sanitised dataset
    """
    tanimoto_coeffs = molecular_descriptors_utils.tanimoto_similarity_smiles_search(
        'O=C(C)Oc1ccccc1C(=O)O',
        sanitised_dataset[SmilesConstants.SMILES]
    )

    """
    6.2 Compound Clustering (Taylor-Butina Clustering)
        Cluster Analysis or Clustering is the task of grouping of a set of objects (usually unlabelled data) in such a 
        way that objects in the same group (a cluster) are more similar to each other than those in other 
        groups/clusters.

        As clustering/cluster analysis is used on unlabelled data (i.e. data where the elements have no distinct 
        identifiers or classifications), it is defined under the branch of unsupervised learning in machine learning. 
        Unsupervised learning aims at gaining insights from unlabelled data points.

        The algorithm we will use to perform our clustering/cluster analysis is called Taylor-Butina Clustering. 6 7

        Taylor-Butina clustering is a unsupervised, non-hierarchical clustering method that guarantees every cluster 
        contains molecules which are within a distance cutoff from to the central molecule. 
            * Non-hierarchical clustering here refers to clustering methods that do not produce a nested series of 
            clusters. Instead they create a single paritioning of the data into a set of distinct clusters. 
            * This contrasts with hierarchical clustering, which creates a tree-like structure of clusters that can be 
            divided or combined at various level
    """
    taylor_butina_clusters = molecular_descriptors_utils.taylor_butina_tanimoto_clustering(
        sanitised_dataset[SmilesConstants.SMILES]
    )

    test = 'test'

if __name__ == '__main__':
    main()
