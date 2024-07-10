#!/usr/bin/env python3

import sys
import os
import time
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import PandasTools
from rdkit.Chem import rdMolDescriptors
from tqdm import tqdm
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.ensemble import RandomForestRegressor
import time
from sklearn.model_selection import ShuffleSplit, cross_validate,train_test_split
from lightgbm import LGBMRegressor

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from utils import avalon_fingerprints_utils


def main():
    # Load data
    homo_lumo_dataset: pd.DataFrame = pd.read_csv(
        os.path.join(
            os.path.dirname(__file__),
            '../../data/raw/orbital-energies-input-data.csv'
        )
    )

    # Add 2D structure columnd
    PandasTools.AddMoleculeColumnToFrame(
        homo_lumo_dataset,
        'SMILES',
        'Structure',
        includeFingerprints=True
    )

    # Generate Avalon fingerprints
    avalon_fpts: np.ndarray = avalon_fingerprints_utils.generate_avalon_fingerprints(
        homo_lumo_dataset['Structure']
    )

    # Insert into DataFrame. Each row represents a molecule's Avalon fingerprint and
    # each column an individual bit
    avalon_fpts_dataset: pd.DataFrame = pd.DataFrame(
        avalon_fpts,
        columns=['Bit {}'.format(bit) for bit in range(avalon_fpts.shape[1])]
    )

    avalon_fpts_dataset.head()


if __name__ == '__main__':
    main()
