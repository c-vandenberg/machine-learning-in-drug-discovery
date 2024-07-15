#!/usr/bin/env python3

import sys
import os
import time
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit.Chem import PandasTools
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import ShuffleSplit, cross_validate, train_test_split
from lightgbm import LGBMRegressor
from typing import List, Dict

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from utils import avalon_fingerprints_utils


def main():
    """
    1. Reactivity Prediction: HOMO-LUMO Energy Gap of Organic Molecules
    """
    """
    1.1 Data Preparation
    """
    # Load HOMO-LUMO energy dataset
    homo_lumo_dataset: pd.DataFrame = pd.read_csv(
        os.path.join(
            os.path.dirname(__file__),
            '../../data/raw/orbital-energies-input-data.csv'
        )
    )

    # Add 2D structure column
    PandasTools.AddMoleculeColumnToFrame(
        homo_lumo_dataset,
        'SMILES',
        'Structure',
        includeFingerprints=True
    )

    # Generate Avalon fingerprints from mol structure
    avalon_fpts: np.ndarray = avalon_fingerprints_utils.generate_avalon_fingerprints_from_mol(
        homo_lumo_dataset['Structure'],
        4096
    )

    # Insert into DataFrame. Each row represents a molecule's Avalon fingerprint and
    # each column an individual bit
    homo_lumo_ava_fpts_dataset: pd.DataFrame = pd.DataFrame(
        avalon_fpts,
        columns=['Bit {}'.format(bit) for bit in range(avalon_fpts.shape[1])]
    )

    """
    1.2 Instantiate LightGBM Regressor and Random Forest Regressor Models
    """
    # Instantiate the `LGBMRegressor` model
    lgbm_regressor: LGBMRegressor = LGBMRegressor(n_estimators=800, random_state=42)

    # Instantiate the `RandomForestRegressor`model
    rf_regressor: RandomForestRegressor = RandomForestRegressor(random_state=42)

    """
    1.3 10-Fold Cross-Validation of LightGBM Regressor Model using Avalon Fingerprints & HOMO-LUMO Energy Gap Dataset
    """
    # Start timer to calculate total execution time
    start_time: float = time.time()

    # Create a cross-validator that randomly shuffles and splits data into training and test sets
    homo_lumo_cross_validator: ShuffleSplit = ShuffleSplit(n_splits=10, test_size=0.3, random_state=42)

    # Define scoring metrics.
    homo_lumo_scoring_metrics: List[str] = ['r2', 'neg_mean_absolute_error']

    # Cross-validation execution.
    homo_lumo_cross_validation_scores: Dict = cross_validate(
        lgbm_regressor,
        homo_lumo_ava_fpts_dataset,
        homo_lumo_dataset.Energygap,
        scoring=homo_lumo_scoring_metrics,
        cv=homo_lumo_cross_validator
    )

    print(homo_lumo_cross_validation_scores)

    end_time: float = time.time()
    execution_time: float = end_time - start_time
    print('\n HOMO-LUMO Energy Gap Model Cross-Validation Execution Time: ', round(execution_time / 60, 2), 'min')

    # Output mean values of R² and MAE for 10-fold cross-validation using LightGBM
    print(
        'HOMO-LUMO Energy Gap Model Cross-Validation R²: ',
        round(np.mean(homo_lumo_cross_validation_scores['test_r2']), 2),
        '\n'
    )
    print(
        'HOMO-LUMO Energy Gap Model Cross-Validation Mean Absolute Error (MAE): ',
        round(np.mean(-homo_lumo_cross_validation_scores['test_neg_mean_absolute_error']), 2)
    )

    """
    1.4 Split Avalon Fingerprints and HOMO-LUMO Energy Gap Dataset Into Training Data & Testing Data to Train & Test 
        LightGBM Regressor Model
    """
    # Split data into training and test sets using `train_test_split()` function
    (homo_lumo_model_x_train_data, homo_lumo_model_x_test_data,
     homo_lumo_model_y_train_data, homo_lumo_model_y_test_data) = train_test_split(
        homo_lumo_ava_fpts_dataset,
        homo_lumo_dataset.Energygap,
        test_size=0.3,
        random_state=42
    )

    # LightGBM model training
    ava_homo_lumo_lgbm_model: LGBMRegressor = lgbm_regressor.fit(
        homo_lumo_model_x_train_data,
        homo_lumo_model_y_train_data
    )

    # Test trained LightGBM model using testing feature matrix (Avalon fingerprints)
    # Outputs a predicted target vector (predicted HOMO-LUMO energy gap)
    homo_lumo_lgbm_model_predict: np.array = ava_homo_lumo_lgbm_model.predict(
        homo_lumo_model_x_test_data
    )

    # Measure model performance with R score between actual and predicted HOMO-LUMO gap values
    # Slice correlation between predicted labels and actual labels from correlation matrix
    homo_lumo_model_predict_r: np.float64 = np.corrcoef(
        homo_lumo_model_y_test_data,
        homo_lumo_lgbm_model_predict
    )[0, 1]

    print(
        '\n LightGBM Regressor HOMO-LUMO Energy Gap Prediction R: ',
        round(homo_lumo_model_predict_r, 2)
    )

    # Measure model performance with R² score between actual and predicted HOMO-LUMO gap values
    homo_lumo_model_predict_r2: float = r2_score(
        homo_lumo_model_y_test_data,
        homo_lumo_lgbm_model_predict
    )

    print(
        '\n LightGBM Regressor HOMO-LUMO Energy Gap Prediction R²: ',
        round(homo_lumo_model_predict_r2, 2)
    )

    # Measure model performance by calculating MAE between actual and predicted HOMO-LUMO gap values
    homo_lumo_model_predict_mae: np.float64 = mean_absolute_error(
        homo_lumo_model_y_test_data,
        homo_lumo_lgbm_model_predict
    )

    print(
        '\n LightGBM Regressor HOMO-LUMO Energy Gap Prediction MAE: ',
        round(homo_lumo_model_predict_mae, 2)
    )

    """
    2. Reactivity Prediction: Buchwald-Hartwig Amination (C-N Cross-Coupling) Reaction Yield
    """
    """
    2.1 Data Preparation
    """
    # Load data
    buchwald_amin_yield_dataset: pd.DataFrame = pd.read_csv(
        os.path.join(
            os.path.dirname(__file__),
            '../../data/raw/buchwald_yield_data.csv'
        )
    )

    # Calculate Avalon fingerprints for all four cross-coupling components
    ligand_ava_fpts: ExplicitBitVect = avalon_fingerprints_utils.generate_avalon_fingerprints_from_smiles(
        buchwald_amin_yield_dataset['Ligand'],
        2048
    )
    additive_ava_fpts: ExplicitBitVect = avalon_fingerprints_utils.generate_avalon_fingerprints_from_smiles(
        buchwald_amin_yield_dataset['Additive'],
        1024
    )
    base_ava_fpts: ExplicitBitVect = avalon_fingerprints_utils.generate_avalon_fingerprints_from_smiles(
        buchwald_amin_yield_dataset['Base'],
        1024
    )
    aryl_halide_ava_fpts: ExplicitBitVect = avalon_fingerprints_utils.generate_avalon_fingerprints_from_smiles(
        buchwald_amin_yield_dataset['Aryl halide'],
        1024
    )

    # Concatenate Avalon fingerprints of the four cross-coupling components
    buchwald_amin_ava_fpts = np.concatenate(
        [ligand_ava_fpts, additive_ava_fpts, base_ava_fpts, aryl_halide_ava_fpts],
        axis=1
    )

    # Insert into DataFrame. Each row represents a molecule's Avalon fingerprint and
    # each column an individual bit
    buchwald_amin_ava_fpts_dataset: pd.DataFrame = pd.DataFrame(
        buchwald_amin_ava_fpts,
        columns=['Col_A_{}'.format(bit + 1) for bit in range(buchwald_amin_ava_fpts.shape[1])]
    )

    """
    2.2 10-Fold Cross-Validation of LightGBM Regressor Model using Avalon Fingerprints of Buchwald-Hartwig 
        Cross-Coupling Components and Cross-Coupling Yields
    """
    # Start timer to calculate total execution time
    start_time: float = time.time()

    # Create a cross-validator that randomly shuffles and splits data into training and test sets
    cn_coupling_yields_cross_validator: ShuffleSplit = ShuffleSplit(
        n_splits=10,
        test_size=0.3,
        random_state=42
    )

    # Define scoring metrics.
    cn_coupling_yields_scoring_metrics: List[str] = ['r2', 'neg_root_mean_squared_error']

    # LightGBM regressor cross-validation execution.
    cn_coupling_yields_cross_validation_scores: Dict = cross_validate(
        lgbm_regressor,  # The LightGBM model to evaluate
        buchwald_amin_ava_fpts_dataset,  # The feature matrix
        buchwald_amin_yield_dataset.Output,  # The target variable (observed C-N cross-coupling yields)
        scoring=cn_coupling_yields_scoring_metrics,  # The list of metrics to evaluate
        cv=cn_coupling_yields_cross_validator  # The cross-validator
    )

    print(cn_coupling_yields_cross_validation_scores)

    end_time: float = time.time()
    execution_time: float = end_time - start_time

    print(
        '\n C-N Cross-Coupling Model Cross-Validation Execution Time: ',
        round(execution_time / 60, 2), 'min'
    )

    # Output mean values of R² and MAE for 10-fold cross-validation of LightGBM regressor
    print(
        'C-N Cross-Coupling Model Cross-Validation R²: ',
        round(np.mean(cn_coupling_yields_cross_validation_scores['test_r2']), 2),
        '\n'
    )

    print(
        'C-N Cross-Coupling Model Cross-Validation Mean Absolute Error: ',
        round(
            np.mean(-cn_coupling_yields_cross_validation_scores['test_neg_root_mean_squared_error']),
            2
        )
    )

    """
    2.3 Split Avalon Fingerprints and Cross-Coupling Yield Dataset Into Training Data & Testing Data to Train & Test 
        LightGBM Regressor Model
    """
    # Split data into training and test sets using `train_test_split()` function
    (cn_coupling_avalon_fpts_x_train, cn_coupling_avalon_fpts_x_test,
     cn_coupling_yields_y_train, cn_coupling_yields_y_test) = train_test_split(
        buchwald_amin_ava_fpts_dataset,
        buchwald_amin_yield_dataset.Output,
        test_size=0.3,
        random_state=42
    )

    # LightGBM model training
    ava_cn_coupling_yields_lgbm_model: LGBMRegressor = lgbm_regressor.fit(
        cn_coupling_avalon_fpts_x_train,
        cn_coupling_yields_y_train
    )

    # Test trained LightGBM model using testing feature matrix (Avalon fingerprints)
    # Outputs a predicted target vector (predicted HOMO-LUMO energy gap)
    cn_coupling_yields_model_predict: np.array = ava_homo_lumo_lgbm_model.predict(
        cn_coupling_avalon_fpts_x_test
    )

    # Measure model performance with R score between actual and predicted HOMO-LUMO gap values
    # Slice correlation between predicted labels and actual labels from correlation matrix
    cn_coupling_yields_predict_r: np.float64 = np.corrcoef(
        cn_coupling_yields_y_test,
        cn_coupling_yields_model_predict
    )[0, 1]
    print(
        '\n LightGBM Regressor C-N Cross-Coupling Reaction Yield Prediction R: ',
        round(cn_coupling_yields_predict_r, 2)
    )

    # Measure model performance with R² score between actual and predicted HOMO-LUMO gap values
    cn_coupling_yields_predict_r2: float = r2_score(
        cn_coupling_yields_y_test,
        cn_coupling_yields_model_predict
    )
    print(
        '\n LightGBM Regressor C-N Cross-Coupling Reaction Yield Prediction R²: ',
        round(cn_coupling_yields_predict_r2, 2)
    )

    # Measure model performance by calculating MAE between actual and predicted HOMO-LUMO gap values
    cn_coupling_yields_predict_mae: np.float64 = mean_absolute_error(
        cn_coupling_yields_y_test,
        cn_coupling_yields_model_predict
    )
    print(
        '\n LightGBM Regressor C-N Cross-Coupling Reaction Yield Prediction MAE: ',
        round(cn_coupling_yields_predict_mae, 2)
    )


if __name__ == '__main__':
    main()
