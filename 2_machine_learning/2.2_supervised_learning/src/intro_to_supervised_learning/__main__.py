#!/usr/bin/env python3
import sys
import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import jax
import rdkit
import rdkit.Chem
import rdkit.Chem.Draw
from matplotlib.axes import Axes
from rdkit.Chem import Mol
from pandas import DataFrame, Index
from typing import List, Callable, Union, Tuple
from jaxlib.xla_extension import ArrayImpl

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from utils import linear_regression_utils

matplotlib.use('TkAgg')


def main():
    """
    1) Load AqSolDB Dataset
    """
    sup_learn_data_dir: str = os.path.join(os.path.dirname(__file__), '../../data')
    aq_sol_dataset: DataFrame = pd.read_csv(os.path.join(sup_learn_data_dir, 'raw/aqsoldb-solubility-dataset.csv'))
    aq_sol_dataset.head()

    """
    2) Data Exploration
    """
    # Plot histogram with density=True to get the probability density
    plt.hist(aq_sol_dataset.Solubility, bins=30, density=True, alpha=0.6, color='g', edgecolor='black')

    # Overlay a KDE line
    sns.kdeplot(aq_sol_dataset.Solubility, color='blue')

    # Add labels and title
    plt.xlabel('Solubility/LogS')
    plt.ylabel('Probability Density')
    plt.title('AqSolDB Dataset Solubility Probability Density and Kernel Density Estimate Overlay')

    # Show the plot
    plt.show()

    # Get 3 lowest and 3 highest solubilities
    aq_sol_dataset_sorted: DataFrame = aq_sol_dataset.sort_values('Solubility')
    three_highest_lowest_sol: DataFrame = pd.concat([aq_sol_dataset_sorted[:3], aq_sol_dataset_sorted[-3:]])

    # Create list of strings for legend text
    legend_text: List = [
        f"{compound.ID}: solubility = {compound.Solubility:.2f}" for compound in three_highest_lowest_sol.itertuples()
    ]

    # Plot compounds on a grid
    three_highest_lowest_sol_mols: List[Mol] = [
        rdkit.Chem.MolFromInchi(inchi) for inchi in three_highest_lowest_sol.InChI
    ]
    rdkit.Chem.Draw.MolsToGridImage(
        three_highest_lowest_sol_mols,
        molsPerRow=3,
        subImgSize=(250, 250),
        legends=legend_text
    )

    """
    3) Feature Correlation
    """
    aq_sol_dataset_features: int = list(aq_sol_dataset.columns).index("MolWt")
    aq_sol_dataset_feature_names: Index = aq_sol_dataset.columns[aq_sol_dataset_features:]

    fig, axes = plt.subplots(nrows=5, ncols=4, sharey=True, figsize=(12, 8), dpi=300)

    # Flatten 5x4 grid of subplots into a 1D array for easier iteration
    axes: np.ndarray[Axes] = axes.flatten()

    for key, feature_name in enumerate(aq_sol_dataset_feature_names):
        ax: Axes = axes[key]
        ax.scatter(
            aq_sol_dataset[feature_name],
            aq_sol_dataset.Solubility,
            s=6,
            alpha=0.4,
            color=f"C{key}"
        )
        if key % 4 == 0:
            ax.set_ylabel('Solubility')
        ax.set_xlabel(feature_name)

    for i in range(len(aq_sol_dataset_feature_names), len(axes)):
        fig.delaxes(axes[i])

    plt.tight_layout()
    plt.show()

    """
    4) Linear Regression
    """
    linear_regression_utils.linear_regression(
        feature_vector=np.array([1, 0, 2.5]),
        weight_vector=np.array([0.2, -0.5, 0.4]),
        bias_term=4.3
    )

    """
    4.1) Calculating weight vector and bias term
    """
    aq_sol_feature_values: np.ndarray = aq_sol_dataset.loc[:, aq_sol_dataset_feature_names].values
    aq_sol_label_values: np.ndarray = aq_sol_dataset.Solubility.values
    aq_sol_feature_dim: int = aq_sol_feature_values.shape[1]

    weight_vector: Union[np.ndarray, ArrayImpl] = np.random.normal(size=aq_sol_feature_dim)
    bias_term: Union[float, ArrayImpl] = 0.0

    linear_regression_utils.mean_squared_error_loss_wrapper(
        weight_vector,
        bias_term,
        aq_sol_feature_values,
        aq_sol_label_values
    )

    """
    4.2) Gradient Descent
    """

    """
    4.2.1) Calculate Gradient with Respect to Weight Vector and Bias Term
    """
    loss_grad: Callable = jax.grad(linear_regression_utils.mean_squared_error_loss_wrapper, (0, 1))

    loss_grad(
        weight_vector,
        bias_term,
        aq_sol_feature_values,
        aq_sol_label_values
    )

    """
    4.3) Training Curve (Learning Curve)
    """
    loss_progress: List = []
    learning_rate: float = 1e-6
    for i in range(100):
        grad: Tuple = loss_grad(
            weight_vector,
            bias_term,
            aq_sol_feature_values,
            aq_sol_label_values
        )

        weight_vector -= learning_rate * grad[0]
        bias_term -= learning_rate * grad[1]

        loss_progress.append(
            linear_regression_utils.mean_squared_error_loss_wrapper(
                weight_vector,
                bias_term,
                aq_sol_feature_values,
                aq_sol_label_values
            )
        )

    plt.plot(loss_progress)
    plt.xlabel("Step")
    plt.yscale("log")
    plt.ylabel("Loss")
    plt.title("AqSolDB Training Curve")
    plt.show()

    """
    4.4) Batching - Stochastic Gradient Descent (SGD)
    """
    weight_vector: Union[np.ndarray, ArrayImpl] = np.random.normal(size=aq_sol_feature_dim)
    bias_term: Union[float, ArrayImpl] = 0.0

    loss_progress: List = []
    learning_rate: float = 1e-6
    batch_size: int = 32
    num_data_points: int = len(aq_sol_label_values)
    num_batch_data: int = num_data_points // batch_size * batch_size

    aq_sol_batched_features: np.ndarray = aq_sol_feature_values[:num_batch_data].reshape(
        (-1, batch_size, aq_sol_feature_dim)
    )
    aq_sol_batched_labels: np.ndarray = aq_sol_label_values[:num_batch_data].reshape(
        (-1, batch_size)
    )

    indices: np.ndarray = np.arange(num_batch_data // batch_size)
    np.random.shuffle(indices)

    for i in indices:
        grad: Tuple = loss_grad(
            weight_vector,
            bias_term,
            aq_sol_batched_features[i],
            aq_sol_batched_labels[i]
        )

        weight_vector -= learning_rate * grad[0]
        bias_term -= learning_rate * grad[1]

        if i % 10 == 0:
            loss_progress.append(
                linear_regression_utils.mean_squared_error_loss_wrapper(
                    weight_vector,
                    bias_term,
                    aq_sol_feature_values,
                    aq_sol_label_values
                )
            )

    plt.plot(np.arange(len(loss_progress)) * 10, loss_progress)
    plt.xlabel("Step")
    plt.yscale("log")
    plt.ylabel("Loss")
    plt.title("Batched AqSolDB Training Curve")
    plt.show()

    """
    4.5) Standardise Features
    """
    aq_sol_feature_std: np.ndarray = np.std(aq_sol_feature_values, axis=0)
    aq_sol_feature_mean: np.ndarray = np.mean(aq_sol_feature_values, axis=0)
    aq_sol_std_features: np.ndarray = (aq_sol_feature_values - aq_sol_feature_mean) / aq_sol_feature_std

    weight_vector: Union[np.ndarray, ArrayImpl] = np.random.normal(scale=0.1, size=aq_sol_feature_dim)
    bias_term: Union[float, ArrayImpl] = 0.0

    num_epochs: int = 3
    loss_progress: List = []
    learning_rate: float = 1e-2
    batch_size: int = 32
    num_data_points: int = len(aq_sol_label_values)
    num_batch_data: int = num_data_points // batch_size * batch_size

    aq_sol_batched_features: np.ndarray = aq_sol_std_features[:num_batch_data].reshape(
        (-1, batch_size, aq_sol_feature_dim)
    )
    aq_sol_batched_labels: np.ndarray = aq_sol_label_values[:num_batch_data].reshape(
        (-1, batch_size)
    )

    indices: np.ndarray = np.arange(num_batch_data // batch_size)

    for epoch in range(num_epochs):
        np.random.shuffle(indices)

        for i in indices:
            grad: Tuple = loss_grad(
                weight_vector,
                bias_term,
                aq_sol_batched_features[i],
                aq_sol_batched_labels[i]
            )

            weight_vector -= learning_rate * grad[0]
            bias_term -= learning_rate * grad[1]

            if i % 50 == 0:
                loss_progress.append(
                    linear_regression_utils.mean_squared_error_loss_wrapper(
                        weight_vector,
                        bias_term,
                        aq_sol_std_features,
                        aq_sol_label_values
                    )
                )

    plt.plot(np.arange(len(loss_progress)) * 50, loss_progress)
    plt.xlabel("Step")
    plt.yscale("log")
    plt.ylabel("Loss")
    plt.title("Standardised Batched AqSolDB Training Curve")
    plt.show()

    """
    4.6) Analysing Model Performance
    """

    """
    4.6.1) Parity Plot
    """
    linear_regression_predicted_labels: ArrayImpl = linear_regression_utils.linear_regression(
        aq_sol_std_features,
        weight_vector,
        bias_term
    )

    plt.plot([-100, 100], [-100, 100])
    plt.scatter(aq_sol_label_values, linear_regression_predicted_labels, s=4, alpha=0.7)
    plt.xlabel(r"Measured Solubility $y$")
    plt.ylabel(r"Predicted Solubility $\hat{y}$")
    plt.xlim(-13.5, 2)
    plt.ylim(-13.5, 2)
    plt.show()

    """
    4.6.2) Correlation Coefficient
    """
    linear_regression_corr_coeff = np.corrcoef(aq_sol_label_values, linear_regression_predicted_labels)[0, 1]


if __name__ == '__main__':
    main()
