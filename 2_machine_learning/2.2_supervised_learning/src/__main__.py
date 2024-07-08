import os
import sys
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import jax.numpy as jnp
import jax
import sklearn.manifold
import sklearn.cluster
import rdkit
import rdkit.Chem
import rdkit.Chem.Draw
from matplotlib.axes import Axes
from rdkit.Chem import Mol
from pandas import DataFrame, Index
from typing import List
from utils import linear_regression_utils
from jax.example_libraries import optimizers

matplotlib.use('TkAgg')


def main():
    """
    1) Load AqSolDB Dataset
    """
    sup_learn_data_dir: str = os.path.join(os.path.dirname(__file__), '../data')
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


if __name__ == '__main__':
    main()
