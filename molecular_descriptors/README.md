# Molecular Descriptors/Representations in Machine Learning Models

## Introduction

Molecular descriptors are an **abstract representations of certain structural features of a molecule**. The majority of molecular descriptors can be classified according to their **"dimensionality"**, which refers to the representation of molecules from which descriptor values are computed **<sup>1</sup>**. This includes:

1. **0D Molecular Descriptors**:
   * This is the simplest molecular representation and is independent of any knowledge concerning the molecular structure
   * Examples include atomic mass, atomic charge, covalent & VDW radii, atomic polarizability, electronegativities & hydrophobic atomic constants **<sup>2</sup>**
2. **1D Molecular Descriptors**:
   * This is a substructure list representation and consists of a list of **structural fragments** of a molecule
   * The list of fragments can be **functional groups**, **substituents of interest** or **fingerprints**. Therefore, complete knowledge of the molecules structure is not required **<sup>2</sup>**
   * A common approach is to **encode** this list of molecular fragments into a **bit string**, which encodes the **presence or absence of a certain structural fragment**. This is called a **1D molecular fingerprint** **<sup>3</sup>** 

        <br>
        <div align="center">
          <img src="https://github.com/c-vandenberg/chemistry-machine-learning/assets/60201356/f043aca7-a03b-42b3-803e-fac864c02d5e", alt="fragment-based-1d-fingerprint" width=500/>
        </div>
        
        **Fig 1** Fragment-based 1D molecular fingerprint
   
   * The most popular 1D molecular fingerprint methods is the **Morgan Fingerprint** **<sup>4</sup>**, also known as the **extended-connectivity fingerprint (ECFP4)**

3. **2D Molecular Descriptors**:
   * 2D molecular descriptors consider how the atoms are connected/define the molecular representation based on the connectivity of atoms
   * The most commonly used 2D molecular representation approach is a **simplified molecular-input line-entry system (SMILES) string** **<sup>5</sup>**
   * A SMILES string is linear notation that encodes connectivity, structural & geometric properties of a molecule

        <br>
        <div align="center">
          <img src="", alt="" width=500/>
        </div>

4. **3D Molecular Descriptors**:
   * 3D molecular descriptors define the molecular representation not only as the connectivity of the atoms, but also as the spatial configuration of the molecule **<sup>2</sup>**
   * A popular 3D descriptor approach is 3D pharmacophore-type representations where features (e.g. hydrophobic regions or hydrogen bond donors) known or thought to be responsible for biological activity are mapped to positions on the molecule **<sup>1</sup>**

## Molecular Fingerprints

## References

**[1]** Xue, L. and Bajorath, J. (2000) ‘Molecular descriptors in Chemoinformatics, computational combinatorial chemistry, and virtual screening’, *Combinatorial Chemistry & High Throughput Screening*, 3(5), pp. 363–372.
**[2]** Todeschini, R., Consonni, V. and Gramatica, P. (2009) ‘Chemometrics in Qsar’, *Comprehensive Chemometrics*, pp. 129–172.
**[3]** Kim, J. et al. (2021) ‘Comprehensive survey of recent drug discovery using Deep Learning’, *International Journal of Molecular Sciences*, 22(18), p. 9983.
**[4]** Morgan, H.L. (1965) ‘The generation of a unique machine description for chemical structures-a technique developed at Chemical Abstracts Service.’, *Journal of Chemical Documentation*, 5(2), pp. 107–113.
**[5]** Weininger, D. (1988) ‘Smiles, a chemical language and information system. 1. introduction to methodology and encoding rules’, *Journal of Chemical Information and Computer Sciences*, 28(1), pp. 31–36.
