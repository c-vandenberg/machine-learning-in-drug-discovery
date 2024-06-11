# Molecular Descriptors/Representations in Machine Learning Models

## Introduction

Molecular descriptors are an **abstract representations of certain structural features of a molecule**. The majority of molecular descriptors can be classified according to their **"dimensionality"**, which refers to the representation of molecules from which descriptor values are computed **<sup>1</sup>**. This includes:

1. **0D Molecular Descriptors**:
   * This is the simplest molecular representation and is independent of any knowledge concerning the molecular structure
   * Examples include atomic mass, atomic charge, covalent & VDW radii, atomic polarizability, electronegativities & hydrophobic atomic constants **<sup>2</sup>**
2. **1D Molecular Descriptors**:
   * This is a substructure list representation and consists of a list of **structural fragments** of a molecule
   * The list of fragments can be **functional groups** or **substituents of interest**. Therefore, complete knowledge of the molecules structure is not required **<sup>2</sup>**
   * A common approach is to **encode** this list of molecular fragments into a **bit string**, which encodes the **presence or absence of a certain structural fragment**. **<sup>3</sup>** 

        <br>
        <div align="center">
          <img src="", alt=""/>
        </div>

3. **2D Molecular Descriptors** - 
4. **3D Molecular Descriptors** - 

## References

**[1]** Xue, L. and Bajorath, J. (2000) ‘Molecular descriptors in Chemoinformatics, computational combinatorial chemistry, and virtual screening’, *Combinatorial Chemistry & High Throughput Screening*, 3(5), pp. 363–372.
**[2]** Todeschini, R., Consonni, V. and Gramatica, P. (2009) ‘Chemometrics in Qsar’, *Comprehensive Chemometrics*, pp. 129–172.
**[3]** Kim, J. et al. (2021) ‘Comprehensive survey of recent drug discovery using Deep Learning’, *International Journal of Molecular Sciences*, 22(18), p. 9983.
Kim, J. et al. (2021) ‘Comprehensive survey of recent drug discovery using Deep Learning’, International Journal of Molecular Sciences, 22(18), p. 9983. doi:10.3390/ijms22189983. 