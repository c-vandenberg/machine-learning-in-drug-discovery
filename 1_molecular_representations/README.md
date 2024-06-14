# 1 Molecular Representations/Descriptors in Machine Learning-Based Drug Development

## Abstract

With the advancement of deep learning (DL) technology, the growth of drug-related data, and the proliferation of user-friendly DL frameworks in popular programming langauages, **<sup>1</sup>** **<sup>2</sup>** methodologies based on machine learning (ML) algorithms are becoming more ubiquitous in all steps of the drug discovery & drug development process. **<sup>3</sup>**

Data quality and data representation have a significant impact on the performance of ML-based predictive models, as both heavily contribute to the quality of pre-training. Therefore, there has been a surge of interest in research on molecular representation. This includes pre-computed/fixed molecular representations such as **molecular graph representions**, **linear notations** (e.g. SMILES & molecular fingerprints), **<sup>4</sup>** **<sup>5</sup>** as well as **learned molecular representations**. **<sup>6</sup>**.

This literature review serves as an introduction to the different molecular representation approaches in ML-based drug development.

## Contents

1.1 **Molecular Graph Theory**<br><br>
1.2 **Molecular Descriptors**:
   * 1.2.1 **Introduction to Molecular Descriptors**:
   * 1.2.2 **Molecular Fingerprints**
   * 1.2.3 **Key-Based Molecular Fingerprints - MACCS Keys**
   * 1.2.4 **Hash-Based Molecular Fingerprints - Daylight Fingerprint & ECFPs**
   * 1.2.5 **Advantages & Applications of Molecular Fingerprints**
   * 1.2.6 **Molecular Fingerprints in Machine Learning**
   * 1.2.7 **References**

Additionally, each section of the literature review will have accompanying Python exercises related to the topic in question.

## 1.2 References

**[1]** Abadi, M. *et al.* (2015) ‘TensorFlow: Large-Scale Machine Learning on Heterogeneous Systems’, https://www.tensorflow.org/, Software available from tensorflow.org. <br><br>
**[2]** Paszke, A. *et al.* (2017) ‘NIPS Autodiff Workshop’. <br><br>
**[3]** Kim, J. *et al.* (2021) ‘Comprehensive survey of recent drug discovery using Deep Learning’, *International Journal of Molecular Sciences*, 22(18), p. 9983. <br><br>
**[4]** Rifaioglu, A.S. et al. (2020) ‘DEEPScreen: High performance drug–target interaction prediction with convolutional neural networks using 2-D structural compound representations’, *Chemical Science*, 11(9), pp. 2531–2557.<br><br>
**[5]** David, L. et al. (2020) ‘Molecular representations in AI-Driven Drug Discovery: A review and practical guide’, *Journal of Cheminformatics*, 12(1).<br><br>
**[6]** Yang, K. et al. (2019) ‘Analyzing learned molecular representations for property prediction’, *Journal of Chemical Information and Modeling*, 59(8), pp. 3370–3388.<br><br>
