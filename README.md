# Machine Learning in Drug Discovery

## Abstract

Advancements in machine learning technology, the exponential growth of drug-related data, and the widespread availability of user-friendly machine learning frameworks in popular programming languages **<sup>1</sup>** **<sup>2</sup>** are making machine learning methodologies increasingly prevalent throughout all stages of the drug discovery and development process. **<sup>3</sup>**

**Data quality** and **representation** significantly impact the performance of machine learning-based predictive models, as both are crucial for **effective pre-training**. As a result, there has been a surge of research interest in molecular representation. This research encompasses **pre-computed** or **fixed molecular representations**, such as **molecular graph representations** and **linear notations** (e.g., SMILES and molecular fingerprints), **<sup>4</sup>** **<sup>5</sup>** as well as **learned molecular representations**. **<sup>6</sup>**.

This literature review provides an overview of the various molecular representation approaches used in machine learning-based drug development and explores their applications in conjunction with machine learning models for predicting molecular properties and reactions.

## Contents
1. [Molecular Representations/Descriptors in Machine Learning-Based Drug Development](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/blob/master/1_molecular_representations/README.md#1-molecular-representationsdescriptors-in-machine-learning-based-drug-development)<br>
  1.1 [Molecular Graph Theory](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#11-molecular-graph-theory)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.1.1 [Introduction To The Molecular Graph Representation](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#111-introduction-to-the-molecular-graph-representation)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.1.2 [Mathematical Defintion of a Graph](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#112-mathematical-defintion-of-a-graph)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.1.3 [Graph Traversal Algorithms](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#113-graph-traversal-algorithms)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.1.4 [Molecular Graph Reprentations](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#114-molecular-graph-reprentations)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.1.5 [Advantages of Molecular Graph Representations](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#115-advantages-of-molecular-graph-representations)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.1.6 [Disadvantages of Molecular Graph Representations](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#116-disadvantages-of-molecular-graph-representations)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.1.7 [Molecular Graphs in AI-Driven Small Molecule Drug Discovery](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#117-molecular-graphs-in-ai-driven-small-molecule-drug-discovery)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.1.8 [References](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#118-references)<br>
  1.2 [Molecular Descriptors](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#12-molecular-descriptors-linear-notations-in-machine-learning-models)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.2.1 [Introduction to Molecular Descriptors](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#121-introduction-to-molecular-descriptors-linear-notations)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.2.2 [Molecular Fingerprints](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#122-molecular-fingerprints)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.2.3 [Key-Based Molecular Fingerprints - MACCS Keys](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#123-key-based-molecular-fingerprints---maccs-keys)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.2.4 [Hash-Based Molecular Fingerprints - Daylight Fingerprint & ECFPs](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#124-hash-based-molecular-fingerprints---daylight-fingerprint--ecfps)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.2.5 [Advantages & Applications of Molecular Fingerprints](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#125-advantages--applications-of-molecular-fingerprints)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.2.6 [Molecular Fingerprints in Machine Learning](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#126-molecular-fingerprints-in-machine-learning)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1.2.7 [References](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#127-references)<br>
2. [Machine Learning-Based Drug Development](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/2_machine_learning#2-machine-learning-based-drug-development)<br>
2.1 [Introduction to Machine Learning](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/2_machine_learning/2.1_introduction_to_machine_learning#21-introduction-to-machine-learning)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;2.1.1 [How does Machine Learning Work?](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/2_machine_learning/2.1_introduction_to_machine_learning#211-how-does-machine-learning-work)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;2.1.2 [Machine Learning Methods](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/2_machine_learning/2.1_introduction_to_machine_learning#212-machine-learning-methods)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;2.1.3 [Machine Learning Notation](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/2_machine_learning/2.1_introduction_to_machine_learning#213-machine-learning-notation)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;2.1.4 [References](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/2_machine_learning/2.1_introduction_to_machine_learning#214-references)<br>
2.2 [Supervised Learning](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/2_machine_learning/2.2_supervised_learning#22-supervised-learning)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;2.2.1 [Classification Algorithms in Supervised Learning](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/2_machine_learning/2.2_supervised_learning/2.2.1_classification_algorithms#221-classification-algorithms-in-supervised-learning)<br>
  &nbsp;&nbsp;&nbsp;&nbsp;2.2.2 [Regression Algorithms in Supervised Learning](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/2_machine_learning/2.2_supervised_learning/2.2.2_regression_algorithms#222-regression-algorithms-in-supervised-learning)<br>

## References

**[1]** Abadi, M. *et al.* (2015) ‘TensorFlow: Large-Scale Machine Learning on Heterogeneous Systems’, https://www.tensorflow.org/, Software available from tensorflow.org. <br><br>
**[2]** Paszke, A. *et al.* (2017) ‘NIPS Autodiff Workshop’. <br><br>
**[3]** Kim, J. *et al.* (2021) ‘Comprehensive survey of recent drug discovery using Deep Learning’, *International Journal of Molecular Sciences*, 22(18), p. 9983. <br><br>
**[4]** Rifaioglu, A.S. et al. (2020) ‘DEEPScreen: High performance drug–target interaction prediction with convolutional neural networks using 2-D structural compound representations’, *Chemical Science*, 11(9), pp. 2531–2557.<br><br>
**[5]** David, L. et al. (2020) ‘Molecular representations in AI-Driven Drug Discovery: A review and practical guide’, *Journal of Cheminformatics*, 12(1).<br><br>
**[6]** Yang, K. et al. (2019) ‘Analyzing learned molecular representations for property prediction’, *Journal of Chemical Information and Modeling*, 59(8), pp. 3370–3388.<br><br>
