# 1. Molecular Representations/Descriptors in Machine Learning-Based Drug Development

## Abstract

With the advancement of deep learning (DL) technology, the growth of drug-related data, and the proliferation of user-friendly DL frameworks in popular programming languages, **<sup>1</sup>** **<sup>2</sup>** methodologies based on machine learning (ML) algorithms are becoming more ubiquitous in all steps of the drug discovery and drug development process. **<sup>3</sup>**

Data quality and data representation have a significant impact on the performance of ML-based predictive models, as both heavily contribute to the quality of pre-training. Therefore, there has been a surge of interest in research on molecular representation. This includes pre-computed/fixed molecular representations such as **molecular graph representations**, **linear notations** (e.g. SMILES and molecular fingerprints), **<sup>4</sup>** **<sup>5</sup>** as well as **learned molecular representations**. **<sup>6</sup>**.

This literature review serves as an introduction to the different molecular representation approaches in ML-based drug development.

## Contents
[Abstract](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations#abstract)<br>
1.1 [Molecular Graph Theory](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#11-molecular-graph-theory)<br>
  1.1.1 [Introduction To The Molecular Graph Representation](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#111-introduction-to-the-molecular-graph-representation)<br>
  1.1.2 [Mathematical Defintion of a Graph](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#112-mathematical-defintion-of-a-graph)<br>
  1.1.3 [Graph Traversal Algorithms](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#113-graph-traversal-algorithms)<br>
  1.1.4 [Molecular Graph Reprentations](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#114-molecular-graph-reprentations)<br>
  1.1.5 [Advantages of Molecular Graph Representations](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#115-advantages-of-molecular-graph-representations)<br>
  1.1.6 [Disadvantages of Molecular Graph Representations](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#116-disadvantages-of-molecular-graph-representations)<br>
  1.1.7 [Molecular Graphs in AI-Driven Small Molecule Drug Discovery](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#117-molecular-graphs-in-ai-driven-small-molecule-drug-discovery)<br>
  1.1.8 [References](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.1_molecular_graph_theory#118-references)<br>
1.2 [Molecular Descriptors](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#12-molecular-descriptors-linear-notations-in-machine-learning-models)<br>
  1.2.1 [Introduction to Molecular Descriptors](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#121-introduction-to-molecular-descriptors-linear-notations)<br>
  1.2.2 [Molecular Fingerprints](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#122-molecular-fingerprints)<br>
  1.2.3 [Key-Based Molecular Fingerprints - MACCS Keys](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#123-key-based-molecular-fingerprints---maccs-keys)<br>
  1.2.4 [Hash-Based Molecular Fingerprints - Daylight Fingerprint & ECFPs](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#124-hash-based-molecular-fingerprints---daylight-fingerprint--ecfps)<br>
  1.2.5 [Advantages & Applications of Molecular Fingerprints](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#125-advantages--applications-of-molecular-fingerprints)<br>
  1.2.6 [Molecular Fingerprints in Machine Learning](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#126-molecular-fingerprints-in-machine-learning)<br>
  1.2.7 [References](https://github.com/c-vandenberg/machine-learning-in-drug-discovery/tree/master/1_molecular_representations/1.2_molecular_descriptors#127-references)<br>

Additionally, each section of the literature review will have accompanying Python exercises related to the topic in question.

## References

**[1]** Abadi, M. *et al.* (2015) ‘TensorFlow: Large-Scale Machine Learning on Heterogeneous Systems’, https://www.tensorflow.org/, Software available from tensorflow.org. <br><br>
**[2]** Paszke, A. *et al.* (2017) ‘NIPS Autodiff Workshop’. <br><br>
**[3]** Kim, J. *et al.* (2021) ‘Comprehensive survey of recent drug discovery using Deep Learning’, *International Journal of Molecular Sciences*, 22(18), p. 9983. <br><br>
**[4]** Rifaioglu, A.S. et al. (2020) ‘DEEPScreen: High performance drug–target interaction prediction with convolutional neural networks using 2-D structural compound representations’, *Chemical Science*, 11(9), pp. 2531–2557.<br><br>
**[5]** David, L. et al. (2020) ‘Molecular representations in AI-Driven Drug Discovery: A review and practical guide’, *Journal of Cheminformatics*, 12(1).<br><br>
**[6]** Yang, K. et al. (2019) ‘Analyzing learned molecular representations for property prediction’, *Journal of Chemical Information and Modeling*, 59(8), pp. 3370–3388.<br><br>
