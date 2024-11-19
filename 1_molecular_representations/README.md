# 1. Molecular Representations/Descriptors in Machine Learning-Based Drug Development

## Abstract

The **representation of molecules** has been a central focus in chemistry throughout history, evolving from **traditional structural diagrams** depicting atoms and bonds, to **advanced machine-readable notations** essential for cheminformatics and drug discovery. 

This review explores various chemical representations, including **structural encodings** such as **molecular graphs**, and **linear notations** like **SMILES**, and **molecular fingerprints**. It examines the strengths and limitations of these representations, and highlights the development of **computer-readable formats** that facilitate **efficient digital storage**, **querying**, and **visualisation** of chemical compounds. **<sup>1</sup>** Additionally, the review presents case studies on widely used molecular representations and descriptors, including **MACCS-keys fingerprints**, **Avalon fingerprints**, and **Morgan fingerprints**, as well as their applications in **compound querying and clustering**, such as **Taylor-Butina clustering**.

## Contents
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

Additionally, each section of the literature review will have accompanying Python exercises related to the topic in question.

## References

**[1]** David, L. et al. (2020) ‘Molecular representations in AI-Driven Drug Discovery: A review and practical guide’, *Journal of Cheminformatics*, 12(1).<br><br>
