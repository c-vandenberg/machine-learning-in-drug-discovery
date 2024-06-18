# 1.1 Molecular Graph Theory

## 1.1.1 Introduction To The Molecular Graph Representation

The vast majority of the molecular representations discussed in this literature review are built on the concept of the molecular graph representation. It is therefore important to introduce the theory behind molecular graphs and gain a solid understanding of it. This will primarily involve the abstract mathematical structure/data structure of the graph itself, however it is importatnt to note that notations and file formats such as **SMILES strings** and **Molfiles** are also build using molecular graphs. **<sup>1</sup>**

Molecular graph representation is based on the concept of **mapping a molecules atoms and bonds into sets of nodes and edges**. Typically, the nodes are represented using letters indicating the atom type (as per the periodic table), or simply using points where the bonds meet (for carbon atoms). **<sup>1</sup>**

A molecular graph representation is formally a **2D object** that can be used to **represent 3D information/spatial relationships**. For example, a molecular graph can contain information such as **atomic coordinates**, **bond angles** and **chirality**. However, because graph nodes only show **pairwise relationships**, this spatial information must be **encoded as node and/or edge attributes**. **<sup>1</sup>**

Typically, in schematic representations of graphs, the nodes are represented representated using circles or spheres, and the edges using lines (**Fig 1**). The 2D and 3D representations of graphs can be visualised by many software packages, such as **ChemDraw**, **<sup>2</sup>** **Avogadro**, **<sup>3</sup>** and **VMD** **<sup>4</sup>**

<br>
<div align="center">
  <img src="https://github.com/c-vandenberg/chemistry-machine-learning/assets/60201356/ee1147a9-f6e9-4c16-9f38-f954a6be8c2e", alt="graph-data-structure" width=500/>
    <p>
      <b>Fig 1</b> Schematic representation of a graph data structure <b><sup>5</sup></b>
    </p>
</div>
<br>

This section of the literature review is primarily based on a review into 'Molecular representations in AI-driven drug discovery' by David *et al* **<sup>1</sup>**

## 1.1.2 Mathematical Defintion of a Graph

A graph is formally/mathematically defined as a **tuple** *G=(V, E)* of a **set of nodes** *V*, and a **set of edges** *E*, where each edge element *e* in set *E* (*e ∈ E*) **connects a pair of nodes v** from set *V* (*v ∈ V*).

Intuitively, we can think of set *V* as being a set of **all atoms in the molecule**, and set *E* being a set of **all bonds in the molecule**. This is typically the convention, but it **does not have to be the case**.

Molecular graphs are typically **undirected**, meaning the **pairs of nodes that define each edge are unordered**.

In order to construct a concrete, programmatic representation of a graph from this abstract mathematical concept, we need to **map the sets of nodes and edges to linear data structures**. Common data structures used to do this are **arrays** or **matrices**. A linear data structure is required in order to **specify the connectivity of the nodes**. Therefore, despite the face that the ordering in the sets is irrelevant, an **artificial node-ordering** in order to **encode a molecule using arrays or matrices**. 

The molecular information to be mapped to the graph can be:
1. How atoms are connected in the molecule
2. The identity of the atoms
3. The identity of the bonds

### How atoms are connected in the molecule - Adjacency Matrix

The most common representation for **how the atoms are connected** in a molecule is an **adjacency matrix** (sometimes known as a **connectivity matrix**, **Fig 2b**). An adjacency matrix works as follows:
1. For adjacency matrix **A**, *a<sub>ij</sub>* is an element of **A**
2. If *a<sub>ij</sub> = 1*, then a **bond exists** between nodes *v<sub>i</sub>* and *v<sub>j</sub>* in molecular graph *G*
3. If *a<sub>ij</sub> = 0*, then a **bond does notexists** between nodes *v<sub>i</sub>* and *v<sub>j</sub>* in molecular graph *G*
4. **N.B** An adjacency matrix only specifies a bonds exists, it does not specify what type of bond exists

### The Identity of Atoms - Node Features Matrix

The most common representation for **the idendity of atoms** in a molecule is a **node features matrix** (**Fig 2c**), which works as follows:
1. For node features matrix **X**, each **row** of **X** corresponds to a node *v<sub>i</sub>* (i.e. an **atom in the molecule**) in molecular graph *G*
2. This row is also referred to as the **node feature vector** *x<sub>i</sub>* for that atom
3. The **length** of *x<sub>i</sub>* corresponds to the **number of atom features** you have chosen to encode
4. For example, in **Fig 2c**, a **one-hot encoding** of atom type and formal charge has been chosen (i.e. the atom types and formal charges are represented as a **binary value** so that they can be **fed into machine-learning algorithms**)

### The Identity of Bonds

The most common representation for **the idendity of bonds** in a molecule is a **edge features matrix** (**Fig 2d**), which works as follows:
1. For edge features matrix **E**, each **row** of **E** corresponds to an edge *e<sub>ij</sub> = (v<sub>i</sub>,v<sub>j</sub>)* (i.e. the **bond between atoms i and j**) in molecular graph *G*
2. This row is also referred to as the **edge feature vector** *e<sub>ij</sub>* for that edge
3. The **length** of *e<sub>ij</sub>* corresponds to the **number of bond features** you have chosen to encode
4. For example, in **Fig 2d**, a **one-hot encoding** of possible bond types {single, double, triple, aromatic} has been chosen

<br>
<div align="center">
  <img src="https://github.com/c-vandenberg/chemistry-machine-learning/assets/60201356/0e8a939c-ee30-41d0-adf0-70720dd7f069", alt="molecular-graphs-types"/>
    <p>
      <b>Fig 2</b> Example graph representation for acetic acid. <b>a</b> Graph representation of acetic acid with nodes numbered from one to four. <b>b</b> Example
adjacency matrix, <b>A</b>, for an acetic acid graph with the corresponding node ordering on the left. <b>c</b> Example node features matrix, <b>X</b>, which one-hot
encodes a few selected properties. <b>d</b> Example edge features matrix, <b>E</b>, where each edge feature vector is a one-hot encoding of single, double, or
triple bonds. “Implicit Hs” stands for the number of implicit hydrogens on a given node <b><sup>1</sup></b>
    </p>
</div>
<br>

## References
**[1]** David, L. *et al*. (2020) ‘Molecular representations in AI-Driven Drug Discovery: A review and practical guide’, *Journal of Cheminformatics*, 12(1).<br><br>
**[2]** ChemDraw. PerkinElmer Informatics.<br><br>
**[3]** Marcus, D.H. (2012) Avogadro: an advanced semantic chemical editor, visualization, and analysis platform. *J Cheminform*. 4:17.<br><br>
**[4]** Humphrey, W., Dalke, A and Schulten K (1996) VMD: visual molecular dynam-ics. *J Mol Graph* 14(1):33–38.<br><br>
**[5]** What is graph data structure? (2023) *GeeksforGeeks*. Available at: https://www.geeksforgeeks.org/what-is-graph-data-structure/ (Accessed: 18 June 2024). 
