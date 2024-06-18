# 1.1 Molecular Graph Theory

## 1.1.1 Introduction To The Molecular Graph Representation

The vast majority of the molecular representations discussed in this literature review are built on the concept of the molecular graph representation. It is therefore important to introduce the theory behind molecular graphs and gain a solid understanding of it. This will primarily involve the abstract mathematical structure/data structure of the graph itself, however it is importatnt to note that notations and file formats suycg as **SMILES strings** and **Molfiles** are also build using molecular graphs. **<sup>1</sup>**

Molecular graph representation is based on the concept of **mapping a molecules atoms and bonds into sets of nodes and edges**. Typically, the nodes are represented using letters indicating the atom type (as per the periodic table), or simply using points where the bonds meet (for carbon atoms). **<sup>1</sup>**

A molecular graph representation is formally a **2D object** that can be used to **represent 3D information/spatial relationships**. For example, a molecular graph can contain information such as **atomic coordinates**, **bond angles** and **chirality**. However, because graph nodes only show **pairwise relationships**, this spatial information must be **encoded as node and/or edge attributes**. **<sup>1</sup>**

Typically, in schematic representations of graphs, the nodes are represented representated using circles or spheres, and the edges using lines (**Fig 1**). The 2D and 3D representations of graphs can be easily visualised by many software packages, such as **ChemDraw**, **<sup>2</sup>** **Avogadro**, **<sup>3</sup>** and **VMD** **<sup>4</sup>**

<br>
<div align="center">
  <img src="https://github.com/c-vandenberg/chemistry-machine-learning/assets/60201356/ee1147a9-f6e9-4c16-9f38-f954a6be8c2e", alt="graph-data-structure" width=500/>
    <p>
      <b>Fig 1</b> Schematic representation of a graph data structure
    </p>
</div>
<br>

## References
**[1]** David, L. *et al*. (2020) ‘Molecular representations in AI-Driven Drug Discovery: A review and practical guide’, *Journal of Cheminformatics*, 12(1).<br><br>
**[2]** ChemDraw. PerkinElmer Informatics.<br><br>
**[3]** Marcus, D.H. (2012) Avogadro: an advanced semantic chemical editor, visualization, and analysis platform. *J Cheminform*. 4:17.<br><br>
**[4]** Humphrey, W., Dalke, A and Schulten K (1996) VMD: visual molecular dynam-ics. *J Mol Graph* 14(1):33–38.<br><br>
