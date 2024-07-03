# 2.2.2 Classification Algorithms in Supervised Learning

## 1) Decision Trees

**Decision Trees** is a supervised learning method for **classification of discrete input data** (although it can also be used in **prediction of continuous input data via regression**).

A decision tree is a **hierarchical, tree data structure** consisting of **a root node**, **branches**, **internal nodes** and **leaf nodes**. 

As illustrated in **Fig 1**, the decision tree **starts with a root node**, which has **no incoming branches/parent nodes**. 
1. The **outgoing branches** from the root node then **feed into the interal nodes**, also known as **decision nodes**.
2. Based on the **available features/input vectors**, both internal decision nodes **conduct evaluations** to form **homogenous data subsets**, called **leaf nodes**.
3. The last leaf nodes along a branch (i.e. one with **no child nodes/leaf nodes below it** is called a **terminal node**.
4. The leaf nodes represent **all the possible outcomes within the dataset**. **<sup>1</sup>**

This decision tree can then **predict the result of an input vector based on decision rules inferred from the features present in the data**. **<sup>2</sup>**

  <br>
    <div align="center" width=300>
      <img src="https://github.com/c-vandenberg/machine-learning-in-drug-discovery/assets/60201356/40ef9fcf-4b24-4eef-bc02-8ccbea5945be", alt="decision-tree-generic"
        width="750"/>
      <p>
        <b>Fig 1</b> Decision tree structure. <b><sup>1</sup></b>
      </p>
    </div>
  <br>

Decision trees are useful because they are **easy to visualise, validate and audit**, so that you can **understand the factors that lead to a result**.

As an example (**Fig 2**), we can imagine trying to assess **whether or not to go surf**. One could use the following decision rules to make a choice:

  <br>
    <div align="center">
      <img src="https://github.com/c-vandenberg/machine-learning-in-drug-discovery/assets/60201356/f93c7440-7af4-42b2-89ee-0e4fae7bbd75", alt="decision-tree-example"
        width="750"/>
      <p>
        <b>Fig 2</b> Decision tree example schematic. <b><sup>1</sup></b>
      </p>
    </div>
  <br>

1. Decision tree learning employs a **divide and conquer strategy** by conducting a **greedy search** to **identify the optimal split points within the tree**.
    * A **greedy search algorithm** follows the **problem-sovling heuristic** of making the **locally optimal choice at each stage**. **<sup>3</sup>**
    * In many cases, a greedy search **does not produce an globally optimal solution**, however they can yield a **locally optimal solution** that can **approximate a globally optimal solution** in a reasonable amount of time.

2. The process of splitting the data is then **repeated in a top-down, recursive manner** until all, or the majority of the input vectors have been **classified under specific class labels**.

3. **Smaller decision trees** are **more easily able to attain pure leaf nodes** (i.e. data points in a **single class**).

4. However, **as the tree grows in size**, it becomes **increasingly difficult to maintain this purity**


There are **two types of models** for decision trees:
1. **Classification Trees** - This is where the independent variables/input vectors are **discrete values** and the **leaves of the decision tree represent class labels**. An example of this is shwon in **Fig 1**.
   
2. **Regression Trees** - This is where the independent variables/input vectors are **continuous values**.

The idea behind decision tree machine learning is that you **use a data set to train the decision tree**, which then **build a model from the data**. This trained decision tree can then be used for **decision-making with unseen data** by **walking down the tree based on the current test vector until a lead is encountered**. **<sup>1</sup>**

There are numerous algorithms for decision tree learning:
1.**Iterative Dichotomiser 3 (ID3)**:
  * One of the earliest algorithms for decision tree learning was **ID3** developed by *Quinlan et al.*, **<sup>4</sup>** which works by **splitting the data set into two separate data sets based on a single field in the vector**.
  * This field is selected by **calculating the entropy (a measure of the distribution) of the values in that field**.
  * The goal is to select a field from the vector that will **result in a decrease in entropy in subsequent splits of the data set** as the **tree is built**.
2. **C4.5**:
  * An extension of ID3 is **C4.5**, also developed *Quinlan et al.*. **<sup>5</sup>**
  * C4.5 builds decision trees in the same way as ID3 (i.e. using the concept of **information entropy** and its **reduction with each data set split**).

## 2) Random Forest

**Random Forest** is another supervised machine learning algorithm that can be used for **both classification and regression** (like decision trees).

The **"forest"** references a **collection of uncorrelated decision trees**, which are then **merged together** to **reduce variance (data spread)** and **create more accurate data predictions**. In this way, random forest is an **ensemble learning method** as it uses **multiple machine learning algorithms (multiple decision trees) to obtain a better predictive performance** than could be obtained from any of the constituent learning algorithms alone. **<sup>6</sup>** **<sup>7</sup>** **<sup>8</sup>**

Randopm

## References
**[1]** What is a decision tree? (2021) *IBM*. Available at: https://www.ibm.com/topics/decision-trees (Accessed: 04 July 2024).
**[2]** Jones, T. (2017) 'Models for machine learning - Decision Trees', *IBM developer*. Available at: https://developer.ibm.com/articles/cc-models-machine-learning/#decision-trees2 (Accessed: 02 July 2024). <br><br>
**[3]** Black, P.E. (2005) 'greedy algorithm', *Dictionary of Algorithms and Data Structures*, (accessed 4th July 2024). Available at: https://www.nist.gov/dads/HTML/greedyalgo.html
**[4]**Quinlan, J.R. (1986) ‘Induction of Decision Trees’, *Machine Learning*, 1(1), pp. 81–106. <br><br>
**[5]** Salzberg, S.L. (1994) ‘C4.5: Programs for Machine Learning by J. Ross Quinlan. Morgan Kaufmann Publishers, Inc., 1993’, *Machine Learning*, 16(3), pp. 235–240. <br><br>
**[6]** Opitz, D. and Maclin, R. (1999) ‘Popular Ensemble Methods: An Empirical Study’, *Journal of Artificial Intelligence Research*, 11, pp. 169–198. <br><br>
**[7]** Polikar, R. (2006) ‘Ensemble based systems in decision making’, *IEEE Circuits and Systems Magazine*, 6(3), pp. 21–45. <br><br>
**[8]** Rokach, L. (2009) ‘Ensemble-based classifiers’, *Artificial Intelligence Review*, 33(1–2), pp. 1–39. <br><br>
