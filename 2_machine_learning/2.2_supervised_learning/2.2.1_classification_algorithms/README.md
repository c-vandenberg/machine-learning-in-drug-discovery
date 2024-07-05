# 2.2.2 Classification Algorithms in Supervised Learning

## 1) Decision Trees

**Decision Trees** are a supervised learning method for **classification of discrete input data** (although it can also be used in **prediction of continuous input data via regression**).

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

### How Decision Trees Work

1. Decision tree learning employs a **divide and conquer strategy** by conducting a **greedy search** to **identify the optimal split points within the tree**.
    * A **greedy search algorithm** follows the **problem-sovling heuristic** of making the **locally optimal choice at each stage**. **<sup>3</sup>**
    * In many cases, a greedy search **does not produce an globally optimal solution**, however they can yield a **locally optimal solution** that can **approximate a globally optimal solution** in a reasonable amount of time. **<sup>1</sup>**

2. The process of splitting the data is then **repeated in a top-down, recursive manner** until all, or the majority of the input vectors have been **classified under specific class labels**. **<sup>1</sup>**

3. **Smaller decision trees** are **more easily able to attain pure leaf nodes** (i.e. data points in a **single class**). **<sup>1</sup>**

4. However, **as the tree grows in size**, it becomes **increasingly difficult to maintain this purity**, and it usually results in **too little data falling within a given subtree**. When this occurs, it is known as **data fragmentation** and can often lead to **overfitting**. **<sup>1</sup>**

5. Therefore, there is often a preference for **smaller decision trees**, in keeping with **"Occam's Razor**; that is, "**entities should not be multipled beyond necessity"**. In other words, **decision trees should add complexity only if necessary**, as the **simplest explanation is often the best**. **<sup>1</sup>**

6. To **reduce complexity and prevent overfitting, pruning is often employed**, whereby **branches are removed that split on features with low importance**. **<sup>1</sup>**

7. The model's fit can then be **evaluated through the process of cross-validation**. **<sup>1</sup>**

8. Another way that decision trees can maintain their accuracy is by **forming an ensemble via a random forest algorithm**; this classifier **predicts more accurate results**, particularly when the individual trees are **uncorrelated with each other**. **<sup>1</sup>**

### Types of Decision Trees

There are **two types of models** for decision trees:
1. **Classification Trees** - This is where the independent variables/input vectors are **discrete values** and the **leaves of the decision tree represent class labels**. An example of this is shwon in **Fig 1**.
   
2. **Regression Trees** - This is where the independent variables/input vectors are **continuous values**.

The idea behind decision tree machine learning is that you **use a data set to train the decision tree**, which then **build a model from the data**. This trained decision tree can then be used for **decision-making with unseen data** by **walking down the tree based on the current test vector until a lead is encountered**. **<sup>1</sup>**

There are numerous algorithms for decision tree learning:
1.**Iterative Dichotomiser 3 (ID3)**:
  * One of the earliest algorithms for decision tree learning was **ID3** developed by *Quinlan et al.*, **<sup>4</sup>** which works by **splitting the data set into two separate data sets based on a single field in the vector**.
  * This field is selected by **calculating the entropy of the values in that field**.
  * The goal is to select a field from the vector that will **result in a decrease in entropy in subsequent splits of the data set** as the **tree is built**.
2. **C4.5**:
  * An extension/iteration of ID3 is **C4.5**, also developed *Quinlan et al.*. **<sup>5</sup>**
  * C4.5 builds decision trees in the same way as ID3 (i.e. using the concept of **entropy and information gain** and its **reduction with each data set split**).
3. **CART**:
  * **CART** is an abbrevation for **"classification and regressino trees** and was introduced by Leo Breiman. **<sup>6</sup>**
  * This algorithm typically utilizes a **Gini impurity** to **identifiy the ideal attribute to split on**.
  * Gini impurity measures **how often a randomly chosen attribute is misclassified**. When evaluating a Gini impurity, a **lower value is more ideal**

### How to Choose the Best Attribute at Each Node

While there multiple ways to **select the best attribute at each node**, two methods in particular act as popular splitting criterion for decision tree models:
1. **Entropy and Information Gain**
2. **Gini Impurity**

#### Entropy and Information Gain

**Entropy** is a concept that stems form **information theory**, which measures the **impurity of the sample values**. It is defined by the following formula:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Ctext%7BEntropy%7D(S)%20%3D%20-%20%5Csum_%7Bc%20%5Cin%20C%7D%20p(c)%20%5Clog_2%20p(c)", alt='entropy-formula'/>
    </div>
<br>

where:
* ![dataset_entropy](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DS(X)) - is the **entropy of the current data set**.<br><br>
* ![summation_of_elements_in_class_c](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Csum_%7Bc%20%5Cin%20C%7D) - is the **summation of all elements/data points in class *C*.** <br><br>
* ![data_point_class_c_proportion](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7Dp(c)) - is the **proportion of elements/data points in class *C* to the number of total elements/data points in the total data set**.<br><br>

Entropy values can **fall between 0 and 1**:
* If **all elements/data points in the data set belong to one class**, then the **entropy value will equal 0 (lowest value)**.
* If **half of the elements/data points are classified as one class and the other half are in another class**, then the **entropy value will equal 1 (highest value)**.

In order to **select the best feature to split on** and **find the optimal decision tree**, the feature with the **attribute/feature with the smallest amount of entropy should be used**
* In the example in **Fig 3**, attributes/features would be "Outlook," "Temperature," "Humidity," and "Wind."

**Information gain** represents the **difference in entropy before and after a split** on a given feature/attribute. The feature/attribute with the **highest information gain will produce the best split** as it is doing the **best job at classifying the training data according to its target classification**. **<sup>1</sup>**

Information gain is usually represented with the following formula:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DIG(X,%20%5Calpha)%20=%20S(X)%20-%20%5Csum_{v%20%5Cin%20%5Ctext%7BValues%7D(%5Calpha)}%20%5Cfrac%7B|X_v|%7D%7B|X|%7D%20S(X_v)", alt='information-gain-formula'/>
    </div>
<br>

where:<br><br>
* ![dataset_entropy](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DS(X)) - is the **entropy of the full data set**<br><br>
* ![data_subset](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DX_v) - is a **subset of the full data set**<br><br>
* ![data_subset_proportion](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cfrac%7B|X_v|%7D%7B|X|%7D) - is the **proportion of the number of elements/data points in the subset to the number of elements/data points in the full dataset**<br><br>
* ![data_subset_entropy](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DS(X_v)) - is the **entropy of the subset**<br><br>

As an worked example:

<br>
    <div align="center">
        <img src="https://github.com/c-vandenberg/machine-learning-in-drug-discovery/assets/60201356/52b7de1d-6fb0-4d8c-b8c2-a380d315595c", alt="decision-tree-example" width="750"/>
      <p>
        <b>Fig 3</b> Decision tree example table. <b><sup>1</sup></b>
      </p>
    </div>
<br>

**1) Entropy Calculation**

For the data set in **Fig 3**, the **entropy is 0.985** for the "Tennis" feature/attribute.

The entropy is calculated by finding the **proportion of days where "Tennis" is "Yes" (4/7)**, and the **proportion of days where "Tennis" is "No" (3/7)**. These values are then plugged into the entropy formula:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20S(Tennis)%20%3D%20-%5Cleft(%5Cfrac%7B4%7D%7B7%7D%5Cright)%20%5Clog_2%5Cleft(%5Cfrac%7B4%7D%7B7%7D%5Cright)%20-%20%5Cleft(%5Cfrac%7B3%7D%7B7%7D%5Cright)%20%5Clog_2%5Cleft(%5Cfrac%7B3%7D%7B7%7D%5Cright)%20%3D%200.985", alt='play_tennis_example_entropy_calculation'/>
    </div>
<br>

**2) Information Gain Calculation for Each Feature/Attribute**

We can then **compute the information gain** for **each of the attributes individually** in order to **decide which should be used for the first split in the decision tree**

For example, the information gain for the attribute **"Wind**" would be involves calculating the **entropy for strong wind** and the **entropy for weak wind**:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20S(Strong%20Wind)%20%3D%20-%5Cleft(%5Cfrac%7B2%7D%7B7%7D%5Cright)%20%5Clog_2%5Cleft(%5Cfrac%7B2%7D%7B7%7D%5Cright)%20%3D%200.516", alt='strong_wind_example_entropy_calculation'/>
    </div>
<br>

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20S(Weak%20Wind)%20%3D%20-%5Cleft(%5Cfrac%7B5%7D%7B7%7D%5Cright)%20%5Clog_2%5Cleft(%5Cfrac%7B5%7D%7B7%7D%5Cright)%20%3D%200.347", alt='weak_wind_example_entropy_calculation'/>
    </div>
<br>

The **informtation gain** for the "Wind" attribute when the "Tennis" attribute is used to split the data can be calculated by:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DIG(Tennis%2C%20Wind)%20%3D%20S(Tennis)%20-%20%5Cleft(%5Cfrac%7B%7CStrongWind%7C%7D%7B%7CWind%7C%7D%5Cright)%20*%20S(SrongWind)%20-%20%5Cleft(%5Cfrac%7B%7CWeakWind%7C%7D%7B%7CWind%7C%7D%5Cright)%20*%20S(WeakWind)", alt='tennis_wind_example_information_gain_calculation'/>
    </div>
<br>

Plugging in the numbers:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DIG(Tennis%2C%20Wind)%20%3D%200.985%20-%20%5Cleft(%5Cfrac%7B2%7D%7B7%7D%5Cright)%20*%200.516%20-%20%5Cleft(%5Cfrac%7B5%7D%7B7%7D%5Cright)%20*%200.347%20%3D%200.590", alt='tennis_wind_example_information_gain_calculation_plugged_numbers'/>
    </div>
<br>

We would then **repeat this calculation for information gain for each attribute in Fig 3**, and **select the attribute with the highest information gain** to be the **first split point in the decision tree**.

In this case, **"Outlook" produces the highest information gain**. From there, the process is **repeated for each subtree**

#### Gini Impurity

Gini impurity is the **probability of incorrectly classifying a random data point in the data set** if it were **labelled based on the class distribution of the data set**.

Similar to entropy, Gini Coefficient values can **fall between 0 and 1**:
* If **all elements/data points in the data set belong to one class**, then the **Gini Coefficient value will equal 0 (lowest value)**.
* If **half of the elements/data points are classified as one class and the other half are in another class**, then the **Gini Coefficient value will equal 1 (highest value)**.


<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20G(X)%20%3D%201%20-%20%5Csum_%7Bi%3D1%7D%5En%20p_i%5E2", alt='gini_coefficient_formula'/>
    </div>
<br>

where:
* ![gini_coefficient](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20G(X)) - is the Gini Coefficient for a **given data set X**.
* ![element_class_probability](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20p_i) - is the **probability of an element being classified into class *i* in set *X***.
* ![number_of_classes](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20n) - is the number of classes.

## 2) Random Forest

**Random Forest** is another supervised machine learning algorithm that can be used for **both classification and regression** (like decision trees).

The **"forest"** references a **collection of uncorrelated decision trees**. The **outputs of these multiple decision tress** are then **merged together** to **reduce variance (data spread)** and **create more accurate data predictions**. 

In this way, random forest is an **ensemble learning method** as it uses **multiple machine learning algorithms (multiple decision trees) to obtain a better predictive performance** than could be obtained from any of the constituent learning algorithms alone. **<sup>7</sup>** **<sup>8</sup>** **<sup>9</sup>**

### Ensemble Learning Methods

Ensemble learning methods are made up of a **set of classifiers** (e.g. decision trees), where the **predictions of each classifier are aggregated** in order to **identify the most popular result**.

The most well-known ensemble methods are **bagging** (also known as **bootstrap bagging**), and **boosting**
1. **Bagging (Boostrap Bagging)**:
   * Bagging, developed by *Breiman et al.*, **<sup>10</sup>** is an ensemble method that involves **training multiple models independently on random subsets of the data**, and **aggregating their predictions through voting (for classification) or averaging (for regression)**.
   * Each model is trained on a random subset of the data set and **sampled with replacement**. This means that **the individual data points can be chosen more than once**.
   * This random subset is known as a **bootstrap sample**.
   * By training models on **different bootstraps**, bagging **reduces variance of the individual models**.
   * It also **avoids overfitting** by **exposing the constituent models to different parts of the data set**.
   * By averaging/voting the predictions from all the sampled models, the **aggregated/ensemble model incorporates the strengths of the individual models and cancels out their errors**
  

<br>
    <div align="center">
        <img src="https://github.com/c-vandenberg/machine-learning-in-drug-discovery/assets/60201356/d1151bfc-1843-437d-bde8-57dbae2b1e59", alt="bootstrap-bagging" width="750"/>
      <p>
        <b>Fig 4</b> Bootstrap bagging schematic. <b><sup>11</sup></b>
      </p>
    </div>
<br>

2. **Boosting**:
* Boosting  **<sup>12</sup>** is an ensemble learning method that **combines a set of weak learners (classifers that are only slightly correlated to the true classification** to give a **strong learner (classifier that is well-correlated with the true classification** differs from bagging in the 

## References
**[1]** What is a decision tree? (2021) *IBM*. Available at: https://www.ibm.com/topics/decision-trees (Accessed: 04 July 2024).
**[2]** Jones, T. (2017) 'Models for machine learning - Decision Trees', *IBM developer*. Available at: https://developer.ibm.com/articles/cc-models-machine-learning/#decision-trees2 (Accessed: 02 July 2024). <br><br>
**[3]** Black, P.E. (2005) 'greedy algorithm', *Dictionary of Algorithms and Data Structures*, (accessed 4th July 2024). Available at: https://www.nist.gov/dads/HTML/greedyalgo.html
**[4]**Quinlan, J.R. (1986) ‘Induction of Decision Trees’, *Machine Learning*, 1(1), pp. 81–106. <br><br>
**[5]** Salzberg, S.L. (1994) ‘C4.5: Programs for Machine Learning by J. Ross Quinlan. Morgan Kaufmann Publishers, Inc., 1993’, *Machine Learning*, 16(3), pp. 235–240. <br><br>
**[6]** Breiman, L. (2017) 'Classification and regression trees'. Abingdon: Routledge. <br><br>
**[7]** Opitz, D. and Maclin, R. (1999) ‘Popular Ensemble Methods: An Empirical Study’, *Journal of Artificial Intelligence Research*, 11, pp. 169–198. <br><br>
**[8]** Polikar, R. (2006) ‘Ensemble based systems in decision making’, *IEEE Circuits and Systems Magazine*, 6(3), pp. 21–45. <br><br>
**[9]** Rokach, L. (2009) ‘Ensemble-based classifiers’, *Artificial Intelligence Review*, 33(1–2), pp. 1–39. <br><br>
**[10]** Breiman, L. (1996) ‘Bagging predictors’, *Machine Learning*, 24(2), pp. 123–140. <br><br>
**[11]** Awan, A.A. (2023) A guide to bagging in machine learning: Ensemble method to reduce variance and improve accuracy, *DataCamp*. Available at: https://www.datacamp.com/tutorial/what-bagging-in-machine-learning-a-guide-with-examples (Accessed: 05 July 2024). <br><br>
**[12]** Schapire, R.E. (1990) ‘The strength of weak learnability’, *Machine Learning*, 5(2), pp. 197–227. <br><br>
