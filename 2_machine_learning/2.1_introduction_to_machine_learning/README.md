# 2.1 Introduction to Machine Learning

## 2.1.1 How does Machine Learning Work?

The concept of the **algorithm** is at the core of computer programming. An algorithm is simply a **sequence of instructions** to be carried out to **transform an input to an output**. Many algorithms exist for a variety of tasks, however for some tasks, an algoirthm has not been developed, despite years of research. Some of these tasks humans can do effortlessly, such as recognising an individual from a photograph, or generating an image from a written description. **<sup>1</sup>**

Machine learning (ML) is a branch of AI and computer science that focuses on using data to develop algorithms for these types of tasks. It is therefore a **method of modelling data**, typically with **predictive functions**.

Generally, how ML works can be divided into **three main parts**: **<sup>2</sup>**
1. **A Decision Process**:
  * In general, ML algorithms are used to **make a prediction or classification**.
  * These algorithms will have **various parameters/weights**
  * Based on an **input data set**, which can be **labelled** or **unlabelled**, the **parameters of the algorithm will be adjusted** so that it produces an **estimate about a pattern in the data**
  * The result of this algorithm training is a **model** which can be thought of as a **representation of the patterns and relationships** that the algorithm has **learned from the data**
2. **An Error Function**:
  * An error function **evaluates the prediction of the model**.
  * If there are **known examples**, an error function can make a **comparison to assess the accuracy of the model**.
3. **A Model Optimisation Process**:
  * After comparing the models prediction to the known examples, if there are any discrepancies, the **parameters/weights are adjusted** to **reduce the discrepancy between the known example and model estimate**.
  * This **"evaluate and optimise"** process is **repeated iteratively**, updating the parameters/weights autonomously until a **threshold of accruacy has been met**.

      <br>
      <div align="center">
        <img src="https://github.com/c-vandenberg/machine-learning-in-drug-discovery/assets/60201356/dcdc9431-2e9f-406d-8a94-7e89b6449d63", alt="typical-flowchart-of-a-machine-learning-algorithm"/>
        <p>
          <b>Fig 1</b> Typical flowchart of a machine learning algorithm
        </p>
      </div>
      <br>

## 2.1.2 Machine Learning Methods

### Supervised Learning
* **Supervised learning** is defined by its use of **labelled data sets** to train ML algorithms to **classify data** or **predict outcomes accurately**.
* As input data is fed into the model, the **model adjusts its parameters/weights** until it has been **fitted appropriately** to the data set.
* As such, supervised learning algorithms consist of **input-output pairs** and the goal is to **learn a mapping from inputs to outputs**.
* This adjusting of the parameters/weights occurs as part of the **"evaluate and optimise" process** to ensure that the model avoids **overfitting** or **underfitting**.
* **Overfitting** occurs when the model **fits too closely** or even **fits exactly** to its training data set. This results in a model that **can't make accurate predictions or conclusions from any data other than the training data**.
* **Underfitting** is a scenario where the model is **unable to capture the relationship between the input and output variables accurately**, generating a **high error rate** on both the training data set, and unseen data set.
* Algorithms used in supervised learning include:
  1. **Neural networks**
  2. **Naïve bayes**
  3. **Linear regression**
  4. **Logistic regression**
  5. **Random forest**
  6. **Support vector machine (SVM)**
 
### Unsupervised Learning
* **Unsupervised learning** uses ML algorithms to **analyse and cluster unlabelled data sets** into **subsets** called **clusters**.
* These algorithms **discover hidden patterns** or **data groupings** without the need for human interventions, making them ideal for **very large data sets**.
* As such, unsupervised learning algorithms **have no input-output pairs** and consist **only of input data, without any labelled output**. The goal therefore is to **identify patterns, groupings, or structures within the data**
* Unsupervised learning is also used to **reduce the number of features ina  model** through the process of **dimensionality reduction** (e.g. using **principal component analysis (PCA)** and **singular value decomposition (SVD)**)
* Algorithms used in unsupervised learning include:
  1. **Neural networks**
  2. **K-means clustering**
  3. **Probabilistic clustering methods**
 
### Semi-Supervised Learning
* **Semi-supervised learning** offers a **happy medium** between supervised and unsupervised learning
* During training, it uses a **smaller labeled data set** to **guide classification and feature extraction** from a **larger, unlabeled data set.**
* Semi-supervised learning can solve the problem of not having enough labeled data for a supervised learning algorithm (e.g. if it is too costly to label enough data). **<sup>2</sup>**

### Reinforcement Learning
* In **reinforcement learning**, the algorithm not only attempts to **map an input to an output** (as with supervised learning), but to **map a series of inputs to outputs with dependencies**
* Reinforcement learning exists in the context of **states in an environment** and the **actions possible at a given state**.
* During the learning process, the algorithm **randomly explores the state-action pairs** within some environment, to build a **state-action pair table**
* Then, 

##  2.X References
**[1]** Alpaydin, E. (2020) *Introduction to machine learning*. Cambridge, MA: The MIT Press. <br><br>
**[2]** What is machine learning (ML)? (2021) IBM. Available at: https://www.ibm.com/topics/machine-learning (Accessed: 01 July 2024). <br><br>
**[3]** Salah, S., Alsamamra, H.R. and Shoqeir, J.H. (2022) ‘Exploring wind speed for energy considerations in eastern Jerusalem-Palestine using machine-learning algorithms’, *Energies*, 15(7), p. 2602. <br><br>
