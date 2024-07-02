# 2.1 Introduction to Machine Learning

## 2.1.1 How does Machine Learning Work?

The concept of the **algorithm** is at the core of computer programming. An algorithm is simply a **sequence of instructions** to be carried out to **transform an input to an output**. Many algorithms exist for a variety of tasks, however for some tasks, an algoirthm has not been developed, despite years of research. Some of these tasks humans can do effortlessly, such as recognising an individual from a photograph, or generating an image from a written description. **<sup>1</sup>**

Machine learning (ML) is a branch of AI and computer science that focuses on using data to develop algorithms for these types of tasks. It is therefore a **method of modelling data**, typically with **predictive functions**.

Generally, how ML works can be divided into **three main parts**: **<sup>2</sup>**
1. **A Decision Process**:
   * In general, ML algorithms are used to **train mapping functions** to **make a prediction or classification**.
   * These mapping functions will have **various parameters/weights**
   * Based on an **input data set**, which can be **labelled** or **unlabelled**, the **parameters of the mapping function will be adjusted** so that it produces an **estimate about a pattern in the data**
   * The result of this training is a **model** which can be thought of as a **representation of the patterns and relationships** that the mapping function has **learned from the data**
2. **An Error Function (Critic)**:
   * An error function **evaluates the prediction of the model**.
   * If there are **known examples**, an error function can make a **comparison to assess the accuracy of the model**.
3. **A Model Optimisation Process**:
   * After comparing the models prediction to the known examples, if there are any discrepancies, the **parameters/weights are adjusted** to **reduce the discrepancy between the known example and model estimate**.
   * This **"evaluate and optimise"** process is **repeated iteratively**, updating the parameters/weights autonomously until a **threshold of accruacy has been met**.

     <br>
      <div align="center">
        <img src="https://github.com/c-vandenberg/machine-learning-in-drug-discovery/assets/60201356/dcdc9431-2e9f-406d-8a94-7e89b6449d63", alt="typical-flowchart-of-a-machine-learning-algorithm"/>
        <p>
          <b>Fig 1</b> Typical flowchart of a machine learning algorithm. <b><sup>3</sup></b>
        </p>
      </div>
     <br>

## 2.1.2 Machine Learning Methods

### Supervised Learning
* **Supervised learning** is defined by its use of **labelled data sets** to train **mapping functions** to **classify data** or **predict outcomes accurately**.
* As input data is fed into the model, the **model adjusts the mapping function's parameters/weights** until it has been **fitted appropriately** to the data set.
* The idea is that once the mapping function is trained, it can then be used on **unseen data to meet some predictive performance**
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
 
     <br>
      <div align="center">
        <img src="https://github.com/c-vandenberg/machine-learning-in-drug-discovery/assets/60201356/c96f4b00-d155-4d53-a4f5-89c64d99c44d", 
         alt="supervised-learning"/>
        <p>
          <b>Fig 2</b> The two phases of building and testing a mapping function with supervised learning. <b><sup>4</sup></b>
        </p>
      </div>
     <br>
      
Generally, you **build and test a models mapping function** in **two phases** in supervised learning (**Fig 2**):
 1. **First Phase**:
    * In the first phase, you **segment the data set** into **two subsets**: **training data** and **test data**.
    * Both subsets contain a **test vector** (the inputs) and **one or more known desired output values (if not unsupervised learning)**.
    * You then **train the mapping function** with the training data set and **alter its parameters/weights** until it **meets some level of performance** (a metric for how accurately the mapping function maps the training data to the associated desired output, evaluated by an **error function/critic**).
    * In this phase, the algorithm will conduct the **Decision Process**, **Error Function (Critic)**, and **Model Optimisation Process** iterative cycle described in **2.1.1
 2. **Second Phase**:
    * In the next phase, you **test the trained algorithm/mapping function against the test data**.
    * The test data represents data that **has not been used for training** and provides a good measure for **how well the mapping function generalises to unseen data**
 
### Unsupervised Learning
* **Unsupervised learning** uses ML algorithms to **analyse and cluster unlabelled data sets** into **subsets** called **clusters**.
* These algorithms **discover hidden patterns** or **data groupings** without the need for human interventions, making them ideal for **very large data sets**.
* As such, unsupervised learning algorithms **have no input-output pairs** and consist **only of input data, without any labelled output**.
* Because **there is no labelled output**, there is **no error function/critic to evaluate its predictions**
* The goal therefore is to **build a mapping function** that **categorises the data by identifying hidden patterns, groupings, or structures within the data**
* Unsupervised learning is also used to **reduce the number of features ina  model** through the process of **dimensionality reduction** (e.g. using **principal component analysis (PCA)** and **singular value decomposition (SVD)**)
* Algorithms used in unsupervised learning include:
  1. **Neural networks**
  2. **K-means clustering**
  3. **Probabilistic clustering methods**
 
     <br>
      <div align="center">
        <img src="https://github.com/c-vandenberg/machine-learning-in-drug-discovery/assets/60201356/eb2dabda-783e-4f0b-bb91-c2178c8e3661", 
         alt="unsupervised-learning"/>
        <p>
          <b>Fig 3</b> The two phases of building and testing a mapping function with unsupervised learning. <b><sup>4</sup></b>
        </p>
      </div>
     <br>
     
In unsupervised learning, the mapping function **segments the data set into classes**. Each **input vector** becomes **part of a class**, but the **algorithm cannot apply labels to those classes**  (**Fig 3**)
 
### Semi-Supervised Learning
* **Semi-supervised learning** offers a **happy medium** between supervised and unsupervised learning
* During training, it uses a **smaller labeled data set** to **guide classification and feature extraction** from a **larger, unlabeled data set.**
* Semi-supervised learning can solve the problem of not having enough labeled data for a supervised learning algorithm (e.g. if it is too costly to label enough data). **<sup>2</sup>**

### Reinforcement Learning
* In **reinforcement learning**, the algorithm not only attempts to **map an input to an output** (as with supervised learning), but to **map a series of inputs to outputs with dependencies**
* Reinforcement learning exists in the context of **states in an environment** and the **actions possible at a given state**.
* During the learning process, the algorithm (referred to as an **agent** in the context of reinforcement learning), **randomly explores the state-action pairs** within some environment, to build a **state-action pair table**
* Over many iterations of the learning process, the agent will **learn state-action pair rewards** in order to **rank the best action for a given state, that will lead to a pre-defined goal state**

     <br>
      <div align="center">
        <img src="https://github.com/c-vandenberg/machine-learning-in-drug-discovery/assets/60201356/a5aa16cd-b3f1-4a32-8def-2a71c85f0b58", 
         alt="reinforcement_learning"/>
        <p>
          <b>Fig 4</b> The reinforcement learning model. <b><sup>4</sup></b>
        </p>
      </div>
     <br>
      
An illustrative example of reinforcement learning is **a simple agent that plays blackjack**:
 1. Here, the **states represent the sum of the cards for the player**.
 2. The **actions represent what a blackjack-playing agent may do** (in this case, **hit or stand**).
 3. Training this agent to play blackjack would involve many hands, where **reward for a given state–action nexus is given for winning or losing**.
 4. For example, the value for a **state of 10** would be 1.0 for hit and 0.0 for stand (indicating that hit is the optimal choice).
 5. For **state 20**, the learned reward would likely be 0.0 for hit and 1.0 for stand.
 6. For a less-straightforward hand, e.g. a **state of 17** the action values could be 0.95 stand and 0.05 hit. This agent would then **probabilistically stand 95 percent of the time** and **hit 5 percent of the time**.
 7. These rewards would be **leaned over many hands**, indicating the **best choice for a given state (or hand).**

### Supervised vs Unsupervised vs Reinforcement Learning

   <br>
    <div align="center">
      <img src="https://github.com/c-vandenberg/machine-learning-in-drug-discovery/assets/60201356/a71447aa-84b0-425f-ad63-773b186dbe63", 
       alt="supervised-vs-unsupervised-vs-reinforcement"/>
      <p>
        <b>Fig 4</b> The reinforcement learning model. <b><sup>4</sup></b>
      </p>
    </div>
   <br>

In summary:
1. In **supervised learning**:
   * The data set **includes its desired outputs (or labels)** such that an **error function/critic can calculate an error for a given prediction**.
   * The **supervision** comes when a **prediction is made and an error produced**, which in turn **alters the mapping function** to **learn the input-output mapping**
2. In **unsupervised learning**:
   * The data set **doesn't include a desired output**; therefore, there is **no way to supervise the mapping function with an error function/critic**.
   * Instead, the mapping function attempts to **segment the data set into "classes"** so that each class contains a portion of the data set **with common features**
3. In **reinforcement learning**:
   * The algorithm attempts to **learn actions for a given set of states that lead to a goal state**.
   * An error is provided **not after each example** (as is the case for supervised learning), but instead **on receiving a reinforcement signal**.
   * This learning is **similar to human learning**, where **feedback isn't necessarily provided for all actions** but when a **reward is warranted**.
   
## 2.1.3 Machine Learning Notation

Below are some **common definitions and notation** in machine learning: **<sup>5</sup>**

1. **Features**:
     * A set of ![N](https://latex.codecogs.com/svg.latex?\color{white}N) vectors ![feature_vector}](https://latex.codecogs.com/svg.latex?\color{white}\{\vec{x}_i\}) each having a dimension ![D](https://latex.codecogs.com/svg.latex?\color{white}D).
     * These vectors can consist of real numbers, integers, or other types of values
2. **Labels**:
     * A set of ![N](https://latex.codecogs.com/svg.latex?\color{white}N) integers or real numbers ![label_yi}](https://latex.codecogs.com/svg.latex?\color{white}\{y_i\}).​![label_y](https://latex.codecogs.com/svg.latex?\color{white}y_i).
     * Each ![label_yi}](https://latex.codecogs.com/svg.latex?\color{white}\{y_i\}).​![label_y](https://latex.codecogs.com/svg.latex?\color{white}y_i) is typically a scalar value representing the target output associated with a feature vector.
3. **Labeled Data**:
     * A set of ![N](https://latex.codecogs.com/svg.latex?\color{white}N) tuples ![labelled_data](https://latex.codecogs.com/svg.latex?\color{white}\{\left(\vec{x}_i,y_i\right)\}), where each tuple consists of a feature vector and its corresponding label.
4. **Unlabeled Data**:
     * A set of ![N](https://latex.codecogs.com/svg.latex?\color{white}N) feature vectors ![unlabelled_data](https://latex.codecogs.com/svg.latex?\color{white}\{\vec{x}_i\})
 that do not have associated labels.
6. **Data Generation Process**:
     * The unknown function ![function](https://latex.codecogs.com/svg.latex?\color{white}f(\vec{x})) that, for a given feature vector, returns a real-valued label (output) ![output](https://latex.codecogs.com/svg.latex?\color{white}{y}).
     * This is the process we aim to model with machine learning.
7. **Model**:
     * A function ![mapping_function](https://latex.codecogs.com/svg.latex?\color{white}f(\vec{x})) that takes a given feature vector and returns a predicted output ![predicted_output](https://latex.codecogs.com/svg.latex?\color{white}\hat{y}).
     * The goal of machine learning is to make this function **map the input to the predicted output as accurately as possible** and make it as close to the data generation process as possible.
     * This is therefore a **mapping function**
8. **Predictions**:
     * The predicted output ![output](https://latex.codecogs.com/svg.latex?\color{white}\hat{y}) for a given input ![input](https://latex.codecogs.com/svg.latex?\color{white}\vec{x}), produced by the model ![mapping_function](https://latex.codecogs.com/svg.latex?\color{white}f(\vec{x})).

## 2.1.4 References
**[1]** Alpaydin, E. (2020) *Introduction to machine learning*. Cambridge, MA: The MIT Press. <br><br>
**[2]** *What is machine learning (ML)?* (2021) IBM. Available at: https://www.ibm.com/topics/machine-learning (Accessed: 01 July 2024). <br><br>
**[3]** Salah, S., Alsamamra, H.R. and Shoqeir, J.H. (2022) ‘Exploring wind speed for energy considerations in eastern Jerusalem-Palestine using machine-learning algorithms’, *Energies*, 15(7), p. 2602. <br><br>
**[4]** Jones, T. (2017) *Models for machine learning*, IBM developer. Available at: https://developer.ibm.com/articles/cc-models-machine-learning (Accessed: 02 July 2024). <br><br>
**[5]** White, A.D. (2022) ‘Deep learning for molecules and materials’, *Living Journal of Computational Molecular Science*, 3(1).
