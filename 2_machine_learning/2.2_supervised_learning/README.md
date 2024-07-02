# 2.2 Supervised Learning

## 2.2.1 Introduction to Supervised Learning

* **Supervised learning** is defined by its use of **labelled data sets** to train **mapping functions** to **classify data** or **predict outcomes accurately**.
* As input data is fed into the model, the **model adjusts the mapping function's parameters/weights** until it has been **fitted appropriately** to the data set.
* The idea is that once the mapping function is trained, it can then be used on **unseen data to meet some predictive performance**
* As such, supervised learning algorithms consist of **input-output pairs** and the goal is to **learn a mapping from inputs to outputs**.
* This adjusting of the parameters/weights occurs as part of the **"evaluate and optimise" process** to ensure that the model avoids **overfitting** or **underfitting**.

    <br>
      <div align="center">
        <img src="https://github.com/c-vandenberg/machine-learning-in-drug-discovery/assets/60201356/c96f4b00-d155-4d53-a4f5-89c64d99c44d", 
         alt="supervised-learning"/>
        <p>
          <b>Fig 1</b> The two phases of building and testing a mapping function with supervised learning. <b><sup>4</sup></b>
        </p>
      </div>
     <br>

Generally, how supervised machine learning works can be divided into **five main parts**: **<sup>2</sup>**
1. **Data Segmentation**:
      * Before any training begins, the data set is **segmented into two subsets**; **training data** and **test data**.
      * oth subsets contain a **test vector** (the inputs) and **one or more known desired output values**.
2. **A Decision Process**:
      * In general, ML algorithms are used to **train mapping functions** to **make a prediction or classification**.
      * These mapping functions will have **various parameters/weights**
      * Based on an **input data set**, which can be **labelled** or **unlabelled**, the **parameters of the mapping function will be adjusted** so that it produces an **estimate about a pattern in the data**
  * The result of this training is a **model** which can be thought of as a **representation of the patterns and relationships** that the mapping function has **learned from the data**
3. **An Error Function (Critic)**:
  * An error function **evaluates the prediction of the model**.
  * If there are **known examples**, an error function can make a **comparison to assess the accuracy of the model**.
4. **A Model Optimisation Process**:
  * After comparing the models prediction to the known examples, if there are any discrepancies, the **parameters/weights are adjusted** to **reduce the discrepancy between the known example and model estimate**.
  * This **"evaluate and optimise"** process (steps 2-4) is **repeated iteratively**, updating the parameters/weights autonomously until a **threshold of accruacy has been met**.
5. **Test Trained Model Against Test Data**
  * Finally, the trained model/mapping function is **tested against the test data**
  * The test data represents data that **has not been used for training** and provides a good measure for **how well the mapping function generalises to unseen data**

The supervised learning algorithms we will discuss include:
1. **Regression**
2. **Random Forest**
3. **Classficiation**
