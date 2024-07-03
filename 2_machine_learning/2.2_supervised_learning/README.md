# 2.2 Supervised Learning

## 2.2.1 Introduction to Supervised Learning

* **Supervised learning** is defined by its use of **labelled data sets** to train **mapping functions** to **classify data** or **predict outcomes accurately**.
* In other words, supervised learning invovles predicting **the output ![output](https://latex.codecogs.com/svg.latex?\color{white}{y})** for a **given input vector ![input_vectors](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%7B%5C%7B%5Cvec%7Bx%7D_i%5C%7D%7D)** 
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

Supervised learning can be **separated into two distinct types**
1. **Classficiation**:
   * **Classification** is a supervised learning method whereby the algorithm is used to **accurately assign test data into specific categories**.
   * The algorithm is trained to **recognise specific entities within the data set** and attempts to **draw some conclusions on how those entities should be labelled or defined**
   * Common classification algorithms include:
       1. **Linear Classifiers**
       2. **Support Vector Machines (SVM)**
       3. **Decision Trees**
       4. **K-Nearest Neighbour**
       5. **Random Forest**
2. **Regression**:
   * **Regression** is a supervised learning method used to **understand the relationship between dependent and independent variables**.
   * In other words, regression uses mathematical methods to **predict a continuous outcome ![output](https://latex.codecogs.com/svg.latex?\color{white}{y})** based on the value of **one or more predictor variables/input vectors ![input_vectors](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%7B%5C%7B%5Cvec%7Bx%7D_i%5C%7D%7D)**
   * In regression, the algorithm is **trained on both input vectors and output labels**. This helps in **establishing a relationship among the variables** by **estimating how one variable affects the other**.
   * Common regression algorithms include:
       1. **Linear Regression**
       2. **Logistic Regression**
       3. **Polynomial Regression**

## 2.2.2 Classification Algorithms in Supervised Learning

## 2.2.3 Regression Algorithms in Supervised Learning

### 1) Linear Regression

**Linear regression** is one of the simplest supervised machine learning algorithms. Linear regression is used to find the **linear relationship between the dependent variable (output ![output](https://latex.codecogs.com/svg.latex?\color{white}{y})) and one or more independent variables (input vectors ![input_vectors](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%7B%5C%7B%5Cvec%7Bx%7D_i%5C%7D%7D)**). It is **typically leveraged to make predictions about future outcomes**.

Below are the different types of linear regression. As with all other types of regression, all linear regression types **seek to plot a line of best fit**, which is calculated through the method of **least squares**. What makes linear regression differ from other types of regression (e.g. **polynomial regression**) is that this line of best fit is a **straight line when plotted on a graph**

### Simple Linear Regression 

**Simple linear regression** is the **simplest form of linear regression** and involves **only one independent feature/variable**, and **one dependent variable/output**.

The equation for simple linear regression is:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7Dy%20%3D%20%5Cbeta_0%20%2B%20%5Cbeta_1%20X", alt='simple_linear_regression'/>
    </div>
<br>

where:
- ![Y](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DY) is the dependent variable
- ![X](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DX) is the independent variable
- ![\beta_0](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_0) is the intercept
- ![\beta_1](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_1) is the slope


### Multiple Linear Regression

**Multiple linear regression** involves **more than one independent feature/variable**, and **one dependent variable/output**.

The equation for multiple linear regression (and also **univariate linear regression**) is:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7Dy%20%3D%20%5Cbeta_0%20%2B%20%5Cbeta_1%20X_1%20%2B%20%5Cbeta_2%20X_2%20%2B%20%5Cldots%20%2B%20%5Cbeta_n%20X_n", 
            alt='multiple_linear_regression'/>
    </div>
<br>

where:
* ![Y](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DY) is the dependent variable
* ![X_1, X_2, \ldots, X_p](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DX_1%2C%20X_2%2C%20%5Cldots%2C%20X_p) are the independent variables
* ![\beta_0](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_0) is the intercept
* ![\beta_1, \beta_2, \ldots, \beta_n](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_1%2C%20%5Cbeta_2%2C%20%5Cldots%2C%20%5Cbeta_n) are the slopes

Multiple linear regression can be represented in **matrix form**, which provides a compact and efficient way to represent and compute the **relationships between multiple dependent and independent variables**:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?\color{white}\mathbf{y}=\mathbf{X}\boldsymbol{\beta}", alt='multiple_linear_regression_matrix'/>
    </div>
<br>

where:

<br>
    <div>
        <img src="https://latex.codecogs.com/svg.latex?\color{white}Y=\begin{bmatrix}Y_{1}\\Y_{2}\\\vdots\\Y_{n}\end{bmatrix}", alt='output_vector_matrix'/>
    </div>
<br>

<br>
    <div>
        <img src="https://latex.codecogs.com/svg.latex?\color{white}X=\begin{bmatrix}x_{11}&x_{12}&\cdots&x_{1n}\\x_{21}&x_{22}&\cdots&x_{2n}\\\vdots&\vdots&\ddots&\vdots\\x_{n1}&x_{n2}&\cdots&x_{nn}\end{bmatrix}", alt='input_vector_matrix'/>
    </div>
<br>

<br>
    <div>
        <img src="https://latex.codecogs.com/svg.latex?\color{white}\beta=\begin{bmatrix}\beta_{0}\\\beta_{1}\\\beta_{2}\\\vdots\\\beta_{n}\end{bmatrix}", alt='beta_vector_matrix'/>
    </div>
<br>


### Multvariate Linear Regression

**Multivariate linear regression** involves **more than one dependent variables/outputs**

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DY_i%20%3D%20%5Cbeta_%7B0%7D%20%2B%20%5Cbeta_%7B1%7Dx_%7Bi%7D%5E%7B%281%29%7D%20%2B%20%5Cbeta_%7B2%7Dx_%7Bi%7D%5E%7B%282%29%7D%20%2B%20%5Cldots%20%2B%20%5Cbeta_%7Bn%7Dx_%7Bi%7D%5E%7B%28n%29%7D", alt='multivariate_linear_regression'/>
    </div>
<br>

where:
* ![Y_i](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DY_i) is the dependent variable for the \(i\)-th observation
* ![\alpha](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_%7B0%7D) is the intercept
* ![\beta_{1}x_{i}^{(1)}](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_%7B1%7Dx_%7Bi%7D%5E%7B%281%29%7D) is the coefficient and the \(i\)-th observation of the first independent variable
* ![\beta_{2}x_{i}^{(2)}](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_%7B2%7Dx_%7Bi%7D%5E%7B%282%29%7D) is the coefficient and the \(i\)-th observation of the second independent variable
* ![\ldots](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cldots)  represents additional terms up to \(n\) independent variables
* ![\beta_{n}x_{i}^{(n)}](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_%7Bn%7Dx_%7Bi%7D%5E%7B%28n%29%7D) is the coefficient and the \(i\)-th observation of the \(n\)-th independent variable


Like with multiple linear regression, multivariate linear regression can be represented in **matrix form**:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?\color{white}\mathbf{y}=\mathbf{X}\boldsymbol{\beta}", alt='multiple_linear_regression_matrix'/>
    </div>
<br>

where:

<br>
    <div>
        <img src="https://latex.codecogs.com/svg.latex?\color{white}Y=\begin{bmatrix}Y_{1}\\Y_{2}\\\vdots\\Y_{n}\end{bmatrix}", alt='output_vector_matrix'/>
    </div>
<br>

<br>
    <div>
        <img src="https://latex.codecogs.com/svg.latex?\color{white}X=\begin{bmatrix}x_{11}&x_{12}&\cdots&x_{1n}\\x_{21}&x_{22}&\cdots&x_{2n}\\\vdots&\vdots&\ddots&\vdots\\x_{n1}&x_{n2}&\cdots&x_{nn}\end{bmatrix}", alt='input_vector_matrix'/>
    </div>
<br>

<br>
    <div>
        <img src="https://latex.codecogs.com/svg.latex?\color{white}\beta=\begin{bmatrix}\beta_{0}\\\beta_{1}\\\beta_{2}\\\vdots\\\beta_{n}\end{bmatrix}", alt='beta_vector_matrix'/>
    </div>
<br>

### 2) Logistic Regression

**Logisitic Regression** is a supervised machine learning algorithm that **accomplishes binary classification tasks** by **predicting the probability of an outcome, event, or observation**. While linear regression is **leveraged when dependent variables are continuous**, logistic regression is selected when the **dependent variables are categorical**, meaning they have **binary outputs** (e.g. 'true' and 'false', 'yes' and 'no' etc.).

While both linear and logistic regression seek to **understand relationships between data inputs** logistic regression is mainly used to solve **binary classification problems** (e.g. spam identification).

Logisitic regression utilises what is known as the **logistic function**, which produces an **S-shaped curve**, also known as a **sigmoid curve**. The logistic function has the general equation:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?\color{white}f(x)=\frac{L}{1+e^{-k(x-x_{0})}}", alt='logistic_function_equation'/>
    </div>
<br>

where:
* ![supremum](https://latex.codecogs.com/svg.latex?\color{white}L) - is the **supremum of the values of the function**. The supremum refers to the **upper bound for the set of values of the function** (i.e. it is the **maximum value that the function can attain**
* ![logistic_growth_rate](https://latex.codecogs.com/svg.latex?\color{white}k) -  is the **logistic growth rate** (i.e. the **steepness/gradient of the curve**)
* ![function_midpoint_x_value](https://latex.codecogs.com/svg.latex?\color{white}x_{0}) - is the **x value of the function's midpoint**


The **standard logistic function**, where ![supremum](https://latex.codecogs.com/svg.latex?\color{white}L) = 1, ![logistic_growth_rate](https://latex.codecogs.com/svg.latex?\color{white}k) = 1, and ![function_midpoint_x_value](https://latex.codecogs.com/svg.latex?\color{white}x_{0}) = 0 has the below equation and is plotted in **Fig 2**:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?\color{white}f(x)=\frac{1}{1+e^{-x}}", alt='standard_logistic_function_equation'/>
    </div>
<br>

<br>
  <div align="center">
    <img src="", 
     alt=""/>
    <p>
      <b>Fig 2</b> . <b><sup>4</sup></b>
    </p>
  </div>
 <br>
