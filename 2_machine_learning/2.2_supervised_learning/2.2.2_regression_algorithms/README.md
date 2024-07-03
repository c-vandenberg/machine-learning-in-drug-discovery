# 2.2.3 Regression Algorithms in Supervised Learning

## 1) Linear Regression

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
- ![simple_linear_regression_dependent_var](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DY) is the **dependent variable/output**
- ![simple_linear_regression_independent_var](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DX) is the **independent variable/input vector**
- ![simple_linear_regression_intercept](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_0) is the **intercept**
- ![simple_linear_regression_gradient](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_1) is the **slope/gradient**


## Multiple Linear Regression

**Multiple linear regression** involves **more than one independent feature/variable**, and **one dependent variable/output**.

The equation for multiple linear regression (and also **univariate linear regression**) is:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7Dy%20%3D%20%5Cbeta_0%20%2B%20%5Cbeta_1%20X_1%20%2B%20%5Cbeta_2%20X_2%20%2B%20%5Cldots%20%2B%20%5Cbeta_n%20X_n", 
            alt='multiple_linear_regression'/>
    </div>
<br>

where:
* ![mulitple_linear_regression_dependent_var](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DY) is the **dependent variable/output**
* ![mulitple_linear_regression_independent_vars](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DX_1%2C%20X_2%2C%20%5Cldots%2C%20X_p) are the **independent variables/input vectors**
* ![mulitple_linear_regression__intercept](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_0) is the **intercept**
* ![mulitple_linear_regression_gradient](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_1%2C%20%5Cbeta_2%2C%20%5Cldots%2C%20%5Cbeta_n) is the **slope/gradient**

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


## Multvariate Linear Regression

**Multivariate linear regression** involves **more than one dependent variables/outputs**.

The equation for multivariate linear regression (and also **univariate linear regression**) is:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DY_i%20%3D%20%5Cbeta_%7B0%7D%20%2B%20%5Cbeta_%7B1%7Dx_%7Bi%7D%5E%7B%281%29%7D%20%2B%20%5Cbeta_%7B2%7Dx_%7Bi%7D%5E%7B%282%29%7D%20%2B%20%5Cldots%20%2B%20%5Cbeta_%7Bn%7Dx_%7Bi%7D%5E%7B%28n%29%7D", alt='multivariate_linear_regression'/>
    </div>
<br>

where:
* ![multivariate_linear_regression_independent_ith_var](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7DY_i) is the **dependent variable/output for the *i*-th observation**
* ![multivariate_linear_regression_intercept](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_%7B0%7D) is the **intercept**
* ![multivariate_linear_regression_1st_independent_var_coeff](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_%7B1%7Dx_%7Bi%7D%5E%7B%281%29%7D) is the **coefficient and the *i*-th observation of the first independent variable**
* ![multivariate_linear_regression_2nd_independent_var_coeff](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_%7B2%7Dx_%7Bi%7D%5E%7B%282%29%7D) is the **coefficient and the *i*-th observation of the second independent variable**
* ![multivariate_linear_regression_n_independent_vars](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cldots)  represents **additional terms up to *n* independent variables**
* ![multivariate_linear_regression_final_independent_var](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_%7Bn%7Dx_%7Bi%7D%5E%7B%28n%29%7D) is the **coefficient and the *i*-th observation of the *n*-th independent variable**


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

## 2) Logistic Regression

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
    <img src="https://github.com/c-vandenberg/machine-learning-in-drug-discovery/assets/60201356/ef0af9b3-c922-4160-89e3-d14949569d01", alt="standard-logistic-function-curve"/>
    <p>
      <b>Fig 2</b> The standard logistic function sigmoid curve. <b><sup>1</sup></b>
    </p>
  </div>
 <br>

At a high level, how logistic regression works is:
1. The logistic function **transforms the linear regression function continous value output into a categorical value output**. 
2. It therefore **maps the set of independent variables/input vectors** into a **predicted dependent variable/output value between 0 and 1** (i.e. a **probability between 0 and 1**).
3. We can see from **Fig 2** that:
    * **![function](https://latex.codecogs.com/svg.latex?\color{white}f({x})) tends towards 1 as ![x_to_inf](https://latex.codecogs.com/svg.latex?\color{white}x\to\infty)**
    * **![function](https://latex.codecogs.com/svg.latex?\color{white}f({x})) tends towards 0 as ![x_to_inf](https://latex.codecogs.com/svg.latex?\color{white}x\to-\infty)**
    * **![function](https://latex.codecogs.com/svg.latex?\color{white}f({x})) is always bounded between 0 and 1**

### 3) Polynomial Regression

**Polynomial regression** is a type of regression analysis employed in machine learning when the relationship between the independent variables/input vectors and the dependent variable(s)/output(s) is **not linear**, but is instead an **nth-degree polynomial**:

<br>
    <div align="center">
        <img src="https://latex.codecogs.com/svg.latex?\color{white}y=\beta_0+\beta_1x+\beta_2x^2+\cdots+\beta_nx^n", alt='polynomial_regression_function_equation'/>
    </div>
<br>

where:
* ![polynomial_regression_independent_var](https://latex.codecogs.com/svg.latex?\color{white}x) is the **independent variable/input vectors**
* ![polynomial_regression_dependent_var](https://latex.codecogs.com/svg.latex?\color{white}y) is the **dependent variable/input vectors**
* ![polynomial_regression_intercept](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cbeta_0) is the **intercept**
* ![polynomial_regression_coeffs](https://latex.codecogs.com/svg.latex?\color{white}\beta_1,\\beta_2,\\ldots,\\beta_n) are the **1st, 2nd to nth independent variable/input vector coefficients**

As stated before, by **including higher degree terms (squared, cubic, quadratic etc.)**, a polynomial regression model can **capture non-linear patterns in the data set**:
1. **The choice of the polynomial degree (n)** is a **crucial aspect** of polynomial regression. A **higher degree** allows the model to **fit the training data more closely**, but it may also **lead to overfitting**, especially if the degree is too high.
2. The polynomial regression model is trained to **find the independent variable/input vector coefficients** that **minimize the difference between the predicted values and the actual values** in the training data.
3. Once the model is trained, it can be used to **make predictions on new, unseen data**. The polynomial equation captures the **non-linear patterns** observed in the training data, allowing the model to **generalize to non-linear relationships in unseen data**.

 ## References
 **[1]** Vadurova, D., Gausland, Yvonne. (2018). "Analysis of consumers' ability to identify fake news". <br><br>
