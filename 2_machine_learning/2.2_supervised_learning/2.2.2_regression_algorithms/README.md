# 2.2.2 Regression Algorithms in Supervised Learning

## 1) Linear Regression

**Linear regression** is one of the simplest supervised machine learning algorithms. Linear regression is used to find the **linear relationship between the dependent variable (output $`y`$)** and **one or more independent variables (input vectors {$`\vec{x}_i`$}**). It is **typically leveraged to make predictions about future outcomes**.

Below are the different types of linear regression. As with all other types of regression, all linear regression types **seek to plot a line of best fit**, which is calculated through the method of **least squares**. What makes linear regression differ from other types of regression (e.g. **polynomial regression**) is that this line of best fit is a **straight line when plotted on a graph**

### Simple Linear Regression 

**Simple linear regression** is the **simplest form of linear regression** and involves **only one independent feature/variable**, and **one dependent variable/output**.

The equation for simple linear regression is:
<br>
<br>

$$y = \beta_0 + \beta_1 X$$

where:
- $`Y`$ is the **dependent variable/output**
- $`X`$ is the **independent variable/input vector**
- $`\beta_0`$ is the **intercept**
- $`\beta_1`$ is the **slope/gradient**


## Multiple Linear Regression

**Multiple linear regression** involves **more than one independent feature/variable**, and **one dependent variable/output**.

The equation for multiple linear regression (and also **univariate linear regression**) is:
<br>
<br>

$$y = \beta_0 + \beta_1 X_1 + \beta_2 X_2 + \ldots + \beta_n X_n$$

where:
* $`Y`$ is the **dependent variable/output**
* $`\beta_0`$ is the **intercept**
* $`\beta_1, \beta_2, \ldots, \beta_n`$ are the **slopes/gradients**
* $`{X_1, X_2, \ldots, X_p}`$ are the **independent variables/input vectors**

Multiple linear regression can be represented in **matrix form**, which provides a compact and efficient way to represent and compute the **relationships between multiple dependent and independent variables**:
<br>
<br>

$$\mathbf{y}=\mathbf{X}\boldsymbol{\beta}$$

where:

$`Y=\begin{bmatrix}Y_{1}\\Y_{2}\\\vdots\\Y_{n}\end{bmatrix}`$

<br>

$`X=\begin{bmatrix}x_{11}&x_{12}&\cdots&x_{1n}\\x_{21}&x_{22}&\cdots&x_{2n}\\\vdots&\vdots&\ddots&\vdots\\x_{n1}&x_{n2}&\cdots&x_{nn}\end{bmatrix}`$

<br>

$`\beta=\begin{bmatrix}\beta_{0}\\\beta_{1}\\\beta_{2}\\\vdots\\\beta_{n}\end{bmatrix}`$


## Multvariate Linear Regression

**Multivariate linear regression** involves **more than one dependent variables/outputs**.

The equation for multivariate linear regression (and also **univariate linear regression**) is:
<br>
<br>

$$\color{white}Y_i = \beta_{0} + \beta_{1}x_{i}^{(1)} + \beta_{2}x_{i}^{(2)} + \ldots + \beta_{n}x_{i}^{(n)}$$

where:
* $`Y_i`$ is the **dependent variable/output for the *i*-th observation**
* $`\beta_{0}`$ is the **intercept**
* $`\beta_{1}x_{i}^{(1)}`$ is the **coefficient and the *i*-th observation of the first independent variable**
* $`\beta_{2}x_{i}^{(2)}`$ is the **coefficient and the *i*-th observation of the second independent variable**
* $`\ldots`$ represents **additional terms up to *n* independent variables**
* $`\beta_{n}x_{i}^{(n)}`$ is the **coefficient and the *i*-th observation of the *n*-th independent variable**


Like with multiple linear regression, multivariate linear regression can be represented in **matrix form**:
<br>
<br>

$$\mathbf{y}=\mathbf{X}\boldsymbol{\beta}$$

where:


$`Y=\begin{bmatrix}Y_{1}\\Y_{2}\\\vdots\\Y_{n}\end{bmatrix}`$

<br>

$`X=\begin{bmatrix}x_{11}&x_{12}&\cdots&x_{1n}\\x_{21}&x_{22}&\cdots&x_{2n}\\\vdots&\vdots&\ddots&\vdots\\x_{n1}&x_{n2}&\cdots&x_{nn}\end{bmatrix}`$

<br>

$`\beta=\begin{bmatrix}\beta_{0}\\\beta_{1}\\\beta_{2}\\\vdots\\\beta_{n}\end{bmatrix}`$

## 2) Logistic Regression

**Logisitic Regression** is a supervised machine learning algorithm that **accomplishes binary classification tasks** by **predicting the probability of an outcome, event, or observation**. While linear regression is **leveraged when dependent variables are continuous**, logistic regression is selected when the **dependent variables are categorical**, meaning they have **binary outputs** (e.g. 'true' and 'false', 'yes' and 'no' etc.).

While both linear and logistic regression seek to **understand relationships between data inputs** logistic regression is mainly used to solve **binary classification problems** (e.g. spam identification).

Logisitic regression utilises what is known as the **logistic function**, which produces an **S-shaped curve**, also known as a **sigmoid curve**. The logistic function has the general equation:

$$f(x)=\frac{L}{1+e^{-k(x-x_{0})}}$$

where:
* $`L`$ - is the **supremum of the values of the function**. The supremum refers to the **upper bound for the set of values of the function** (i.e. it is the **maximum value that the function can attain**
* $`k`$ - is the **logistic growth rate** (i.e. the **steepness/gradient of the curve**)
* $`x_{0}`$ - is the **x value of the function's midpoint**


The **standard logistic function**, where $`L`$ = 1, $`k`$ = 1, and $`x_{0}`$ = 0 has the below equation and is plotted in **Fig 2**:
<br>
<br>

$$f(x)=\frac{1}{1+e^{-x}}$$

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
    * **$`f({x})`$ tends towards 1 as $x \to \infty$**
    * **$`f({x})`$ tends towards 0 as $x \to -\infty$**
    * **$`f({x})`$ is always bounded between 0 and 1**

### 3) Polynomial Regression

**Polynomial regression** is a type of regression analysis employed in machine learning when the relationship between the independent variables/input vectors and the dependent variable(s)/output(s) is **not linear**, but is instead an **nth-degree polynomial**:
<br>
<br>

$$y=\beta_0+\beta_1x+\beta_2x^2+\cdots+\beta_nx^n$$

where:
* $`Y`$ is the **dependent variable/input vectors**
* $`\beta_0`$ is the **intercept**
* $`\beta_1,\ \beta_2,\ \ldots,\ \beta_n`$ are the **1st, 2nd to nth independent variable/input vector coefficients**
* $`X`$ are the **independent variable/input vectors**

As stated before, by **including higher degree terms (squared, cubic, quadratic etc.)**, a polynomial regression model can **capture non-linear patterns in the data set**:
1. **The choice of the polynomial degree (n)** is a **crucial aspect** of polynomial regression. A **higher degree** allows the model to **fit the training data more closely**, but it may also **lead to overfitting**, especially if the degree is too high.
2. The polynomial regression model is trained to **find the independent variable/input vector coefficients** that **minimize the difference between the predicted values and the actual values** in the training data.
3. Once the model is trained, it can be used to **make predictions on new, unseen data**. The polynomial equation captures the **non-linear patterns** observed in the training data, allowing the model to **generalize to non-linear relationships in unseen data**.

## References
 **[1]** Vadurova, D., Gausland, Yvonne. (2018). "Analysis of consumers' ability to identify fake news". <br><br>
