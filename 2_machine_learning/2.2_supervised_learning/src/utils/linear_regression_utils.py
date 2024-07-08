import jax.numpy as jnp
import numpy as np


def linear_regression(feature_vector: np.array, weight_vector: np.array, bias_term: float) -> np.array:
    """
    Calculates predicted output vector `y` using the linear regression equation.

    Parameters
    ----------
    feature_vector : np.array
        The input feature matrix where each row represents a data sample and each column represents a feature.
    weight_vector : np.array
        The weight vector (parameters) for the linear regression model.
    bias_term : float
        The bias term (intercept) for the linear regression model.

    Returns
    -------
    np.array
        The predicted output vector `y` obtained by applying the linear regression model.
    """
    return jnp.dot(feature_vector, weight_vector) + bias_term


def mean_squared_error_loss(model_prediction: np.array, labels: np.array) -> np.array:
    """
    A mean squared error (MSE) loss function. It takes the predicted output vector `y` and the input vector labels (true
    output values) and returns the average of the squared differences between them.

    Parameters
    ----------
    model_prediction : np.array
        The predicted output vector `y` obtained from the linear regression model.
    labels : np.array
        The true output values (labels) corresponding to the input feature vectors.

    Returns
    -------
    float
        The mean squared error (MSE) calculated as the average of the squared differences between the predicted and
        true values.
    """
    return jnp.mean((model_prediction - labels) ** 2)
