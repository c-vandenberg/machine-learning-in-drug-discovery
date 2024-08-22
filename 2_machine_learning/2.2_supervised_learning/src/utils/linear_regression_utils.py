import jax.numpy as jnp
import numpy as np
from typing import Union
from jaxlib.xla_extension import ArrayImpl
from pandas import DataFrame


def linear_regression(feature_vector: np.array, weight_vector: Union[np.ndarray, ArrayImpl],
                      bias_term: Union[float, ArrayImpl]) -> ArrayImpl:
    """
    Calculates predicted output vector `y` using the linear regression equation.

    Parameters
    ----------
    feature_vector : np.array
        The input feature matrix where each row represents a data sample and each column represents a feature.
    weight_vector : Union[np.ndarray, ArrayImpl]
        The weight vector (parameters) for the linear regression model.
    bias_term : Union[float, ArrayImpl]
        The bias term (intercept) for the linear regression model.

    Returns
    -------
    ArrayImpl
        The predicted output vector `y` obtained by applying the linear regression model.
    """
    return jnp.dot(feature_vector, weight_vector) + bias_term


def mean_squared_error_loss_wrapper(weight_vector: Union[np.ndarray, ArrayImpl], bias_term: Union[float, ArrayImpl],
                                    feature_values: np.ndarray, label_values: np.ndarray) -> ArrayImpl:
    """
    Calculates the mean squared error (MSE) loss for a linear regression model.

    Parameters
    ----------
    weight_vector : Union[np.ndarray, ArrayImpl]
        The weight vector (parameters) for the linear regression model.
    bias_term : Union[float, ArrayImpl]
        The bias term (intercept) for the linear regression model.
    feature_values : np.array
        The input feature matrix where each row represents a data sample and each column represents a feature.
    label_values : np.array
        The true output values (labels) corresponding to the input feature vectors.

    Returns
    -------
    ArrayImpl
        The mean squared error (MSE) calculated as the average of the squared differences between the predicted and
        true values.
    """
    linear_regress_y: ArrayImpl = linear_regression(
        feature_values,
        weight_vector,
        bias_term
    )
    return _mean_squared_error_loss(linear_regress_y, label_values)


def _mean_squared_error_loss(model_prediction: np.array, labels: np.array) -> ArrayImpl:
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
    ArrayImpl
        The mean squared error (MSE) calculated as the average of the squared differences between the predicted and
        true values.
    """
    return jnp.mean((model_prediction - labels) ** 2)
