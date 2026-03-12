from scipy.interpolate import interp1d
import numpy as np


def vectorize(func_orig):
    """
    A function that takes a function and
    returns another that can fun on lists and arrays
    :param func_orig: any functions
    :return: vectorized function
    """
    def func(x_new):
        if isinstance(x_new, list):
            return [func_orig(xx) for xx in x_new]

        if isinstance(x_new, np.ndarray):
            return np.array([func_orig(xx) for xx in x_new])

        return func_orig(x_new)
    return func


def get_interp_extrapolate_functions(x, base_model, linear_deviations):
    """
    Get the three interp/extrapolation model functions:
        base function, deviates function, total model function
    :param x: the x data
    :param base_model: model model cvxpy expression
    :param linear_deviations: list of completed linear_deviations objects
    :return:  base function, deviates function, total model function
    """

    # TODO: this requires mapping to be given, make it work with matrix only
    interp_base_model_func = interp1d(x, base_model.value, fill_value="extrapolate")

    def func_base(x_new):
        return interp_base_model_func(x_new)

    def func_deviates(x_new):
        linear_dev_value = 0.0
        for lin_dev in linear_deviations:
            index = lin_dev['mapping'](x_new)
            var = lin_dev['variable'].value
            value = var[index]
            linear_dev_value += value

        return linear_dev_value

    def func(x_new):
        return func_base(x_new) + func_deviates(x_new)

    return vectorize(func_base), vectorize(func_deviates), vectorize(func)
