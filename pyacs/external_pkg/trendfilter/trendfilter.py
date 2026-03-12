import numpy as np
import cvxpy
from trendfilter.extrapolate import get_interp_extrapolate_functions
from trendfilter.derivatives import second_derivative_matrix_nes, \
    first_derv_nes_cvxpy
from trendfilter.linear_deviations import complete_linear_deviations
from cvxpy import CLARABEL

default_solver = CLARABEL
# In newer versions of cvxpy, to use ECOS, need to install the ecos package
# Clarabel is thew new built-in solver


def trend_filter(x, y, y_err=None, alpha_1=0.0,
                 alpha_2=0.0, l_norm=2,
                 constrain_zero=False, monotonic=False,
                 positive=False,
                 linear_deviations=None,
                 solver=default_solver):
    """
    :param x: The x-value, numpy array
    :param y: The y variable, numpy array
    :param y_err: The y_err variable, numpy array
        Default to 1
     :param alpha_1: Regularization against non-zero slope (first derivative)
        Setting this very high will result in stair step model (if L1)
    :param alpha_2: Regularization against (second derivative or changing slope)
        Setting this very high will result in piecewise linear model (if L1)
     :param l_norm: 1 or 2 to use either L1 or L2 norm
    :param constrain_zero: If True constrains the model to be zero at origin
        Default False
    :param monotonic: If set to True, will result in a monotonically
        increasing function. Default is False.
    :param positive: If set to True, base model will be positive.
        Default False
    :param linear_deviations: list of linear deviation objects
    :param solver: solver_name, check cvxpy.installed_solvers()
        for list of installed solvers
    :return: The fit model information
    """

    # ensure input is numpy array not list or tuple
    x = np.array(x)
    y = np.array(y)

    if linear_deviations is None:
        linear_deviations = []

    linear_deviations = complete_linear_deviations(linear_deviations, x)

    assert l_norm in [1, 2]
    n = len(x)

    # get the y_err is not supplied
    assert len(y) == n
    if y_err is None:
        y_err = np.ones(n)
    else:
        assert len(y_err) == n

    # the objective function
    result = get_obj_func_model(y, y_err=y_err,
                                positive=positive,
                                linear_deviations=linear_deviations)

    derv_1 = first_derv_nes_cvxpy(x, result['base_model'])

    # the regularization
    reg_sum, regs = get_reg(x, result['base_model'], derv_1, l_norm, alpha_1, alpha_2,
                            linear_deviations=linear_deviations)

    # the total objective function with regularization
    obj_total = result['objective_function'] + reg_sum

    # The objective
    obj = cvxpy.Minimize(obj_total)

    # Get the constraints if any
    constraints = []
    if constrain_zero:
        constraints.append(result['model'][0] == 0)

    if monotonic:
        # TODO:
        # should this be derv of base or model?
        constraints.append(derv_1 >= 0)

    # define and solve the problem
    problem = cvxpy.Problem(obj, constraints=constraints)
    problem.solve(solver=solver)

    func_base, func_deviates, func = \
        get_interp_extrapolate_functions(x, result['base_model'], linear_deviations)

    tf_result = {'x': x,
                 'y': y,
                 'y_err': y_err,
                 'function': func,
                 'function_base': func_base,
                 'function_deviates': func_deviates,
                 'model': result['model'],
                 'base_model': result['base_model'],
                 'objective_model': result['objective_function'],
                 'regularization_total': reg_sum,
                 'regularizations': regs,
                 'objective_total': obj,
                 'y_fit': result['model'].value,
                 'constraints': constraints,
                 'linear_deviations': linear_deviations}

    return tf_result


def get_reg(x, base_model, derv_1, l_norm, alpha_1, alpha_2, linear_deviations=None):
    """
    Get the regularization term
    :param x: The x-value, numpy array
    :param base_model: The y variable, cvxpy.Variable(n)
    :param derv_1: the first derivative cvxpy expression from first_derv_nes_cvxpy
    :param l_norm: 1 or 2 to use either L1 or L2 norm
    :param alpha_1: Regularization against non-zero slope (first derivative)
        Setting this very high will result in stair step model (if L1)
    :param alpha_2: Regularization against (second derivative or changing slope)
        Setting this very high will result in piecewise linear model (if L1)
    :param linear_deviations: list of linear deviation objects
    :return: (sum of regs, list of regs)
    """
    d2 = second_derivative_matrix_nes(x, scale_free=True)

    if l_norm == 2:
        norm = cvxpy.sum_squares
    else:
        norm = cvxpy.norm1

    reg_1 = alpha_1 * norm(derv_1)
    reg_2 = alpha_2 * norm(d2 @ base_model)
    regs = [reg_1,  reg_2]

    for lin_dev in linear_deviations:
        reg = lin_dev['alpha'] * norm(lin_dev['variable'])
        regs.append(reg)

    reg_sum = sum(regs)

    return reg_sum, regs


def get_obj_func_model(y, y_err=None, positive=False, linear_deviations=None):
    """
    Get the objective function and the model as cvxpy expressions
    :param y: The y variable, numpy array
    :param y_err: The y_err variable, numpy array
        Default to 1
    :param positive: If set to True, will result in an always positive base model.
        Default is False.
    :param linear_deviations: List of completed linear deviation objects
    :return: objective function and the model
    """

    if linear_deviations is None:
        linear_deviations = []

    n = len(y)
    if y_err is None:
        y_err = np.ones(n)
    else:
        assert len(y_err) == n
        y_err = np.array(y_err)

    # prevent divide by zero with buffer
    buff = 0.01 * np.median(abs(y))
    buff_2 = buff ** 2
    isig = 1 / np.sqrt(buff_2 + y_err ** 2)

    base_model = cvxpy.Variable(n, pos=positive)

    model = base_model

    for lin_dev in linear_deviations:
        model += lin_dev['model_contribution']

    diff = cvxpy.multiply(isig, model - y)
    obj_func = cvxpy.sum(cvxpy.huber(diff))

    result = {'base_model': base_model,
              'model': model,
              'objective_function': obj_func}

    return result
