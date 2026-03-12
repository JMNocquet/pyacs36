import numpy as np
import cvxpy
from scipy.sparse import dok_matrix


def get_model_deviation_matrix(x, deviation_mapping, n_deviates):
    num_x = len(x)
    matrix = dok_matrix((num_x, n_deviates), dtype=np.float32)
    for i, xx in enumerate(x):
        j = deviation_mapping(xx)
        assert isinstance(j, int)
        assert j < n_deviates
        matrix[i, j] = 1.0

    return matrix


def complete_linear_deviation(linear_deviation, x, default_name):
    assert 'n_vars' in linear_deviation
    lin_dev = linear_deviation.copy()

    if 'name' not in lin_dev:
        lin_dev['name'] = default_name

    if 'alpha' not in lin_dev:
        lin_dev['alpha'] = 1e-3

    assert lin_dev['alpha'] >= 0.0

    if 'matrix' not in linear_deviation:
        assert 'mapping' in lin_dev
        matrix = get_model_deviation_matrix(x, lin_dev['mapping'], lin_dev['n_vars'])
        lin_dev['matrix'] = matrix

    if 'positive' not in lin_dev:
        lin_dev['positive'] = False

    lin_dev['variable'] = cvxpy.Variable(lin_dev['n_vars'], pos=lin_dev['positive'])
    lin_dev['model_contribution'] = lin_dev['matrix'] @ lin_dev['variable']

    return lin_dev


def complete_linear_deviations(linear_deviations, x):
    completed_devs = []
    names = set()
    for i, linear_deviation in enumerate(linear_deviations):
        default_name = 'linear_deviation_%s' % i
        completed = complete_linear_deviation(linear_deviation, x, default_name)

        # ensure unique names
        assert completed['name'] not in names
        names.add(completed['name'])

        completed_devs.append(completed)

    return completed_devs
