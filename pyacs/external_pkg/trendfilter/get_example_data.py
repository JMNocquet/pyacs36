import numpy as np


def get_example_data():
    noise = 0.2
    np.random.seed(420)

    # make a set of x, y points
    # y = sqrt(x) plus noise

    x = np.linspace(0, 10, 80)
    n = len(x)
    y = np.sqrt(x)
    y_noisy = y + noise * np.random.randn(n)

    # add an outlier points
    y_noisy[20] += 3
    y_noisy[60] += 2

    return x, y_noisy


def deviation_mapping(index):
    # assume data points are monthly values
    # maps the x data point to the right linear
    index = int(index)
    # deviation variable
    return index % 12


def get_example_data_seasonal():
    noise = 0.3
    # noise = 0.0
    np.random.seed(400)

    # make a set of x, y points
    # signal plus noise, plus seasonality

    x = np.arange(100)
    n = len(x)
    y = 4.0 + np.sqrt(x) - 0.1 * x

    # choose 12 random deviates
    period = 10.0
    seasonal_values = np.random.randn(12) + 1.1 * np.sin(2 * np.pi * np.arange(12)/period)
    seasonal_values *= 0.3
    seasonal_values[1] = seasonal_values[1] * 3
    seasonal_indices = np.array([deviation_mapping(i) for i in x])
    seasonal_offsets = seasonal_values[seasonal_indices]

    noise_term = noise * np.random.randn(n)

    y_noisy = y + seasonal_offsets + noise_term

    # add an outlier points
    y_noisy[20] += 6
    y_noisy[60] += 5

    return x, y_noisy
