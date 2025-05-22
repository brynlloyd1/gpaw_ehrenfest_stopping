import numpy as np


def append_to_dict(dictionary, key, value):
    """
    appends to a dictionary, where the value is a list of elements

    Paramters:
    dictionary (dict): dictionary
    key (string): key
    value: (anything??): value to append to the list
    """

    if key not in dictionary:
        dictionary[key] = [value]
    else:
        dictionary[key].append(value)

    return dictionary


def calculate_r2(x, y, y_fit):
    ss_res = np.sum((y - y_fit)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot != 0 else -np.inf
    return r2


def sliding_fit(x, y, min_window_size, max_window_size=None):
    if max_window_size is None:
        max_window_size = len(x)

    # initialize variables for tracking the best fit
    best_r2 = -np.inf
    best_fit = None
    best_cov = None
    best_x = None
    best_y = None

    # iterates over window size and starting position to find the optimum fit
    for window_size in range(min_window_size, len(x)):
        for i in range(len(x) - window_size):
            x_window = x[i: i+window_size]
            y_window = y[i: i+window_size]

            fit, cov = np.polyfit(x_window, y_window, 1, cov=True)
            pfit = np.poly1d(fit)
            y_fit = pfit(x_window)

            r2 = calculate_r2(x_window, y_window, y_fit)

            if r2 > best_r2:
                best_r2 = r2
                best_fit = fit
                best_cov = cov
                best_x = x_window
                best_y = y_window

    return best_fit, best_cov, best_x, best_y
