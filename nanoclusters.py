import numpy as np
import pandas as pd
import scipy as scp


def right_quantile_weight(arr, q):
    # Numerator of the right quantile weight
    num = np.quantile(arr, (1+q)/2) + np.quantile(arr, 1 -(q)/2) - 2*np.quantile(arr, 0.75)
    # Denominator of the right quantile weight
    den = np.quantile(arr, (1 + q)/2) - np.quantile(arr, 1 - q/2)

    return num/den


def left_quantile_weight(arr ,p):
    # Numerator of the left quantile weight
    num = np.quantile(arr, (1 - p)/2) + np.quantile(arr, p/2) - 2*np.quantile(arr, 0.25)
    # Denominator of the right quantile weight
    den = np.quantile(arr, (1 - p)/2) - np.quantile(arr, p/2)

    return -num/den