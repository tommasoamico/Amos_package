
import numpy as np


def min_method(array):
    return(np.max(abs(array)))

def info_min_method():
    print('''
    This method returns the maximum absolute value of the array provided in input
    (thaught to choose the line that returns the minimum value of this method)
    Takes in input,
    1) The array of residues
    ''')

def sum_method(array):
    return np.sum(abs(array))


def info_sum_method():
    print('''
    This method returns the sum of the absolute values of the array provided in input
    Takes in input,
    1) The array of residues
    ''')