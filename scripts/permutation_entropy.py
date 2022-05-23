import numpy as np
# import itertools
#
#
# def permutation_entropy(time_series, m, delay):
#     """Calculate the Permutation Entropy
#
#     Args:
#         time_series: Time series for analysis
#         m: Order of permutation entropy
#         delay: Time delay
#
#     Returns:
#         Vector containing Permutation Entropy
#
#     Reference:
#         [1] Massimiliano Zanin et al. Permutation Entropy and Its Main Biomedical and Econophysics Applications:
#             A Review. http://www.mdpi.com/1099-4300/14/8/1553/pdf
#         [2] Christoph Bandt and Bernd Pompe. Permutation entropy — a natural complexity
#             measure for time series. http://stubber.math-inf.uni-greifswald.de/pub/full/prep/2001/11.pdf
#         [3] http://www.mathworks.com/matlabcentral/fileexchange/37289-permutation-entropy/content/pec.m
#     """
#     n = len(time_series)
#     permutations = np.array(list(itertools.permutations(range(m))))
#     c = [0] * len(permutations)
#
#     for i in range(n - delay * (m - 1)):
#         # sorted_time_series =    np.sort(time_series[i:i+delay*m:delay], kind='quicksort')
#         sorted_index_array = np.array(np.argsort(time_series[i:i + delay * m:delay], kind='quicksort'))
#         for j in range(len(permutations)):
#             if abs(permutations[j] - sorted_index_array).any() == 0:
#                 c[j] += 1
#
#     c = [element for element in c if element != 0]
#     p = np.divide(np.array(c), float(sum(c)))
#     pe = -sum(p * np.log(p)) / np.log(np.math.factorial(m))
#     return pe
from math import factorial


def perm_entropy(x, order=3, delay=1, normalize=False):
    ''' Permutation Entropy.

    Parameters
    ----------
    x : list or np.array
        One-dimensional time series of shape (n_times)
    order : int
        Order of permutation entropy. Default is 3.
    delay : int
        Time delay (lag). Default is 1.
    normalize : bool
        If True, divide by log2(order!) to normalize the entropy between 0
        and 1. Otherwise, return the permutation entropy in bit.

    Returns
    -------
    pe : float
        Permutation Entropy.

    Notes
    -----
    The permutation entropy is a complexity measure for time-series first
    introduced by Bandt and Pompe in 2002.

    The permutation entropy of a signal :math:`x` is defined as:

    .. math:: H = -\\sum p(\\pi)\\log_2(\\pi)

    where the sum runs over all :math:`n!` permutations :math:`\\pi` of order
    :math:`n`. This is the information contained in comparing :math:`n`
    consecutive values of the time series. It is clear that
    :math:`0 ≤ H (n) ≤ \\log_2(n!)` where the lower bound is attained for an
    increasing or decreasing sequence of values, and the upper bound for a
    completely random system where all :math:`n!` possible permutations appear
    with the same probability.

    The embedded matrix :math:`Y` is created by:

    .. math::
        y(i)=[x_i,x_{i+\\text{delay}}, ...,x_{i+(\\text{order}-1) *
        \\text{delay}}]

    .. math:: Y=[y(1),y(2),...,y(N-(\\text{order}-1))*\\text{delay})]^T

    References
    ----------
    Bandt, Christoph, and Bernd Pompe. "Permutation entropy: a
    natural complexity measure for time series." Physical review letters
    88.17 (2002): 174102.


    '''
    x = np.array(x)
    ran_order = range(order)
    hashmult = np.power(order, ran_order)
    # Embed x and sort the order of permutations
    sorted_idx = _embed(x, order=order, delay=delay).argsort(kind='quicksort')
    # Associate unique integer to each permutations
    hashval = (np.multiply(sorted_idx, hashmult)).sum(1)
    # Return the counts
    _, c = np.unique(hashval, return_counts=True)
    # Use np.true_divide for Python 2 compatibility
    p = np.true_divide(c, c.sum())
    pe = -np.multiply(p, np.log2(p)).sum()
    if normalize:
        pe /= np.log2(factorial(order))
    return pe