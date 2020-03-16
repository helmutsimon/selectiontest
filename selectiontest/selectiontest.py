# coding=utf-8


import numpy as np
import sys
from bisect import bisect
from math import log, factorial
from scipy.special import binom
from scipy.stats import multinomial, expon
from scipy.stats import dirichlet


__author__ = "Helmut Simon"
__copyright__ = "© Copyright 2020, Helmut Simon"
__license__ = "BSD-3"
__version__ = "0.1.0"
__maintainer__ = "Helmut Simon"
__email__ = "helmut.simon@anu.edu.au"
__status__ = "Test"



def get_ERM_matrix(n):
    ERM_matrix = np.zeros((n - 1, n - 1))
    for m in range(n - 1):
        for k in range(n - 1):
            ERM_matrix[m, k] = (k + 2) * binom(n - m - 2, k) / binom(n - 1, k + 1)
    return ERM_matrix


def generate_wf_variates(n, reps, random_state=None):
    """
    Calculate variates for the probability distribution Q under Wright Fisher model.

    Parameters
    ----------
    n: int
        Sample size
    reps: int
        Number of variates to generate if default is used.
    random_state: int
        Value used to seed local RandomState instance for generating variates from exponential distribution.

    Returns
    -------
    numpy.ndarray
         Array of variates (reps, n-1)

    """
    erm = get_ERM_matrix(n)
    kvec = np.arange(2, n + 1, dtype=int)
    branch_lengths = expon.rvs(scale=1 / binom(kvec, 2), size=(reps, n - 1), random_state=random_state)
    total_branch_lengths = branch_lengths @ kvec
    rel_branch_lengths = list()
    for row, total_length in zip(branch_lengths, total_branch_lengths):
        rel_row = row / total_length
        rel_branch_lengths.append(rel_row)
    rel_branch_lengths = np.array(rel_branch_lengths)
    variates = (erm @ rel_branch_lengths.T).T
    return variates


def generate_uniform_variates(n, reps, random_state=None):
    """
    Calculate variates for the uniform probability distribution Q.

    Parameters
    ----------
    n: int
        Sample size
    reps: int
        Number of variates to generate if default is used.
    random_state: int
        Value used to seed local RandomState instance for generating variates from Dirichlet distribution.

    Returns
    -------
    numpy.ndarray
         Array of variates (reps, n-1)

    """
    j_n = np.diag(1 / np.arange(2, n + 1))
    erm = get_ERM_matrix(n)
    avge_mx = erm.dot(j_n)
    sample = dirichlet.rvs(np.ones(n - 1), size=reps, random_state=random_state)
    variates = avge_mx @ sample.T
    return variates.T


def multinomial_pmf(counts, probs):
    """
    Calculate PMF of multinomial distribution. Number of draws is the sum of counts.
    probs can be a 2D array, with each row totalling 1.

    """
    xx = list()
    for row in probs:
        new_row = np.log(row) * counts
        xx.append(new_row)
    xx = np.array(xx)
    counts = np.array(counts)
    xx[:, np.where(counts == 0.)] = 0.   # otherwise x will contain NaN even if count = 0
    x = np.sum(xx, axis=1)
    y = np.sum([log(factorial(i)) for i in counts])
    z = np.sum(np.log(np.arange(1, sum(counts) + 1)))
    mult_prob = np.exp(x + z - y)
    return mult_prob


def test_neutrality(sfs, variates0=None, variates1=None, reps=10000):
    """
    Calculate :math:`\\rho`, the log odds ratio of the data for the distribution given by variates0 over
    the distribution given by variates1.

    Parameters
    ----------
    sfs: string
        Site frequency spectrum separated by commas, e.g. 1,3,0,2,1
    variates0: numpy array
        Array of variates from null hypothesis distribution. Default uses Wright-Fisher model.
    variates1: numpy array
        Array of variates from null hypothesis distribution. Default uses \`uniform\' model.
    reps: int
        Number of variates to generate if default is used.

    Returns
    -------
    float
        :math:`\\rho` (value of log odds ratio)

    """
    n = len(sfs) + 1
    if variates0 is None:
        variates0 = generate_wf_variates(n, reps)
    if variates1 is None:
        variates1 = generate_uniform_variates(n, reps)
    h0 = np.mean(multinomial_pmf(sfs, variates0))
    h1 = np.mean(multinomial_pmf(sfs, variates1))
    if h0 == 0 or h1 == 0:
        print(sfs, 'h0 = ', h0, 'h1 = ', h1)
        if h0 != 0:
            h1 = sys.float_info.min
    return np.log10(h1) - np.log10(h0)


def pi_calc(sfs):
    """
    Calculate the mean number of pairwise differences from a site frequency spectrum.

    """
    sfs = np.array(sfs)
    n = len(sfs) + 1
    g1 = np.arange(1, n)
    g2 = n - g1
    g3 = g1 * g2 * sfs
    pi = np.sum(g3) * 2 / (n * (n - 1))
    return pi


def calculate_D(sfs):
    """
    Calculate Tajima's D from a site frequency spectrum.

    Parameters
    ----------
    sfs: string
        Site frequency spectrum separated by commas, e.g. 1,3,0,2,1

    Returns
    -------
    float
        Value of Tajima\'s D.

    """
    seg_sites = np.sum(sfs)
    pi = pi_calc(sfs)
    n = len(sfs) + 1
    n_seq = np.arange(1, n)
    a1 = (1 / n_seq).sum()
    a2 = ((1 / n_seq) ** 2).sum()
    b1 = (n + 1) / (3 * (n - 1))
    b2 = (2 * (n ** 2 + n + 3)) / (9 * n * (n - 1))
    c1 = b1 - (1 / a1)
    c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / a1 ** 2)
    e1 = c1 / a1
    e2 = c2 / (a1 ** 2 + a2)
    tajD = (pi - (seg_sites / a1)) / np.sqrt(e1 * seg_sites + e2 * seg_sites * (seg_sites - 1))
    return tajD


def mul(seg_sites):
    def multinom(p):
        return multinomial.rvs(seg_sites, p)

    return multinom


def generate_sfs_array(n, seg_sites, reps=10000):
    """
    Sample SFS values for Wright-Fisher model for given sample size n and conditioned on the
    number of segregating sites.

    """
    variates = generate_wf_variates(n, reps)
    sfs_array = np.apply_along_axis(mul(seg_sites), 1, variates)
    return sfs_array


def compute_threshold(n, seg_sites, reps=10000, fpr=0.02):
    """
    Calculate threshold value of :math:`\\rho` corresponding to a given false positive rate (FPR).
    For values of :math:`\\rho` below the threshold we reject the
    null (by default neutral) hypothesis.

    Parameters
    ----------
    n: int
        Sample size
    seg_sites: int
        Number of segregating sites in sample.
    reps: int
        Number of variates to generate if default is used.
    fpr: float
        Selected FPR tolerance.

    Returns
    -------
    float
        Threshold value (upper) for log odds ratio

    """

    sfs_array = generate_sfs_array(n, seg_sites, reps)
    results = np.apply_along_axis(test_neutrality, 1, sfs_array)
    results = np.sort(results)
    results = results[~np.isnan(results)]
    return results[int(len(results) * fpr)]


def calc_breakpoints(pop_sizes, timepoints):
    y = 0
    result = [0]
    for i in range(0, len(timepoints) - 1):
        y += (pop_sizes[0] / pop_sizes[i]) * (timepoints[i + 1] - timepoints[i])
        result.append(y)
    return result


def calc_branch_length2(pop_sizes, timepoints):
    def calc_branch_length3(branch):
        breakpoints = calc_breakpoints(pop_sizes, timepoints)
        i = bisect(breakpoints, branch)
        return timepoints[i - 1] + (pop_sizes[i - 1] / pop_sizes[0]) * (branch - breakpoints[i - 1])

    return calc_branch_length3


def piecewise_constant_variates(n, timepoints, pop_sizes, reps=10000):
    """
    Generate variates corresponding to a piecewise constant demographic ghistory.

    Parameters
    ----------
    n: int
        Sample size
    timepoints: float
        Times at which population changes (in generations, backward from the present).
    pop_sizes: float
        Population sizes between timepoints (only relative sizes matter.)
    reps: int
        Number of variates to generate.

    Returns
    -------
    numpy.ndarray
         Array of variates

    """
    variates = generate_wf_variates(n, reps)
    branches = np.flip(variates, axis=1)
    s_k = np.cumsum(branches, axis=1)
    func1 = calc_branch_length2(pop_sizes, timepoints)
    vfunc = np.vectorize(func1)
    coal_times = vfunc(s_k)
    temp = np.diff(coal_times, axis=1)
    col2 = coal_times[:, 0]
    col2.shape = (col2.shape[0], 1)
    temp1 = np.hstack([col2, temp])
    return np.flip(temp1, axis=1)




