# coding=utf-8


import numpy as np
from bisect import bisect
from scipy.special import binom
from collections import Counter
import re
from scipy.special import xlogy, gammaln
from more_itertools import locate
import functools
from selectiontest.__init__ import __version__


__author__ = "Helmut Simon"
__copyright__ = "© Copyright 2020, Helmut Simon"
__license__ = "BSD-3"
__version__ = __version__
__maintainer__ = "Helmut Simon"
__email__ = "helmut.simon@anu.edu.au"
__status__ = "Test"


@functools.lru_cache(maxsize=512)
def replace_zero(s, j):
    """Replace jth 0 in string s with a 1. First character is given by j=0 etc."""
    i = list(locate(s, lambda x: x == '0'))[j]
    return s[:i] + '1' + s[i + 1:]


@functools.lru_cache(maxsize=512)
def make_row(s):
    """Convert string of zeros (+) and ones (,) to matrix row, i.e. counting partitions by size."""
    c = re.split('1', s)
    return [int(len(x) + 1) for x in c]


def derive_tree_matrix(f):
    """Derive tree matrix from the list f. The first element of f is an integer in [0, n-1], the second
    in [0, n-2] and so on (Lehmer code)."""
    n = len(f) + 1
    s = '0' * (n - 1)
    result = list()
    for j in f:
        s = replace_zero(s, j)
        orow = make_row(s)
        urow = np.bincount(orow, minlength=n + 1)
        urow = urow[1:-1]
        result.append(urow)
    mx = np.stack(result, axis=0)
    return mx


def sample_matrix(n, reps):
    "A generator returning reps matrices for sample size n, according to the ERM distribution."
    for i in range(reps):
        f = list()
        for j in range(1, n):
            f.append(np.random.choice(j))
        f = f[::-1]
        mx = derive_tree_matrix(f)
        yield mx


def sample_branch_lengths(n, reps):
    "A generator returning reps sets of relative branch lengths for sample size n, according to the Wright-Fisher model."
    kvec = np.arange(2, n + 1, dtype=int)
    branch_lengths = np.random.exponential(scale=1 / binom(kvec, 2), size=(reps, n - 1))
    total_branch_lengths = branch_lengths @ kvec
    for row, total_length in zip(branch_lengths, total_branch_lengths):
        rel_row = row / total_length
        yield rel_row


def sample_wf_distribution(n, reps):
    """
    Calculate variates for the probability distribution Q under Wright Fisher model.

    Parameters
    ----------
    n: int
        Sample size
    reps: int
        Number of variates to generate if default is used.

    Yields
    -------
    numpy.ndarray
         Array of variates (n-1)

    """
    for mx, rel_branch_length in zip(sample_matrix(n, reps), sample_branch_lengths(n, reps)):
        variate = (mx.T).dot(rel_branch_length)
        err = 1 - np.sum(variate)
        variate[np.argmax(variate)] += err
        yield variate


def get_ERM_matrix(n):
    ERM_matrix = np.zeros((n - 1, n - 1))
    for m in range(n - 1):
        for k in range(n - 1):
            ERM_matrix[m, k] = (k + 2) * binom(n - m - 2, k) / binom(n - 1, k + 1)
    return ERM_matrix


def sample_uniform_distribution(n, reps):
    """
    Calculate variates for the uniform probability distribution Q.

    Parameters
    ----------
    n: int
        Sample size
    reps: int
        Number of variates to generate if default is used.

    Returns
    -------
    numpy.ndarray
         Array of variates (reps, n-1)

    """
    j_n = np.diag(1 / np.arange(2, n + 1))
    erm = get_ERM_matrix(n)
    avge_mx = erm.dot(j_n)
    sample = np.random.dirichlet(np.ones(n - 1), size=reps)
    variates = avge_mx @ sample.T
    return variates.T


def multinomial_pmf(x, n, p):
    x = np.array(x)
    logpmf = gammaln(n + 1) + np.sum(xlogy(x, p) - gammaln(x + 1), axis=-1)
    return np.exp(logpmf)


def test_neutrality(sfs, variates0=None, variates1=None, reps=10000):
    """
    Calculate :math:`\\rho`, the log odds ratio of the data for the distribution given by variates0 over
    the distribution given by variates1.

    Parameters
    ----------
    sfs: list
        Site frequency spectrum, e.g. [1, 3, 0, 2, 1]
    variates0: numpy array
        Array of variates from null hypothesis distribution. Default uses Wright-Fisher model.
    variates1: numpy array
        Array of variates from alternative distribution. Default uses \`uniform\' model.
    reps: int
        Number of variates to generate if default is used.

    Returns
    -------
    numpy.float64
        :math:`\\rho` (value of log odds ratio). Values can include inf, -inf or nan if one or both probabilities
        are zero due to underflow error.

    """
    n = len(sfs) + 1
    segsites = sum(sfs)
    if variates0 is None:
        variates0 = np.empty((reps, n - 1), dtype=float)
        for i, q in enumerate(sample_wf_distribution(n, reps)):
            variates0[i] = q
    if variates1 is None:
        variates1 = sample_uniform_distribution(n, reps)
    h0 = np.mean(multinomial_pmf(sfs, segsites, variates0))
    h1 = np.mean(multinomial_pmf(sfs, segsites, variates1))
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
    sfs: list
        Site frequency spectrum, e.g. [1, 3, 0, 2, 1]

    Returns
    -------
    numpy.float64
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


def generate_sfs_array(n, seg_sites, reps=10000):
    """
    Sample SFS values for Wright-Fisher model for given sample size n and conditioned on the
    number of segregating sites.

    """
    for p in sample_wf_distribution(n, reps):
        yield np.random.multinomial(seg_sites, p)


def compute_threshold(n, seg_sites, sreps=10000, wreps=10000, fpr=0.02):
    """
    Calculate threshold value of :math:`\\rho` corresponding to a given false positive rate (FPR).
    For values of :math:`\\rho` above the threshold we reject the
    null (by default neutral) hypothesis.

    Parameters
    ----------
    n: int
        Sample size
    seg_sites: int
        Number of segregating sites in sample.
    sreps: int
        Number of SFS configs and of uniform variates to generate if default is used.
    wreps: int
        Number of Wright-Fisher variates to generate if default is used.
    fpr: float
        Selected FPR tolerance.

    Returns
    -------
    numpy.float64
        Threshold value for log odds ratio

    """
    variates0 = np.empty((wreps, n - 1), dtype=float)
    for i, q in enumerate(sample_wf_distribution(n, wreps)):
        variates0[i] = q
    variates1 = sample_uniform_distribution(n, sreps)
    num_wf_vars = variates0.shape[0]
    results = list()
    for sfs in generate_sfs_array(n, seg_sites, sreps):
        a = sfs > 0
        b = variates0[:, a] > 0
        c = np.all(b > 0, axis=1)
        compat_vars = variates0[c, :]
        if compat_vars.shape[0] == 0:
            h0 = 0
        else:
            h0 = np.sum(multinomial_pmf(sfs, seg_sites, compat_vars)) / num_wf_vars
        h1 = np.mean(multinomial_pmf(sfs, seg_sites, variates1))
        rho = np.log10(h1) - np.log10(h0)
        results.append(rho)
    results = np.array(results)
    print("Count -inf: ", np.sum(np.isneginf(results)))
    print("Count  inf: ", np.sum(np.isinf(results)))
    print("Count  nan: ", np.sum(np.isnan(results)))
    results = results[~np.isnan(results)]
    results = np.sort(results)
    return results[int(len(results) * (1 - fpr))]


def calc_breakpoints(pop_sizes, timepoints):
    y = 0
    assert timepoints[0] == 0, 'Timepoints must start with zero.""'
    result = [0]
    for i in range(len(timepoints) - 1):
        y += (pop_sizes[0] / pop_sizes[i]) * (timepoints[i + 1] - timepoints[i])
        result.append(y)
    return result


def transform_coal_times(coal_times, timepoints, pop_sizes):
    """
    Transform an array of times.
    """
    breakpoints = calc_breakpoints(pop_sizes, timepoints)
    results = list()
    for ct in coal_times:
        i = bisect(breakpoints, ct)
        y = (timepoints[i - 1] + (pop_sizes[i - 1] / pop_sizes[0]) * (ct - breakpoints[i - 1]))
        results.append(y)
    return np.array(results)


def branch_lengths_to_coal_times(branch_lengths):
    blsf = np.flip(branch_lengths)
    return np.cumsum(blsf)


def coal_times_to_branch_lengths(coal_times):
    ctsf1 = np.diff(coal_times)
    ctsf = np.insert(ctsf1, 0, coal_times[0])
    return np.flip(ctsf)


def piecewise_constant_variates(n, timepoints, pop_sizes, reps=10000):
    """
    Generate variates corresponding to a piecewise constant demographic history.

    Parameters
    ----------
    n: int
        Sample size
    timepoints: array-like
        Times at which population changes (in generations, backward from the present).
    pop_sizes: array-like
        Population sizes between timepoints (only relative sizes matter.)
    reps: int
        Number of variates to generate.

    Yields
    -------
    numpy.ndarray
         Variates

    """
    zipped = zip(sample_matrix(n, reps), sample_branch_lengths(n, reps))
    for mx, branch_lengths in zipped:
        coal_times = branch_lengths_to_coal_times(branch_lengths)
        transf_coal_times = transform_coal_times(coal_times, timepoints, pop_sizes)
        transf_branch_lengths = coal_times_to_branch_lengths(transf_coal_times)
        kvec = np.arange(2, n + 1, dtype=int)
        total_branch_length = np.sum(transf_branch_lengths * kvec)
        rel_branch_lengths = transf_branch_lengths / total_branch_length
        variate = (mx.T).dot(rel_branch_lengths)
        err = 1 - np.sum(variate)
        variate[np.argmax(variate)] += err
        yield variate, transf_coal_times



def vcf2sfs(vcf_file, panel, coord, start, end, select_chr=True):
    """
    Get SFS from vcf data for given population and sequence. The panel file is used to select samples.

    Parameters
    ----------
    vcf_file: pyvcf class: Reader (https://pyvcf.readthedocs.io/en/latest/)
        Variant details

    panel: pandas DataFrame
        Proband details

    coord: str
        Coordinate (e.g. chromosome).

    start: int
        Start position of sequence.

    end: int
        End position of sequence.

    select_chr: bool
        If True, sample first chromosome. If False, use both.

    Returns
    -------
    list
        Site frequency spectrum

    int
        Sample size

    list
        Names of variants common to all elements of the sample.

    """
    n = panel.shape[0]
    if not select_chr:
        n = 2 * n
    snps = vcf_file.fetch(str(coord), start, end)
    count, anc_count = 0, 0
    allele_counts = list()
    non_seg_snps = list()
    for record in snps:
        allele_count = 0
        if record.is_snp:
            count += 1
            # Test the ancestral is one of the alleles
            if record.INFO['AA'][0] not in [record.REF, record.ALT]:
                continue
            anc_count += 1
            for sample in record.samples:
                if sample.sample in panel.index:
                    gt = sample.gt_alleles
                    if select_chr:
                        allele_count += int(gt[0])
                    else:
                        allele_count += int(gt[0]) + int(gt[1])
            if allele_count < n:    #Some SNPs may not segregate in some subpopulations.
                allele_counts.append(allele_count)
            else:
                non_seg_snps.append(record.ID)
    sfs_c = Counter(allele_counts)
    del sfs_c[0]
    sfs = np.zeros(n - 1, int)
    for i in sfs_c.keys():
        sfs[i - 1] = sfs_c[i]
    return sfs, n, non_seg_snps








