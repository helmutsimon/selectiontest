import click
import numpy as np
from selectiontest import selectiontest


@click.group()
def selectiontestcli():
    pass


@selectiontestcli.command()
@click.argument('sfs', nargs=-1)
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
    sfs = np.asarray(sfs, dtype=float)
    tajd = selectiontest.calculate_D(sfs)
    click.echo(tajd)


@selectiontestcli.command()
@click.argument('sfs', nargs=-1)
@click.option('-r', '--reps', default=10000, help='Number of samples for Monte Carlo integration.')
def test_neutrality(sfs, reps):
    """
    Calculate the log odds ratio of the data for the distribution given by variates0 over
    the distribution given by variates1.

    Parameters
    ----------

    sfs: int
        Site frequency spectrum, e.g. 1 3 0 2 1

    reps: int
        Number of variates to generate if default is used.

    Returns
    -------

    numpy.float64
        Value of log odds ratio. Values can include inf, -inf or nan if one or both probabilities are zero due to underflow error.

    """
    sfs = np.asarray(sfs, dtype=float)
    rho = selectiontest.test_neutrality(sfs, variates0=None, variates1=None, reps=reps)
    click.echo(rho)


@selectiontestcli.command()
@click.argument('n', type=int)
@click.argument('seg_sites', type=int)
@click.option('-r', '--reps', default=10000, help='Number of test samples.')
@click.option('-f', '--fpr', default=0.02, help='Value of false positive rate set for threshold.')
def compute_threshold(n, seg_sites, reps, fpr):
    """
    Calculate threshold value of :math:`\\rho` corresponding to a given false positive rate (FPR).
    For values of :math:`\\rho` above the threshold we reject the
    null (by default neutral) hypothesis.

    Parameters
    ----------
    n: int
        Sample size.

    seg_sites: int
        Number of segregating sites in sample.

    reps: int
        Number of variates to generate if default is used.

    fpr: float
        Selected FPR tolerance.

    Returns
    -------
    numpy.float64
        Threshold value for log odds ratio
    """
    threshold = selectiontest.compute_threshold(n, seg_sites, reps, fpr)
    click.echo(threshold)