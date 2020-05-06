import click
import numpy as np
from selectiontest import selectiontest


@click.group()
def selectiontestcli():
    pass


@selectiontestcli.command()
@click.argument('sfs', nargs=-1)
def calculate_D(sfs):
    """Calculate Tajima's D from SFS (integers separated by spaces)."""
    sfs = np.asarray(sfs, dtype=float)
    tajd = selectiontest.calculate_D(sfs)
    click.echo(tajd)


@selectiontestcli.command()
@click.argument('sfs', nargs=-1)
@click.option('-r', '--reps', default=10000, help='Number of samples for Monte Carlo integration.')
def test_neutrality(sfs, reps):
    """Calculate the log odds ratio of the data SFS (integers separated by spaces)."""
    sfs = np.asarray(sfs, dtype=float)
    rho = selectiontest.test_neutrality(sfs, variates0=None, variates1=None, reps=reps)
    click.echo(rho)


@selectiontestcli.command()
@click.argument('n', type=int)
@click.argument('seg_sites', type=int)
@click.option('-r', '--reps', default=10000, help='Number of test samples.')
@click.option('-f', '--fpr', default=0.02, help='Value of false positive rate set for threshold.')
def compute_threshold(n, seg_sites, reps, fpr):
    """Calculate threshold value of log odds ratio for sample size N and numver of segregating sites SEG_SITES."""
    threshold = selectiontest.compute_threshold(n, seg_sites, reps, fpr)
    click.echo(threshold)