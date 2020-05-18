import click
import numpy as np
import pandas as pd
import pysam
from vcf import Reader        # https://pypi.org/project/PyVCF/
from selectiontest import selectiontest
from selectiontest.__init__ import __version__


__author__ = "Helmut Simon"
__copyright__ = "© Copyright 2020, Helmut Simon"
__license__ = "BSD-3"
__version__ = __version__
__maintainer__ = "Helmut Simon"
__email__ = "helmut.simon@anu.edu.au"
__status__ = "Test"


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


@selectiontestcli.command()
@click.argument('vcf_name', type=click.Path('rb'))
@click.argument('panel_name',  type=click.Path('rb'))
@click.argument('chrom', type=int)
@click.argument('start', type=int)
@click.argument('end', type=int)
@click.option('--sel', type=(str, str), default=(None, None), help='Select from panel column (column name, value)')
@click.option('-r', '--reps', default=10000, help='Number of samples for Monte Carlo integration.')
@click.option('-s', '--select_chr', default=True, type=bool)
def test_neutrality_from_vcf(vcf_name, panel_name, chrom, start, end, sel, reps, select_chr):
    """Calculate the log odds ratio of the data computer from pPyVCF file VCF_NAME, proband details PANEL_NAME and
    region defined by CHROM, START and END."""
    vcf_file = Reader(filename=vcf_name, compressed=True, encoding='utf-8')
    panel = pd.read_csv(panel_name, sep=None, engine='python', skipinitialspace=True, index_col=0)
    if sel != (None, None):
        panel = panel[panel[sel[0]] == sel[1]]
    sfs, n, non_seg_snps = selectiontest.vcf2sfs(vcf_file, panel, chrom, start, end, select_chr)
    rho = selectiontest.test_neutrality(sfs, variates0=None, variates1=None, reps=reps)
    click.echo(rho)