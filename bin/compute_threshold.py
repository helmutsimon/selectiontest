#!/usr/bin/env python

"""A command line wrapper for calculation of neutrality test statistic rho."""


from selectiontest import selectiontest
import click

@click.command()
@click.argument('n')
@click.argument('seg_sites')
@click.option('-r', '--reps', type=int, default=10000)
@click.option('-f', '--fpr', type=float, default=0.02)
def main(n, seg_sites, reps, fpr):
    print(selectiontest.compute_threshold(n, seg_sites, reps, fpr))


if __name__ == "__main__":
    main()