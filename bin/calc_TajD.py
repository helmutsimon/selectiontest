# coding=utf-8

"""A command line wrapper for calculation of neutrality test statistic rho."""


from selectiontest import selectiontest
import click

@click.command()
@click.argument('sfs')
def main(sfs):
    sfs = sfs.split(',')
    sfs = [int(x) for x in sfs]
    print(selectiontest.calculate_D(sfs))


if __name__ == "__main__":
    main()