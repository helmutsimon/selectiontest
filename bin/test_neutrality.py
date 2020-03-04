# coding=utf-8

"""A command line wrapper for calculation of neutrality test statistic rho."""


from selectiontest import selectiontest
import pickle, gzip
import click

@click.command()
@click.argument('sfs')
@click.option('-q0', '--q0file', type=click.Path(), default=None)
@click.option('-q1', '--q1file', type=click.Path(), default=None)
@click.option('-r', '--reps', type=int, default=10000)
def main(sfs, q0file, q1file, reps):
    sfs = sfs.split(',')
    sfs = [int(x) for x in sfs]
    if q0file is not None:
        with gzip.open(q0file, 'rb') as variates0:
            variates0 = pickle.load(variates0)
    else:
        variates0 = None
    if q1file is not None:
        with gzip.open(q1file, 'rb') as variates1:
            variates1 = pickle.load(variates1)
    else:
        variates1 = None
    print(selectiontest.test_neutrality(sfs, variates0=variates0, variates1=variates1, reps=reps))


if __name__ == "__main__":
    main()