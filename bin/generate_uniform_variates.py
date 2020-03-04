#!/usr/bin/env python


"""A command line wrapper for sampling from the distribution Q when Q is a uniform distribution."""


from selectiontest import selectiontest
import pickle, gzip
import click


@click.command()
@click.argument('n', type=int)
@click.argument('outfile', type=click.Path())
@click.option('-r', '--reps', type=int, default=10000)
def main(n, outfile, reps):
    results = selectiontest.generate_uniform_variates(n, reps, random_state=None)
    with gzip.open(outfile, 'wb') as outfile:
        pickle.dump(results, outfile)


if __name__ == "__main__":
    main()