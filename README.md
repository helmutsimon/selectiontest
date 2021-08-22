# selectiontest
## Introduction
This repository contains python modules and command line scripts to support testing for selective neutrality using relative likelihood. It enables application of the methods described in Simon and Huttley 2021 *A New Likelihood-based Test for Natural Selection* bioRxiv doi = 10.1101/2021.07.04.451068. For other software supporting that paper, but not required for applications, see https://github.com/helmutsimon/NeutralityTest.
Functions available include:
* calculate statistic for relative neutrality, &rho;, which is a relative likelihood of two models;
* generate variates corresponding to the Wright-Fisher model and the `uniform distribution' model;
* calculate Tajima's D (for comparison purposes);
* calibrate &rho;, that is, find the threshold corresponding to a desired false positive (Type I error) rate; and
* generate variates corresponding to a piece-wise constant demographic history.


## Installation
Installation is currently from test.pypi. Use:
python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps selectiontest-helmutsimon --upgrade
## Documentation
See https://selectiontest.readthedocs.io/en/latest/index.html
