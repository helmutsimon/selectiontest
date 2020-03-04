# selectiontest
## Introduction
This repository contains python modules and command line scripts to support testing for selective neutrality using relative likelihood. It enables application of the methods described in Simon and Huttley A New Statistical Test Provides Evidence of Selection Against Deleterious Mutations in Genes Promoting Disease Resistance (in preparation). For other software supporting that paper, but not required for applications, see https://github.com/helmutsimon/NeutralityTest.

Functions available include:
* calculate statistic for relative neutrality, &rho;, which is a relative likelihood of two models;
* generate variates corresponding to the Wright-Fisher model and the `uniform distribution' model;
* calculate Tajima's D (for comparison purposes);
* calibrate &rho;, that is, find th threshold corresponding to a desired false positive (Type I error) rate; and
* generate variates corresponding to a piece-wise constant demographic history.


## Installation
Installation is currently from test.pypi. Use:
python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps selectiontest-helmutsimon --upgrade
## Documentation
To be advised.
