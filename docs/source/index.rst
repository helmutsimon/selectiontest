

Welcome to selectiontest's documentation!
=========================================

Overview
--------
This repository contains python modules and command line scripts to support testing for selective neutrality using relative likelihood. It enables application of the methods described in Simon and Huttley 2021 *A New Likelihood-based Test for Natural Selection* bioRxiv doi = 10.1101/2021.07.04.451068. For other software supporting that paper, but not required for applications, see https://github.com/helmutsimon/NeutralityTest.

Functions available include:


* calculate statistic for relative neutrality, :math:`{\rho }`, which is a relative likelihood of two models;
* generate variates corresponding to the Wright-Fisher model and the \`uniform distribution\' model;
* calculate Tajima's D (for comparison purposes);
* calibrate :math:`{\rho }`, that is, find the threshold corresponding to a desired false positive (Type I error) rate; and
* generate variates corresponding to a piece-wise constant demographic history.
* compute SFS from variant data in vcf format (pyvcf class: Reader (https://pyvcf.readthedocs.io/en/latest/)

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   modules
   cli
   license

Source
------

Source code is hosted at https://github.com/helmutsimon/selectiontest.

Citation
--------
For the present, cite github repository. Bibtex:

::

    @misc{simon_selection,
    author = {Simon, Helmut and Huttley, Gavin},
    title = {selectiontest 0.1.10},
    year = {2020},
    url = {https://github.com/helmutsimon/selectiontest},
    note = {[https://github.com/helmutsimon/selectiontest]}}

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
