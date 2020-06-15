Command line interface
======================


The following python modules can be run from the bash command line, using the prefix command ``st'' as shown:

::

    st test-neutrality sfs -r reps

    st calculate-d sfs

    st compute-threshold n seg_sites -r reps -f fpr

For the command line interface, sfs is entered as a string of integers separated by spaces (variadic argument), e.g. 1 3 0 2 1.

For more information, use ``st --help'', ``st test-neutrality --help'' etc.

In addition, there is a command st test-neutrality-from-vcf to calculate :math:`{\rho }` directly from VCF data, combining the modules vcf2sfs and test_neutrality.

For details use command ``st test-neutrality-from-vcf --help''.

This requires vcf data and sample panel details in the format used by the 1000 Genomes Project at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/.

