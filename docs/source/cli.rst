Command line interface
======================


The python modules can be run from the comand line as follows:

::

    test_neutrality sfs -q0 variates0 -q1 variates1 -r reps

    calculate_D sfs

    sample_wf_distribution n outfile -r reps

    sample_uniform_distribution n outfile -r reps

    compute_threshold n seg_sites -r reps -f fpr

``sample_wf_distribution`` and ``sample_uniform_distribution`` save output as gzipped pickle files with pathname ``outfile``. These pathnames can be supplied to ``test_neutrality`` by using the options -q0 and -q1.

Note
----

At present, there is no command line script for the python module piecewise_constant_variates.



