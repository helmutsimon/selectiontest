Command line interface
======================


The python modules can be run from the comand line as follows:

::

    test_neutrality sfs -q0 variates0 -q1 variates1 -r reps

    calculate_D sfs

    generate_wf_variates n outfile -r reps

    generate_uniform_variates n outfile -r reps

``generate_wf_variates`` and ``generate_uniform_variates`` save output as gzipped pickle files with pathname ``outfile``.

These pathnames can be suppplied to ``test_neutrality`` by using the options -q0 and -q1.

Note
----

Command line scripts for other modules will be supplied.



