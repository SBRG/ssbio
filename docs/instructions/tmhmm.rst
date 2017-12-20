.. _tmhmm:

*****
TMHMM
*****


Description
===========

* `TMHMM home page`_
* `TMHMM download page`_
* `TMHMM installation instructions`_

TMHMM is a program to predict the location of transmembrane helices in proteins, directly from sequence. *ssbio* provides a wrapper to execute and parse the "long" output format of TMHMM.


Instructions (Unix)
===================

.. note:: These instructions were created on an Ubuntu 17.04 system.

#. Register for the software (academic license only) at the `TMHMM download page`_
#. Receive instructions to download the software at your email address
#. Download the file *tmhmm-2.0c.Linux.tar.gz*
#. Extract it to a place where you store software
#. Install it according to the `TMHMM installation instructions`_, repeated and annotated below...
    
    #. Insert the correct path for perl 5.x in the first line of the scripts bin/tmhmm and bin/tmhmmformat.pl (if not ``/usr/local/bin/perl``). Use ``which perl`` and ``perl -v`` in the terminal to help find the correct path.
    #. Make sure you have an executable version of *decodeanhmm* in the bin directory.
    #. Include the directory containing tmhmm in your path (how do I add something to my :ref:`dummiesunix-path`?)
    #. Read the TMHMM2.0.guide.html
    #. Run the program by doing the following:
       
        .. code-block:: console

            tmhmm my_sequences.fasta


FAQs
====

* How do I cite TMHMM?

    - Krogh A, Larsson B, von Heijne G & Sonnhammer EL (2001) Predicting transmembrane protein topology with a hidden Markov model: application to complete genomes. J. Mol. Biol. 305: 567–580 Available at: http://dx.doi.org/10.1006/jmbi.2000.4315
      
* I'm having issues running TMHMM...

    - See the `ssbio wiki`_ for (hopefully) some solutions - or add yours in when you find the answer!


API
===

.. automodule:: ssbio.protein.sequence.properties.tmhmm
    :members:


.. Links
.. _TMHMM home page: http://www.cbs.dtu.dk/services/TMHMM-2.0/
.. _TMHMM download page: http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm
.. _TMHMM installation instructions: http://www.cbs.dtu.dk/services/doc/tmhmm-2.0c.readme
.. _ssbio wiki: https://github.com/SBRG/ssbio/wiki/Troubleshooting