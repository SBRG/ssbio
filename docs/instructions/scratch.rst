.. _scratch:

**********************************
SCRATCH
**********************************

.. image:: ../assets/ssbio_protein_props_scratch.png
    :align: center
    :alt: Secondary structure
    :scale: 30 %


Description
===========

* `SCRATCH home page`_
* `SCRATCH download page (for SSpro and ACCpro)`_

SCRATCH is a suite of tools to predict many types of structural properties directly from sequence. *ssbio* contains wrappers to execute and parse results from *SSpro*/*SSpro8* - predictors of secondary structure, and *ACCpro*/*ACCpro20* - predictors of solvent accessibility.


Installation instructions (Ubuntu)
==================================

.. note:: These instructions were created on an Ubuntu 17.04 system.

#. Download the source and install it using the perl script:

    .. code-block:: console

        mkdir /path/to/my/software/scratch
        cd /path/to/my/software/scratch
        wget http://download.igb.uci.edu/SCRATCH-1D_1.1.tar.gz
        tar -zxf SCRATCH-1D_1.1.tar.gz
        cd SCRATCH-1D_1.1
        perl install.pl

#. To run it from the command line directly:
    
    .. code-block:: console

        

#. *ssbio* also provides command line wrappers to run it and parse the results, see  for details.


Program execution
=================

In the shell
------------

To run the program on its own in the shell...
   
    .. code-block:: console

        /path/to/my/software/scratch/SCRATCH-1D_1.1/bin/run_SCRATCH-1D_predictors.sh  input_fasta  output_prefix  [num_threads]

With *ssbio*
------------

To run the program using the ssbio Python wrapper, see: :func:`ssbio.protein.sequence.properties.scratch.SCRATCH.run_scratch`
  

FAQs
====

* How do I cite SCRATCH?

    - Cheng J, Randall AZ, Sweredoski MJ & Baldi P (2005) SCRATCH: a protein structure and structural feature prediction server. Nucleic Acids Res. 33: W72–6 Available at: http://dx.doi.org/10.1093/nar/gki396
      
* I'm having issues running STRIDE...

    - See the `ssbio wiki`_ for (hopefully) some solutions - or add yours in when you find the answer!


API
===

.. automodule:: ssbio.protein.sequence.properties.scratch
    :members:


.. Links
.. _SCRATCH home page: http://scratch.proteomics.ics.uci.edu/
.. _SCRATCH download page (for SSpro and ACCpro): http://download.igb.uci.edu/#sspro
.. _ssbio wiki: https://github.com/SBRG/ssbio/wiki/Troubleshooting