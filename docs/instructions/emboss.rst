.. _emboss:

******
EMBOSS
******


Description
===========

* `EMBOSS home page`_
* `EMBOSS source code`_

EMBOSS is the European Molecular Biology Open Software Suite. EMBOSS contains a wide array of general purpose bioinformatics programs. For the GEM-PRO pipeline, we mainly need the *needle* pairwise alignment tool (although this can be replaced with Biopython's built-in pairwise alignment function), and the *pepstats* protein sequence statistics tool.


Installation instructions (Ubuntu)
==================================

.. note:: These instructions were created on an Ubuntu 17.04 system.

#. Install the EMBOSS package which contains many programs

    .. code-block:: console
    
        sudo apt-get install emboss

#. And then once that installs, try running the ``needle`` program:

    .. code-block:: console
    
        needle


Installation instructions (Mac OSX, other Unix)
===============================================

#. Just install after downloading the `EMBOSS source code`_
   
    .. code-block:: console

       ./configure
       make
       sudo make install


Program execution
=================

In the shell
------------

To run the program on its own in the shell...
     
    .. code-block:: console

        needle    

With *ssbio*
------------

To run the program using the ssbio Python wrapper, see: 
        
    - :func:`ssbio.protein.sequence.properties.residues.emboss_pepstats_on_fasta`


FAQs
====

* How do I cite EMBOSS?

    - Rice P, Longden I & Bleasby A (2000) EMBOSS: the European Molecular Biology Open Software Suite. Trends Genet. 16: 276–277 Available at: http://www.ncbi.nlm.nih.gov/pubmed/10827456
      
* I'm having issues running EMBOSS programs...

    - See the `ssbio wiki`_ for (hopefully) some solutions - or add yours in when you find the answer!


API
===

.. automodule:: ssbio.protein.sequence.properties.residues
    :members:


.. Links
.. _EMBOSS home page: http://emboss.sourceforge.net/
.. _EMBOSS source code: http://emboss.sourceforge.net/download/#Stable
.. _ssbio wiki: https://github.com/SBRG/ssbio/wiki/Troubleshooting