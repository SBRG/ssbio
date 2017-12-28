.. _freesasa:

**********************************
FreeSASA
**********************************

.. image:: ../assets/ssbio_protein_props_sasa.png
    :align: center
    :alt: SASA
    :scale: 30 %


Description
===========

* `FreeSASA home page`_
* `FreeSASA Github`_

FreeSASA is an open source library written in C for calculating solvent accessible surface areas of a protein. FreeSASA also contains Python bidings, and the plan is to include these bindings with *ssbio* in the future.


Instructions (Unix)
===================

.. note:: These instructions were created on an Ubuntu 17.04 system with a Python installation through Anaconda3.

.. note:: FreeSASA Python bindings are slightly difficult to install with Python 3 - ssbio provides wrappers for the command line executable instead


#. Download the latest tarball (see `FreeSASA home page`_), expand it and run
    
    .. code-block:: console

        ./configure --enable-python-bindings CFLAGS="-fPIC -O2"
        make

#. If you have a user-specific Python executable (ie. through Anaconda), edit the freesasa-2.0/bindings/Makefile, lines 805, 809, 815 to change:

    .. code-block:: console
        
        python setup.py [...]

    to (type `which python` to get the path to enter): 

    .. code-block:: console

        /path/to/your/anaconda/python setup.py [...]

#. Install with

    .. code-block:: console
    
        sudo make install


FAQs
====

* How do I cite FreeSASA?

    - Mitternacht S (2016) FreeSASA: An open source C library for solvent accessible surface area calculations. F1000Res. 5: 189 Available at: http://dx.doi.org/10.12688/f1000research.7931.1

* I'm having issues running FreeSASA...

    - See the `ssbio wiki`_ for (hopefully) some solutions - or add yours in when you find the answer!


API
===

.. automodule:: ssbio.protein.structure.properties.freesasa
    :members:


.. Links
.. _FreeSASA home page: http://freesasa.github.io/
.. _FreeSASA Github: https://github.com/mittinatten/freesasa
.. _ssbio wiki: https://github.com/SBRG/ssbio/wiki/Troubleshooting