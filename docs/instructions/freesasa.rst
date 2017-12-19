.. _freesasa:

**********************************
freesasa
**********************************

.. image:: ../assets/ssbioStructPropssasa.png
    :align: center
    :alt: SASA
    :scale: 60 %


Description
===========

* `freesasa home page`_
* `freesasa Github`_

freesasa 


Instructions
============

.. note:: These instructions were created on an Ubuntu 17.04 system with a Python installation through Anaconda3.

.. note:: freesasa Python bindings do not work with Python 3 - ssbio provides wrappers for the command line executable


#. Download the latest tarball (see home page), expand it and run
    
    .. code-block:: console

        ./configure --enable-python-bindings CFLAGS="-fPIC -O2"
        make

#. Edit the freesasa-2.0/bindings/Makefile, lines 805, 809, 815 to change

    .. code-block:: console
        
        python setup.py [...]

    to

    .. code-block:: console

        /path/to/anaconda/python setup.py [...]

#. Install with

    .. code-block:: console
    
        sudo make install


FAQs
====

* How do I cite freesasa?

    - Mitternacht S (2016) FreeSASA: An open source C library for solvent accessible surface area calculations. F1000Res. 5: 189 Available at: http://dx.doi.org/10.12688/f1000research.7931.1

* I'm having issues running freesasa...

    - See the `ssbio wiki`_ for (hopefully) some solutions - or add yours in when you find the answer!


API
===

.. automodule:: ssbio.protein.structure.properties.freesasa
    :members:


.. Links
.. _freesasa home page: http://freesasa.github.io/
.. _freesasa Github: https://github.com/mittinatten/freesasa
.. _ssbio wiki: https://github.com/SBRG/ssbio/wiki/Troubleshooting