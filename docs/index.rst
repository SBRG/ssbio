*************************************************
ssbio: A Framework for Structural Systems Biology
*************************************************

Introduction
============

This Python package provides a collection of tools for people with questions in the realm of structural systems biology. The main goals of this package are to:

#. Provide an easy way to map proteins to sequences and structures
#. Directly link structures to genome-scale SBML models
#. Prepare structures for downstream analyses, such as their use in molecular modeling software

Example questions you can answer with this package:

- How can I determine the number of protein structures available for my list of genes?
- What is the best, representative structure for my protein?
- Where, in a metabolic network, do these proteins work?
- Where do popular mutations show up on a protein?
- How can I compare the structural features of entire proteomes?
- and more...

Installation
============
Clone this repository to any location and then install it.

**Cloning**

.. code-block:: bash

    git clone https://github.com/SBRG/ssbio.git

**Installation**

First install NGLview using pip:

.. code-block:: bash

    pip install nglview

Then install ssbio:

.. code-block:: bash

    cd ssbio
    python setup.py develop --user

**Updating**

.. code-block:: bash

    cd ssbio
    git pull

**Uninstalling**

.. code-block:: bash

    pip uninstall ssbio


Dependencies
------------

See: `Software Installations <https://github.com/SBRG/ssbio/wiki/Software-Installations>`_ for additional programs to install.


Tutorials
=========

Check out some Jupyter notebook tutorials at :ref:`protein` and :ref:`gempro`.


Citation
========

Currently, use of this package can be cited by our 2016 paper in BMC Systems Biology [1]_, which details the GEM-PRO pipeline. The manuscript for the ``ssbio`` package itself is in preparation at this moment.

.. [1] Brunk, E.*, Mih, N.*, Monk, J., Zhang, Z., O’Brien, E. J., Bliven, S. E., Bourne, P. E., Palsson, B. O. (2016). Systems biology of the structural proteome. BMC Systems Biology, 10(1), 26. http://doi.org/10.1186/s12918-016-0271-6. *Authors contributed equally.


Index
========

.. toctree::
   :glob:
   :maxdepth: 3

   getting_started
   structure
   protein
   gempro
   atlas
   python_api

* :ref:`genindex`
