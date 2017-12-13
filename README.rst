*************************************************
ssbio: A Framework for Structural Systems Biology
*************************************************

Introduction
============

This Python package provides a collection of tools for people with questions in the realm of structural systems biology. The main goals of this package are to:

#. Provide an easy way to map genes to their encoded proteins sequences and structures
#. Directly link structures to genome-scale SBML models
#. Prepare structures for downstream analyses, such as their use in molecular modeling software
#. Demonstrate fully-featured Python scientific analysis environments in Jupyter notebooks

Example questions you can (start to) answer with this package:

- How can I determine the number of protein structures available for my list of genes?
- What is the best, representative structure for my protein?
- Where, in a metabolic network, do these proteins work?
- Where do popular mutations show up on a protein?
- How can I compare the structural features of entire proteomes?
- How can I zoom in and visualize the interactions happening in the cell at the molecular level?
- How do structural properties correlate with my experimental datasets?
- How can I improve the contents of my model with structural data?
- and more...

Installation
============

First install NGLview using pip:

.. code-block:: bash

    pip install nglview

Then install ssbio:

.. code-block:: bash

    pip install ssbio

**Updating**

.. code-block:: bash

    pip install ssbio --upgrade

**Uninstalling**

.. code-block:: bash

    pip uninstall ssbio


Dependencies
------------

See: `Software Installations`_ for additional programs to install. Most of these additional programs are used to predict or calculate properties of proteins.


Tutorials
=========

Check out some Jupyter notebook tutorials for a single Protein_ and or for many in a GEM-PRO_ model.


Citation
========

The manuscript for the ``ssbio`` package can be found and cited at [1]_.

.. [1] Mih N, Brunk E, Chen K, Catoiu E, Sastry A, Kavvas E, et al. ssbio: A Python Framework for Structural Systems Biology. bioRxiv. 2017. p. 165506. doi:10.1101/165506


.. Links
.. _Software Installations: https://github.com/SBRG/ssbio/wiki/Software-Installations
.. _Protein: http://ssbio.readthedocs.io/en/latest/protein.html
.. _GEM-PRO: http://ssbio.readthedocs.io/en/latest/gempro.html