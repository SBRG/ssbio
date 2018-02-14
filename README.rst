*************************************************
ssbio: A Framework for Structural Systems Biology
*************************************************


Introduction
============

This Python package provides a collection of tools for people with questions in the realm of structural systems biology. The main goals of this package are to:

#. Provide an easy way to map hundreds or thousands of genes to their encoded protein sequences and structures
#. Directly link protein structures to genome-scale metabolic models
#. Demonstrate fully-featured Python scientific analysis environments in Jupyter notebooks

Example questions you can (start to) answer with this package:

- How can I determine the number of protein structures available for my list of genes?
- What is the best, representative structure for my protein?
- Where, in a metabolic network, do these proteins work?
- Where do popular mutations show up on a protein?
- How can I compare the structural features of entire proteomes?
- How do structural properties correlate with my experimental datasets?
- How can I improve the contents of my metabolic model with structural data?


Try it without installing
=========================

.. image:: https://mybinder.org/badge.svg
    :target: https://mybinder.org/v2/gh/SBRG/ssbio/master?filepath=Binder.ipynb


.. note:: Binder notebooks are still in beta, but they mostly work! Third-party programs are also preinstalled in the Binder notebooks except I-TASSER and TMHMM due to licensing restrictions.


Installation
============

First install NGLview using pip, then install ssbio

.. code-block:: console

    pip install nglview
    jupyter-nbextension enable nglview --py --sys-prefix
    pip install ssbio

Updating
--------

.. code-block:: console

    pip install ssbio --upgrade

Uninstalling
------------

.. code-block:: console

    pip uninstall ssbio

Dependencies
------------

See: Software_ for a list of external programs to install, along with the functionality that they add. Most of these additional programs are used to predict or calculate properties of proteins, and are only required if you desire to calculate the described properties.


Tutorials
=========

Check out some Jupyter notebook tutorials for a single Protein_ and or for many in a GEM-PRO_ model. See a list of all Tutorials_.


Citation
========

The manuscript for the *ssbio* package can be found and cited at [1]_.

.. [1] Mih N, Brunk E, Chen K, Catoiu E, Sastry A, Kavvas E, Monk JM, Zhang Z, Palsson BO. 2018. ssbio: A Python Framework for Structural Systems Biology. Bioinformatics. https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty077/4850940.

.. Links
.. _Protein: http://ssbio.readthedocs.io/en/latest/protein.html
.. _GEM-PRO: http://ssbio.readthedocs.io/en/latest/gempro.html
.. _Software: http://ssbio.readthedocs.io/en/latest/software.html
.. _Tutorials: http://ssbio.readthedocs.io/en/latest/tutorials.html