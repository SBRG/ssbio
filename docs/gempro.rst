********************
The GEM-PRO Pipeline
********************

Introduction
============
The GEM-PRO pipeline is focused on annotating genome-scale models with protein structure information. Any SBML model can be used as input to
the pipeline, although it is not required to have a one. Here are the possible starting points for using the pipeline:

* An SBML model in `.xml`, `.mat`, or `.json` formats
* A list of gene IDs (`['b0001', 'b0002', ...]`)
* A dictionary of gene IDs and their sequences (`{'b0001':'MSAVEVEEAP..', 'b0002':'AERAPLS', ...}`)

Tutorials
=========

.. toctree::
   :glob:
   :maxdepth: 2

   notebooks/GEM-PRO*