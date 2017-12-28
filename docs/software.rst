.. _software:

********
Software
********

This section provides a simple list of external software that may be required to carry out specific computations on a protein sequence or structure. This list only contains software that is wrapped with *ssbio* -- there may be other programs that carry out these same functions, and do it better (or worse)!

Tables describing functionalities of these software packages in relation to their input, as well as links to internal wrappers and parses are found on :ref:`sequence` and :ref:`structure` pages.


Protein structure predictions
=============================

Homology modeling
-----------------

.. toctree::
    :maxdepth: 1

    instructions/itasser


Transmembrane orientations
--------------------------

.. toctree::
    :maxdepth: 1

    instructions/opm
               
- Kinetic folding rate

    + :doc:`instructions/foldrate`    

Calculations
------------

- Secondary structure

    + :doc:`instructions/dssp` 
    + :doc:`instructions/stride`
              
- Solvent accessibilities

    + :doc:`instructions/dssp`
    + :doc:`instructions/freesasa`
              
- Residue depths

    + :doc:`instructions/msms`
              
- Structural similarity

    + :doc:`instructions/fatcat`
              
- Various structure properties

    + :mod:`structure residues module <ssbio.protein.structure.properties.residues>`
              
- Quality

    + :func:`set_representative_structure function <ssbio.core.protein.Protein.set_representative_structure>`
              
- Structure cleaning, mutating

    + :mod:`cleanpdb module <ssbio.protein.structure.utils.cleanpdb>`
    + :mod:`mutatepdb module <ssbio.protein.structure.utils.mutatepdb>`
          
Protein sequence
================

Predictions
-----------

- Secondary structure

    + :doc:`instructions/scratch`
              
- Solvent accessibilities

    + :doc:`instructions/scratch`
              
- Thermostability

    + :mod:`thermostability module <ssbio.protein.sequence.properties.thermostability>` 
              
- Transmembrane domains

    + :doc:`instructions/tmhmm`
              
- Aggregation propensity

    + :doc:`instructions/amylpred`

Calculations
------------

- Various sequence properties

    + :doc:`instructions/emboss` *pepstats*
              
- Sequence alignment

    + :doc:`instructions/emboss` *needle*