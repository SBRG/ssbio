.. _software:

********
Software
********

This section provides a simple list of external software that may be required to carry out specific computations on a protein sequence or structure. This list only contains software that is wrapped with *ssbio* - there may be other programs that carry out these same functions, and do it better (or worse)!

Tables describing functionalities of these software packages in relation to their input, as well as links to internal wrappers and parses are found on :ref:`sequence` and :ref:`structure` pages.


-------------------------------------------------------------------------------------------------------------------


Protein structure predictions
=============================

Homology modeling
^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    instructions/itasser


Transmembrane orientations
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    instructions/opm
               
Kinetic folding rate
^^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    instructions/foldrate

Protein structure calculations
==============================

Secondary structure
^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    instructions/dssp 
    instructions/stride
              
Solvent accessibilities
^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    instructions/dssp
    instructions/freesasa
              
Residue depths
^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    instructions/msms
              
Structural similarity
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    instructions/fatcat
              
Various structure properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- :mod:`structure residues module <ssbio.protein.structure.properties.residues>`
              
Quality
^^^^^^^

- :func:`set_representative_structure function <ssbio.core.protein.Protein.set_representative_structure>`
              
Structure cleaning, mutating
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- :mod:`cleanpdb module <ssbio.protein.structure.utils.cleanpdb>`
- :mod:`mutatepdb module <ssbio.protein.structure.utils.mutatepdb>`
  

-------------------------------------------------------------------------------------------------------------------

          
Protein sequence predictions
============================

Secondary structure
^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    instructions/scratch
              
Solvent accessibilities
^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    instructions/scratch
              
Thermostability
^^^^^^^^^^^^^^^

- :mod:`thermostability module <ssbio.protein.sequence.properties.thermostability>` 
              
Transmembrane domains
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    instructions/tmhmm
              
Aggregation propensity
^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    instructions/amylpred

Protein sequence calculations
=============================

Various sequence properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    EMBOSS pepstats <instructions/emboss>
              
Sequence alignment
^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    EMBOSS needle <instructions/emboss>