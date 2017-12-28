.. _software:

********
Software
********


.. role:: raw-html(raw)
   :format: html


Analysis level: network model or a set of proteins
--------------------------------------------------

+---------------+---------+----------------------------------------------------------+---------------------------------------------------+
| Function type | Name    | Function                                                 | Internal Python class used and functions provided |
+===============+=========+==========================================================+===================================================+
| Pipeline      | GEM-PRO | Pipeline to automatically map gene IDs, protein          | :doc:`gempro`                                     |
|               |         | sequences, or GEMs to available experimental structures. |                                                   |
|               |         | Enables streamlined analysis for all functions described |                                                   |
|               |         | below for individual proteins.                           |                                                   |
+---------------+---------+----------------------------------------------------------+---------------------------------------------------+


Analysis level: protein sequence
--------------------------------

+----------------+------------------------------+-----------------------------------------------------------+--------------------------------------------------------------------------------------------------------+------------------------------+------------------------------+---------------------------------------+
| Function type  | Name                         | Function                                                  | Internal Python class used and functions provided                                                      | External software to install | Web server                   | Alternate external software           |
+================+==============================+===========================================================+========================================================================================================+==============================+==============================+=======================================+
| Sequence-based | Various sequence  properties | Basic properties of the sequence, such aspercent of       | - `Biopython ProteinAnalysis`_  &                                                                      |                              |                              | :doc:`instructions/emboss` *pepstats* |
| calculation    |                              | polar, non-polar, hydrophobic or hydrophilic residues.    | - :mod:`ssbio residues module <~ssbio.protein.sequence.properties.residues>`                           |                              |                              |                                       |
+----------------+------------------------------+-----------------------------------------------------------+--------------------------------------------------------------------------------------------------------+------------------------------+------------------------------+---------------------------------------+
|                | Sequence alignment           | Basic functions to run pairwise or multiple sequence      | `Biopython pairwise2`_  &                                                                              |                              |                              | :doc:`instructions/emboss` *needle*   |
|                |                              | alignments                                                | :mod:`ssbio alignment module <~ssbio.protein.sequence.utils.alignment>`                                |                              |                              |                                       |
+----------------+------------------------------+-----------------------------------------------------------+--------------------------------------------------------------------------------------------------------+------------------------------+------------------------------+---------------------------------------+
| Sequence-based | Aggregation propensity       | Consensus method to predict the aggregation propensity of | :mod:`ssbio aggregation_propensity module <~ssbio.protein.sequence.properties.aggregation_propensity>` |                              | :doc:`instructions/amylpred` |                                       |
| prediction     |                              | proteins, specifically the number of aggregation-prone    |                                                                                                        |                              |                              |                                       |
|                |                              | segments on an unfolded protein sequence                  |                                                                                                        |                              |                              |                                       |
+----------------+------------------------------+-----------------------------------------------------------+--------------------------------------------------------------------------------------------------------+------------------------------+------------------------------+---------------------------------------+
|                | Secondary structure and      | Predictions of secondary structure and relative solvent   | :mod:`ssbio scratch module <~ssbio.protein.sequence.properties.scratch>`                               | :doc:`instructions/scratch`  |                              |                                       |
|                | solvent accessibilities      | accessibilities per residue                               |                                                                                                        |                              |                              |                                       |
+----------------+------------------------------+-----------------------------------------------------------+--------------------------------------------------------------------------------------------------------+------------------------------+------------------------------+---------------------------------------+
|                | Thermostability              | Free energy of unfolding (ΔG), adapted from Oobatake      | :mod:`ssbio thermostability module <~ssbio.protein.sequence.properties.thermostability>`               |                              |                              |                                       |
|                |                              | (Oobatake & Ooi 1993) and Dill (Dill et al. 2011)         |                                                                                                        |                              |                              |                                       |
+----------------+------------------------------+-----------------------------------------------------------+--------------------------------------------------------------------------------------------------------+------------------------------+------------------------------+---------------------------------------+
|                | Transmembrane domains        | Prediction of transmembrane domains from sequence         | :mod:`ssbio tmhmm module <~ssbio.protein.sequence.properties.tmhmm>`                                   | :doc:`instructions/tmhmm`    |                              |                                       |
+----------------+------------------------------+-----------------------------------------------------------+--------------------------------------------------------------------------------------------------------+------------------------------+------------------------------+---------------------------------------+


Analysis level: protein structure
---------------------------------

+-----------------+-------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
| Function type   | Name                          | Function                                                  | Internal Python class used and functions provided                | External software to install | Web server                   | Alternate external software  |
+=================+===============================+===========================================================+==================================================================+==============================+==============================+==============================+
| Sequence-based  | Homology modeling             | Preparation scripts and parsers for executing homology    | - :mod:`~ssbio.protein.structure.homology.itasser.itasserprep`   | :doc:`instructions/itasser`  |                              |                              |
| prediction      |                               | modeling algorithms                                       | - :mod:`~ssbio.protein.structure.homology.itasser.itasserprop`   |                              |                              |                              |
+-----------------+-------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
| Structure-based | Kinetic folding rate          | Prediction of protein folding rates from amino acid       | :mod:`~ssbio.protein.sequence.properties.kinetic_folding_rate`   |                              | :doc:`instructions/foldrate` |                              |
| prediction      |                               | sequence                                                  |                                                                  |                              |                              |                              |
+-----------------+-------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                 | Transmembrane orientation     | Prediction of transmembrane domains and orientation in a  | :mod:`~ssbio.protein.structure.properties.opm`                   |                              | :doc:`instructions/opm`      |                              |
|                 |                               | membrane                                                  |                                                                  |                              |                              |                              |
+-----------------+-------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
| Structure-based | Secondary structure           | Calculations of secondary structure                       | `Biopython Structure`_                                           | :doc:`instructions/dssp`     |                              | :doc:`instructions/stride`   |
| calculation     |                               |                                                           | :mod:`~ssbio.protein.structure.properties.dssp`                  |                              |                              |                              |
|                 |                               |                                                           | :mod:`~ssbio.protein.structure.properties.stride`                |                              |                              |                              |
+-----------------+-------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                 | Solvent accessibilities       | Calculations of per-residue absolute and relative solvent | `Biopython Structure`_                                           | :doc:`instructions/dssp`     |                              | :doc:`instructions/freesasa` |
|                 |                               | accessibilities                                           | :mod:`~ssbio.protein.structure.properties.dssp`                  |                              |                              |                              |
|                 |                               |                                                           | :mod:`~ssbio.protein.structure.properties.freesasa`              |                              |                              |                              |
+-----------------+-------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                 | Residue depths                | Calculations of residue depths                            | `Biopython Structure`_                                           | :doc:`instructions/msms`     |                              |                              |
|                 |                               |                                                           | :mod:`~ssbio.protein.structure.properties.msms`                  |                              |                              |                              |
+-----------------+-------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                 | Structural similarity         | Pairwise calculations of 3D structural similarity         | :mod:`~ssbio.protein.structure.properties.fatcat`                | :doc:`instructions/fatcat`   |                              |                              |
+-----------------+-------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                 | Quality                       | Custom functions to allow ranking of structures by        | :func:`~ssbio.core.protein.Protein.set_representative_structure` |                              |                              |                              |
|                 |                               | percent identity to a defined sequence, structure         |                                                                  |                              |                              |                              |
|                 |                               | resolution, and other structure quality metrics           |                                                                  |                              |                              |                              |
+-----------------+-------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                 | Various structure  properties | Basic properties of the structure, such as distance       | `Biopython Structure`_                                           |                              |                              |                              |
|                 |                               | measurements between residues or number of disulfide      | :mod:`~ssbio.protein.structure.properties.residues`              |                              |                              |                              |
|                 |                               | bridges                                                   |                                                                  |                              |                              |                              |
+-----------------+-------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
| Structure-based | Structure cleaning,  mutating | Custom functions to allow for the preparation of          | `Biopython Structure`_                                           |                              | AmberTools_                  |                              |
| function        |                               | structure files for molecular modeling, with options to   | :mod:`~ssbio.protein.structure.utils.cleanpdb`                   |                              |                              |                              |
|                 |                               | remove hydrogens/waters/heteroatoms, select specific      | :mod:`~ssbio.protein.structure.utils.muatatepdb`                 |                              |                              |                              |
|                 |                               | chains, or mutate specific residues.                      |                                                                  |                              |                              |                              |
+-----------------+-------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+


.. raw:: html
   :file: table_test.html

.. _Biopython Structure: http://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
.. _Biopython ProteinAnalysis: http://biopython.org/wiki/ProtParam
.. _Biopython pairwise2: http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
.. _Biopython DSSP: http://biopython.org/DIST/docs/api/Bio.PDB.DSSP%27-module.html
.. _Biopython ResidueDepth: http://biopython.org/DIST/docs/api/Bio.PDB.ResidueDepth%27-module.html
.. _Biopython Struct: http://biopython.org/wiki/Struct
.. _Biopython Select: http://biopython.org/DIST/docs/api/Bio.PDB.PDBIO%27.Select-class.html
.. _AmberTools: http://ambermd.org/#AmberTools