.. _software:

********
Software
********

.. tabularcolumns:: |\Y{0.2}|\Y{0.4}|\Y{0.4}|\Y{0.4}|\Y{0.2}|\Y{0.4}|\Y{0.4}|\Y{0.4}|

+-------------------+-----------------+------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
| Analysis level    | Function type   | Name                         | Function                                                  | Internal Python class used and functions provided                | External software to install | Web server                   | Alternate external software  |
+===================+=================+==============================+===========================================================+==================================================================+==============================+==============================+==============================+
| Network model or  | Pipeline        | GEM-PRO                      | Pipeline to automatically map gene IDs, protein           | :doc:`gempro`                                                    |                              |                              |                              |
| a set of proteins |                 |                              | sequences, or GEMs to available experimental structures.  |                                                                  |                              |                              |                              |
|                   |                 |                              | Enables streamlined analysis for all functions described  |                                                                  |                              |                              |                              |
|                   |                 |                              | below for individual proteins.                            |                                                                  |                              |                              |                              |
+-------------------+-----------------+------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
| Protein sequence  | Sequence-based  | Various sequence properties  | Basic properties of the sequence, such as percent of      | `Biopython ProteinAnalysis`_,                                    | :doc:`instructions/emboss`   |                              |                              |
|                   | calculation     |                              | polar, non-polar, hydrophobic or hydrophilic residues.    | :mod:`~ssbio.protein.sequence.properties.residues`               |                              |                              |                              |
+                   +                 +------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   |                 | Sequence alignment           | Basic functions to run pairwise or multiple sequence      | `Biopython pairwise2`_,                                          | :doc:`instructions/emboss`   |                              |                              |
|                   |                 |                              | alignments                                                | :mod:`~ssbio.protein.sequence.utils.alignment`                   |                              |                              |                              |
+                   +-----------------+------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   | Sequence-based  | Aggregation propensity       | Consensus method to predict the aggregation propensity of | :mod:`~ssbio.protein.sequence.properties.aggregation_propensity` |                              | :doc:`instructions/amylpred` |                              |
|                   | prediction      |                              | proteins, specifically the number of aggregation-prone    |                                                                  |                              |                              |                              |
|                   |                 |                              | segments on an unfolded protein sequence                  |                                                                  |                              |                              |                              |
+                   +                 +------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   |                 | Secondary structure and      | Predictions of secondary structure and relative solvent   | :mod:`~ssbio.protein.sequence.properties.scratch`                | :doc:`instructions/scratch`  |                              |                              |
|                   |                 | solvent accessibilities      | accessibilities per residue                               |                                                                  |                              |                              |                              |
+                   +                 +------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   |                 | Thermostability              | Free energy of unfolding (ΔG), adapted from Oobatake      | :mod:`~ssbio.protein.sequence.properties.thermostability`        |                              |                              |                              |
|                   |                 |                              | (Oobatake & Ooi 1993) and Dill (Dill et al. 2011)         |                                                                  |                              |                              |                              |
+                   +                 +------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   |                 | Transmembrane domains        | Prediction of transmembrane domains from sequence         | :mod:`~ssbio.protein.sequence.properties.tmhmm`                  | :doc:`instructions/tmhmm`    |                              |                              |
+-------------------+-----------------+------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
| Protein structure | Sequence-based  | Homology modeling            | Preparation scripts and parsers for executing homology    | :mod:`~ssbio.protein.structure.homology.itasser.itasserprep`,    | :doc:`instructions/itasser`  |                              |                              |
|                   | prediction      |                              | modeling algorithms                                       | :mod:`~ssbio.protein.structure.homology.itasser.itasserprop`     |                              |                              |                              |
+                   +-----------------+------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   | Structure-based | Kinetic folding rate         | Prediction of protein folding rates from amino acid       | :mod:`~ssbio.protein.sequence.properties.kinetic_folding_rate`   |                              | :doc:`instructions/foldrate` |                              |
|                   | prediction      |                              | sequence                                                  |                                                                  |                              |                              |                              |
+                   +                 +------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   |                 | Transmembrane orientation    | Prediction of transmembrane domains and orientation in a  | :mod:`~ssbio.protein.structure.properties.opm`                   |                              | :doc:`instructions/opm`      |                              |
|                   |                 |                              | membrane                                                  |                                                                  |                              |                              |                              |
+                   +-----------------+------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   | Structure-based | Secondary structure          | Calculations of secondary structure                       | `Biopython Structure`_,                                          | :doc:`instructions/dssp`     |                              | :doc:`instructions/stride`   |
|                   | calculation     |                              |                                                           | :mod:`~ssbio.protein.structure.properties.dssp`,                 |                              |                              |                              |
|                   |                 |                              |                                                           | :mod:`~ssbio.protein.structure.properties.stride`                |                              |                              |                              |
+                   +                 +------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   |                 | Solvent accessibilities      | Calculations of per-residue absolute and relative solvent | `Biopython Structure`_,                                          | :doc:`instructions/dssp`     |                              | :doc:`instructions/freesasa` |
|                   |                 |                              | accessibilities                                           | :mod:`~ssbio.protein.structure.properties.dssp`,                 |                              |                              |                              |
|                   |                 |                              |                                                           | :mod:`~ssbio.protein.structure.properties.freesasa`              |                              |                              |                              |
+                   +                 +------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   |                 | Residue depths               | Calculations of residue depths                            | `Biopython Structure`_,                                          | :doc:`instructions/msms`     |                              |                              |
|                   |                 |                              |                                                           | :mod:`~ssbio.protein.structure.properties.msms`                  |                              |                              |                              |
+                   +                 +------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   |                 | Structural similarity        | Pairwise calculations of 3D structural similarity         | :mod:`~ssbio.protein.structure.properties.fatcat`                | :doc:`instructions/fatcat`   |                              |                              |
+                   +                 +------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   |                 | Quality                      | Custom functions to allow ranking of structures by        | :func:`~ssbio.core.protein.Protein.set_representative_structure` |                              |                              |                              |
|                   |                 |                              | percent identity to a defined sequence, structure         |                                                                  |                              |                              |                              |
|                   |                 |                              | resolution, and other structure quality metrics           |                                                                  |                              |                              |                              |
+                   +                 +------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   |                 | Various structure properties | Basic properties of the structure, such as distance       | `Biopython Structure`_,                                          |                              |                              |                              |
|                   |                 |                              | measurements between residues or number of disulfide      | :mod:`~ssbio.protein.structure.properties.residues`              |                              |                              |                              |
|                   |                 |                              | bridges                                                   |                                                                  |                              |                              |                              |
+                   +-----------------+------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+
|                   | Structure-based | Structure cleaning, mutating | Custom functions to allow for the preparation of          | `Biopython Structure`_,                                          |                              | AmberTools_                  |                              |
|                   | function        |                              | structure files for molecular modeling, with options to   | :mod:`~ssbio.protein.structure.utils.cleanpdb`,                  |                              |                              |                              |
|                   |                 |                              | remove hydrogens/waters/heteroatoms, select specific      | :mod:`~ssbio.protein.structure.utils.muatatepdb`                 |                              |                              |                              |
|                   |                 |                              | chains, or mutate specific residues.                      |                                                                  |                              |                              |                              |
+-------------------+-----------------+------------------------------+-----------------------------------------------------------+------------------------------------------------------------------+------------------------------+------------------------------+------------------------------+


.. list-table:: Example table
   :header-rows: 1

   * - First header
     - Second header
     - Third header
   * - Some text
     - Some text
     - A list:
         * foo
         * bar
         * baz
   * - Second row
     - More cells
     - etc.
   * - ...
     - ...
     - ...


.. _Biopython Structure: http://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
.. _Biopython ProteinAnalysis: http://biopython.org/wiki/ProtParam
.. _Biopython pairwise2: http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
.. _AmberTools: http://ambermd.org/#AmberTools