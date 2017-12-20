.. _foldrate:

*********
FOLD-RATE
*********

Description
===========

* `FOLD-RATE home page`_

This module provides a function to predict the **kinetic folding rate** (k\ :sup:`f`) given an amino acid sequence and its structural classficiation (alpha/beta/mixed).


Instructions
============

#. Obtain your protein's sequence
#. Determine the main secondary structure composition of the protein (``all-alpha``, ``all-beta``, ``mixed``, or ``unknown``)
#. Input the sequence and secondary structure composition into the function ``get_foldrate``


FAQs
====

* What is the main secondary structure composition of my protein?

    * ``all-alpha`` = dominated by α-helices; α > 40% and β < 5%
    * ``all-beta`` = dominated by β-strands; β > 40% and α < 5%
    * ``mixed`` = contain both α-helices and β-strands; α > 15% and β > 10%

* What is the kinetic folding rate?

    * Protein folding rate is a measure of slow/fast folding of proteins from the unfolded state to native
      three-dimensional structure.

* What units is it in?

    * Number of proteins folded per second

* How can I install FOLD-RATE?

    * FOLD-RATE is only available as a web server. *ssbio* provides a wrapper for the web server and allows you to
      submit protein sequences to it along with caching the output files.

* How do I cite FOLD-RATE?

    * Gromiha MM, Thangakani AM & Selvaraj S (2006) FOLD-RATE: prediction of protein folding rates from amino acid
      sequence. Nucleic Acids Res. 34: W70–4 Available at: http://dx.doi.org/10.1093/nar/gkl043

* How can this parameter be used on a genome-scale?

    * See: Chen K, Gao Y, Mih N, O'Brien EJ, Yang L & Palsson BO (2017) Thermosensitivity of growth is determined by
      chaperone-mediated proteome reallocation. Proceedings of the National Academy of Sciences 114: 11548–11553
      Available at: http://www.pnas.org/content/114/43/11548.abstract


API
---
.. automodule:: ssbio.protein.sequence.properties.kinetic_folding_rate
    :members:


.. Links
.. _FOLD-RATE home page: http://www.iitm.ac.in/bioinfo/fold-rate/
