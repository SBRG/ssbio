.. _amylpred:

**********************************
AMYLPRED2 - Aggregation Propensity
**********************************

Description
===========

Home page: AMYLPRED2_

This module provides a function to predict the aggregation propensity of proteins, specifically the number
of aggregation-prone segments on an unfolded protein sequence. AMYLPRED2 is a consensus method of different methods.
In order to obtain the best balance between sensitivity and specificity, we follow the author's guidelines to consider
every 5 consecutive residues agreed among at least 5 methods contributing 1 to the aggregation propensity.


Instructions
============

#. Create an account on the webserver at the `AMYLPRED2 registration link`_.
#. Create a new AMYLPRED object with your email and password initialized along with it.
#. Run ``get_aggregation_propensity`` on a protein sequence.


FAQs
====

* What is aggregation propensity?

    * The number of aggregation-prone segments on an unfolded protein sequence.

* How can I install AMYLPRED2?

    * AMYLPRED2 is only available as a web server. *ssbio* provides a wrapper for the web server and allows you to
      submit protein sequences to it along with caching the output files.

* How do I cite AMYLPRED2?

    * Tsolis AC, Papandreou NC, Iconomidou VA & Hamodrakas SJ (2013) A consensus method for the prediction of
      'aggregation-prone' peptides in globular proteins. PLoS One 8: e54175 Available at:
      http://dx.doi.org/10.1371/journal.pone.0054175

* How can this parameter be used on a genome-scale?

    * See: Chen K, Gao Y, Mih N, O'Brien EJ, Yang L & Palsson BO (2017) Thermosensitivity of growth is determined by
      chaperone-mediated proteome reallocation. Proceedings of the National Academy of Sciences 114: 11548–11553
      Available at: http://www.pnas.org/content/114/43/11548.abstract


API
===

.. automodule:: ssbio.protein.sequence.properties.aggregation_propensity
    :members:


.. Links
.. _AMYLPRED2: http://aias.biol.uoa.gr/AMYLPRED2/
.. _AMYLPRED2 registration link: http://aias.biol.uoa.gr/AMYLPRED2/register.php
