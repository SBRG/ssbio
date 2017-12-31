.. _fatcat:

*******
FATCAT
*******


Description
===========

* `FATCAT home page`_
* `jFATCAT Java version`_
* `jFATCAT download page`_

FATCAT is a structural alignment tool that allows you to determine the similarity of a pair of protein structures.

.. warning:: Parsing FATCAT results is currently incomplete and will only return TM-scores as of now - but TM-scores only show up in development versions of jFATCAT


Installation instructions
=========================

.. note:: These instructions were created on an Ubuntu 17.04 system.

#. Download the Java port of FATCAT from the `jFATCAT download page`_, under the section "Older file downloads" with the filename "protein-comparison-tool\_<DATE>.tar.gz"
#. Extract it to a place where you store software


Program execution
=================

In the shell
------------

* To run the program on its own in the shell...
   
    .. code-block:: console

        <code>

With *ssbio*
------------

* To run the program using the ssbio Python wrapper, see: :func:`ssbio.protein.structure.properties.fatcat.run_fatcat`. Run it on two structures, pointing to the path of the *runFATCAT.sh* script.


FAQs
====

* How do I cite FATCAT?

    - Ye Y & Godzik A (2003) Flexible structure alignment by chaining aligned fragment pairs allowing twists. Bioinformatics 19 Suppl 2: ii246–55 Available at: https://www.ncbi.nlm.nih.gov/pubmed/14534198
      
* I'm having issues running FATCAT...

    - See the `ssbio wiki`_ for (hopefully) some solutions - or add yours in when you find the answer!


API
===

.. automodule:: ssbio.protein.structure.properties.fatcat
    :members:


.. Links
.. _FATCAT home page: http://fatcat.sanfordburnham.org/
.. _jFATCAT Java version: http://source.rcsb.org/
.. _jFATCAT download page: http://source.rcsb.org/download.jsp
.. _ssbio wiki: https://github.com/SBRG/ssbio/wiki/Troubleshooting