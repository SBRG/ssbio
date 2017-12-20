.. _opm:

***
OPM
***


Description
===========

* `OPM home page`_
* `OPM web server`_
* `OPM web server instructions and description of results`_

OPM is a program to predict the location of transmembrane planes in protein structures, utilizing the atomic coordinates. *ssbio* provides a wrapper to submit PDB files to the web server, cache, and parse the results


Instructions
============

#. Use the function :meth:`ssbio.protein.structure.properties.opm.run_ppm_server` to upload a PDB file to the PPM server.


FAQs
====

* How can I install OPM?

    * OPM is only available as a web server. *ssbio* provides a wrapper for the web server and allows you to
      submit protein structures to it along with caching the output files.

* How do I cite OPM?

    - Lomize MA, Pogozheva ID, Joo H, Mosberg HI & Lomize AL (2012) OPM database and PPM web server: resources for positioning of proteins in membranes. Nucleic Acids Res. 40: D370–6 Available at: http://dx.doi.org/10.1093/nar/gkr703
      
* I'm having issues running OPM...

    - See the `ssbio wiki`_ for (hopefully) some solutions - or add yours in when you find the answer!


API
===

.. automodule:: ssbio.protein.structure.properties.opm
    :members:


.. Links
.. _OPM home page: http://opm.phar.umich.edu/
.. _OPM web server: http://opm.phar.umich.edu/server.php
.. _OPM web server instructions and description of results: http://sunshine.phar.umich.edu/instruction.htm
.. _ssbio wiki: https://github.com/SBRG/ssbio/wiki/Troubleshooting