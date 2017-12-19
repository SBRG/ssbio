.. _stride:

**********************************
STRIDE
**********************************

.. image:: ../assets/ssbioStructPropssecstruct.png
    :align: center
    :alt: Secondary structure
    :scale: 60 %


Description
===========

* `STRIDE home page`_
* `STRIDE download page`_

STRIDE is...


Instructions
============

.. note:: These instructions were created on an Ubuntu 17.04 system.

* `Instructions for installing on Mac`_
* `Instructions for installing on Mac (alternate)`_

#. Download the source from the `STRIDE download page`_

#. Build the program from source:

    .. code-block:: console

        mkdir stride
        cp stride.tar.gz stride
        cd stride
        tar -zxf stride.tar.gz
        make
        cp stride /usr/local/bin

#. Then you should be able to run ``stride`` in your terminal


FAQs
====

* How do I cite STRIDE?

    - Frishman D & Argos P (1995) Knowledge-based protein secondary structure assignment. Proteins 23: 566–579 Available at: http://dx.doi.org/10.1002/prot.340230412


API
===

.. automodule:: ssbio.protein.structure.properties.stride
    :members:


.. Links
.. _STRIDE home page: http://webclu.bio.wzw.tum.de/stride/
.. _STRIDE download page: http://webclu.bio.wzw.tum.de/stride/install.html
