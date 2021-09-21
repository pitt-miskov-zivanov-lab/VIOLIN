Formatting (:py:mod:`VIOLIN.formatting`)
========================================

This page details the formatting functions of VIOLIN, used during model and reading input. 

The formatting step is important, as it:

- identifies duplicate interactions in the reading output,
- counts the number of times an interaction was found in the reading (:ref:`scoring:Evidence Score`),
- converts the variable representation of the model regulators into the common names

The formatting functions are also responsible for inputting models and machine
reading output which are not in the BioRECIPES or REACH format (respectively).

Functions
---------

.. currentmodule:: VIOLIN.formatting
.. autofunction:: evidence_score

.. currentmodule:: VIOLIN.formatting
.. autofunction:: add_regulator_names_id

.. currentmodule:: VIOLIN.formatting
.. autofunction:: convert_to_biorecipes

.. currentmodule:: VIOLIN.formatting
.. autofunction:: convert_reading


Dependencies
------------
**Python**: `pandas <https://pandas.pydata.org/>`_ and 
`NumPy <https://numpy.org/>`_ libraries, as well as the
`os.path <https://docs.python.org/2/library/os.path.html>`_ module

**VIOLIN**: none

Usage
-----
This module is used in during file input in the `input/output` module. For an example of using the `convert` functions, see :ref:`tutorials:Tutorial 4: Alternative Input`.

