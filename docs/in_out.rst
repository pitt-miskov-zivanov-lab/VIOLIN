Input and Output Functions (:py:mod:`violin.in_out`)
=======================================================

This page details the functions which handle the input files and output of VIOLIN.

For more information on the types of accepted inputs, see :doc:`files`.

Functions
---------

.. currentmodule:: in_out
.. autofunction:: input_biorecipes

.. currentmodule:: in_out
.. autofunction:: input_reading

.. currentmodule:: in_out
.. autofunction:: output


Dependencies
------------
**Python**: `pandas <https://pandas.pydata.org/>`_  and
`NumPy <https://numpy.org/>`_ libraries, and
`os.path <https://docs.python.org/2/library/os.path.html>`_ module

**VIOLIN**: ``formatting`` and ``network`` modules.

Defaults
--------
Default Reading Columns

.. literalinclude:: ../src/violin/in_out.py
    :language: python
    :lines: 34-37
    :lineno-start: 34

Default Model Columns (From BioRECIPES format)

.. literalinclude:: ../src/violin/in_out.py
    :language: python
    :lines: 38-40
    :lineno-start: 38
