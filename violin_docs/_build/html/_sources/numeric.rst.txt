Numerical Functions (:py:mod:`VIOLIN.numeric`)
==============================================

This page details the numeric operators of VIOLIN.

These functions discretize qualitative operations within VIOLIN:

#. searching for an element in the reading output,
#. comparing attributes, indentifying whether a given attribute

   * matches exactly an attribute in a corresponding model interaction (MI),
   * is missing where a MI attribute is present,
   * is present where a MI attribute is missing,
   * differs from an attribute in a corresponding MI.

Both functions return numerical values to represent the outcome of the function.


Functions
---------

.. currentmodule:: VIOLIN.numeric
.. autofunction:: find_element

.. currentmodule:: VIOLIN.numeric
.. autofunction:: compare


Dependencies
------------
**Python**: `pandas <https://pandas.pydata.org/>`_ library

**VIOLIN**: none

Usage
-----


This example searches for the index of protein family AMPK in the model spreadsheet.::

    find_element_index("name","ampk","protein family",model_df)
    >> 2

This example compares the *location* of the reading interaction to the *location* of its counterpart interaction in the model.::
    
    reading_att = "nan"
    model_att = "GO:0005737"
    compare(model_att,reading_att)
    >> 2
