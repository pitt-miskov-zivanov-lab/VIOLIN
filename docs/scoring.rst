.. _scoring:

Scoring (:py:mod:`violin.scoring`)
==================================

This page details the scoring functions of VIOLIN

Match Score
-----------

The Match Score (S\ :sub:`M`\) measures how many new nodes are found in the reading with respect to the model.
For an interaction from the reading **A → B**, where **A** is the regulator and **B** is the
regulated node, this calculation considers 4 cases which determine the scoring outcome:

#. Both **A** and **B** are in the model
#. **A** is in the model, **B** is not
#. **B** is in the model, **A** is not
#. Neither **A** nor **B** are in the model

Deafult Match Level scores are given for the assumption that the user wants to extend a given model without
adding new nodes which may not be useful to the network. Thus, new regulators and new edges between model nodes are
considered most important.

Kind Score
----------

The Kind Score (S\ :sub:`K`\) measures the edges of a reading interaction (LEE) with respect to the model (MI).
The Kind Score easily identifies the classification of an interaction, as well as
searching for paths between nodes in the model when the reading interaction is identified as indirect.
Using the same assumption from the Match Level calculation, the Kind Score represents the following
scenarios:

+----------------------+--------------------------------------------------------+
|    Classification    |                       Definition                       |
+======================+========================================================+
|    Corroboration     |                        LEE matches MI                  |
+----------------------+--------------------------------------------------------+
|      Extension       |       LEE contains information not found in model      |
+----------------------+--------------------------------------------------------+
|     Contradiction    |               LEE disputes information in MI           |
+----------------------+--------------------------------------------------------+
|        Flagged       |                 Must be judged manually                |
+----------------------+--------------------------------------------------------+

And within each classification, there are further sub-classifications.
These subclassifications allow for more detailed scoring, if the user wishes.

Corroborations
^^^^^^^^^^^^^^
    Strong Corroboration: LEE matches MI exactly

    Weak Corroboration Type 1: LEE matches direction, sign, connection type, and node type, of a model interaction
    but is missing additional attributes

    Weak Corroboration Type 2: an indirect LEE matches direction and sign of direct model interaction
    with non-contradictory attributes

    Weak Corroboration Type 3: an indrect LEE matches the direction and sign of a *path* in the model
    with non-contradictory attributes

Extensions
^^^^^^^^^^
    Full Extension: Neither source nor target of the LEE is in the model

    Hanging Extension: The target of the LEE is in the model

    Internal Extension: Both the source and target of the LEE are in the model,
    but there is no model interaction between them

    Specification: LEE contains more information (attributes) than MI, or
    shows a direct relationship compared to Model Path

Contradictions
^^^^^^^^^^^^^^
    Direction Contradiction: The target and source of the LEE correspond to
    the source and target of the model interaction, respectively

    Sign Contradiction: The regulation sign of the LEE is opposite of the corresponding
    model interaction (e.g. the LEE shows a positive regulation where the model interaction shows negative)

    Attribute Contradiction: One or more of the LEE node attributes differs from that found
    in the corresponding model interaction

Flagged
^^^^^^^
    Flagged Type 1: Mismatched Direction and non-contradictory Other
    Attributes with a Direct connection type in the model

    Flagged Type 2: An LEE with a corresponding path which has one or
    more Mismatched Attributes

    Flagged Type 3: An LEE which is a self-regulation based on the definition
    of model element
    (e.g. LEE has caspase-8 --> caspase-3, but the model considers cas-8 and cas-3 to be the same element)



Evidence Score
--------------
The Evidence Score (S\ :sub:`E`\) is a measure of how many times an LEE is found in the machine reading output. In the :py:func:`violin.formatting.evidence_score` function, column names
are defined to determine how the function determines duplicates. For example, the Evidence Score can be calculated by comparing all LEE attributes and all machine readings spreadsheet columns.
So only an exact match between LEEs will be counted as a duplicate. However, the user can also define fewer attributes, creating a more coarse-grained Evidence Score calculation.

Epistemic Value
---------------

In the NLP output, we sometimes receive an Epistemic Value (S\ :sub:`B`\), which is a measure
of the believability of an interaction in the LEI. Zero, Low, Moderate, and High
believability correspond to numerical scores of 0.0, 0.33, 0.67, and 1.0, respectively.

Total Score
-----------
The total score (S\ :sub:`T`\) is calculated by

.. math:: S_T = [S_K + (S_E*S_M)]*S_B


Functions
---------

.. currentmodule:: scoring
.. autofunction:: match_score

.. currentmodule:: scoring
.. autofunction:: kind_score

.. currentmodule:: scoring
.. autofunction:: epistemic_value

.. currentmodule:: scoring
.. autofunction:: score_reading


Dependencies
------------
**Python**: `pandas <https://pandas.pydata.org/>`_ library

**VIOLIN**: ``network`` and ``numeric`` modules.

Defaults
--------
Default Match Score values

.. literalinclude:: ../src/violin/scoring.py
    :language: python
    :lines: 28-31
    :lineno-start: 28

Default Kind Score values

.. literalinclude:: ../src/violin/scoring.py
    :language: python
    :lines: 14-27
    :lineno-start: 14

Usage
-----
*scoring.score_reading* scores the reading output in the following manner:


.. literalinclude:: ../src/violin/scoring.py
    :language: python
    :lines: 406-414
    :lineno-start: 406
