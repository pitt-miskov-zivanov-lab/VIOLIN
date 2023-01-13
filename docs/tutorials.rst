#########
Tutorials
#########

Test files
----------
Input files used to create publication data can be found at [link]

Tutorial 1: Default Use (general extension)
-------------------------------------------
This tutorial assumes that the input files are a model in the standard BioRECIPES format
and an LEE set from REACH.

.. raw:: html

    <div style="position: relative; padding-bottom: 10.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe width="560" height="315" src="https://www.youtube.com/embed/CIwezpI689c" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
    </div>

Tutorial 2: Specifying Score and Attributes
-------------------------------------------

.. raw:: html

    <div style="position: relative; padding-bottom: 10.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe width="560" height="315" src="https://www.youtube.com/embed/jCSfBH3vzQw" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
    </div>

Tutorial 3: Using VIOLIN at the terminal
----------------------------------------
.. raw:: html

    <div style="position: relative; padding-bottom: 10.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe width="560" height="315" src="https://www.youtube.com/embed/8kFBzVdhvvg" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
    </div>

This tutorial assumes that the input files are a model in the standard BioRECIPES format
and an LEE set from REACH. This tutorial also assumes that the user wants to run VIOLIN
for basic extension using VIOLIN's default values, and visualization is for the total output

The `use_violin_script.py` script is included in the  `violin_tutorial` folder. The input for this script allows for four classification schemes:

#. 'extend' - default Kind and Match Score values for general extension
#. 'extend subcategories' - general extension values with subcategories specified in Kind Score values
#. 'corroborate' - Kind and Match Score values for general corroboration (preference towards strong corroborations, weak corroborations, contradictions)
#. 'corroborate subcategories' - general extension values with subcategories specified in Kind Score values

as well as the same filtering options from :doc:`visualization`

.. currentmodule:: use_violin_script
.. autofunction:: use_violin

To run `use_violin_script.py` at the command line: ::

    python use_violin_script.py test_input/ModelA.csv test_input/RA2_reading.xlsx output extend 50%



Tutorial 4: Alternative Input
------------------------------
.. raw:: html

    <div style="position: relative; padding-bottom: 10.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe width="560" height="315" src="https://www.youtube.com/embed/gkHnisreKWo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
    </div>
