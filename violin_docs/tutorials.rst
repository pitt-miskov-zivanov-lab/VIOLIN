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

    <div style="position: relative; padding-bottom: 20.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe width="560" height="315" src="https://www.youtube.com/embed/CIwezpI689c" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
    </div>

Tutorial 2: Specifying Score and Attributes
-------------------------------------------

.. raw:: html

    <div style="position: relative; padding-bottom: 20.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe width="560" height="315" src="https://www.youtube.com/embed/jCSfBH3vzQw" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
    </div>

Tutorial 3: Using VIOLIN at the terminal
----------------------------------------
This tutorial assumes that the input files are a model in the standard BioRECIPES format
and an LEE set from REACH. This tutorial also assumes that the user wants to run VIOLIN
for basic extension using VIOLIN's default values, and visualization is for the total output

The `use_violin_script.py` script is included in the  `violin_tutorial` folder. The input for this script allows for four classification schemes:

#. 'extend' - default Kind and Match Score values for general extension
#. 'extend subcategories' - general extension values with subcategories specified in Kind Score values
#. 'corroborate' - Kind and Match Score values for general corroboration (preference towards strong corroborations, weak corroborations, contradictions)
#. 'corroborate subcategories' - general extension values with subcategories specified in Kind Score values

as well as the same filtering options from :doc:`visualization` 

.. currentmodule:: violin_tutorial.use_violin_script
.. autofunction:: use_violin

To run `use_violin_script.py` at the command line: ::

    python use_violin_script.py test_input/ModelA.csv test_input/RA2_reading.xlsx output extend 50%

.. raw:: html

    <div style="position: relative; padding-bottom: 20.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe width="560" height="315" src="https://www.youtube.com/embed/8kFBzVdhvvg" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
    </div>
    




Tutorial 4: Alternative Input
------------------------------

Machine Reading Output
^^^^^^^^^^^^^^^^^^^^^^

REACH has a specific output format, where the source nodes are separated by their regulation sign - positive regulators are stored 
separately from negative regulators. In other machine readers, source nodes are not separated, and instead the regulation sign is stored
as a node attribute (e.g. increase/decrease or activate/inhibit):

The :py:func:`VIOLIN.formatting.convert_reading` can automatically separate the regulators into the specific positive/negative 
distinctions, along with any associated attributes.

If the user has the following machine reading output spreadsheet:

.. list-table:: LEE input spreadsheet
   :widths: 10 10 10 10 10 10 10 10 10 10
   :header-rows: 0   
   :class: tight-table   
 
   * - **Target Name**
     - **Target ID**
     - **Target Type**
     - **Source Name**
     - **Source ID**
     - **Source Type**
     - **Source Location**
     - **Regulation**
     - **PMCID**
     - **Evidence**
   * - Name
     - ID
     - Type
     - regulator name
     - reg ID
     - reg type
     - reg location
     - increase/decrease
     - PMCID
     - evidence text

Running the following function: ::

    reading_df = convert_reading(reading, 'separate', atts = ['location'])

Would produce the new spreadsheet (stored as a pandas dataframe)

.. list-table:: LEE input spreadsheet
   :widths: 10 10 10 10 10 10 10 10 10 10 10 10 10
   :header-rows: 0   
   :class: tight-table   
 
   * - **Target Name**
     - **Target ID**
     - **Target Type**
     - **Positive Source Name**
     - **Positive Source ID**
     - **Positive Source Type**
     - **Positive Source Location**
     - **Negative Source Name**
     - **Negative Source ID**
     - **Negative Source Type**
     - **Negative Source Location**
     - **PMCID**
     - **Evidence**
   * - Name
     - ID
     - Type
     - Pos regulator
     - reg ID
     - reg type
     - reg location
     - Neg regulator
     - reg ID
     - reg type
     - reg location
     - PMCID
     - evidence text

Model Representation
^^^^^^^^^^^^^^^^^^^^
While the BioRECIPES format is necessary for the DySE modeling framework, models are frequently presented as node-edge lists.
In this case, the model will need to be converted to the BioRECIPE format using the :py:func:`VIOLIN.formatting.nconvert_to_biorecipes` function: ::

    model_df = formatting.convert_to_biorecipes(model, att_list=['location','organism'], separate=False)

In a single step, VIOLIN separates the regulators into positive and negative (if needed) and condenses interactions into lists of regulators
for each element, as shown in :ref:`files:Input`.
