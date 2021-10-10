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

A `violin_tutorial` folder containing an example ipy notebook, input files, and expected output
is included with the VIOLIN package. This tutorial recreates some of the data presented in 
[Paper TBA]

Tutorial 2: Specifying Score and Attributes
-------------------------------------------
If the user wants to compare attributes in addition to the default, the additional attributes must be defined
before inputting the reading, and then call them in the input_reading function::

    attributes = ['Location ID']
    reading_df, reading_cols = input_reading(reading_file, evidence_score_cols,atts=attributes)

If the user wants to use a different scoring scheme than the default (e.g. to account for the subclassifications in judgement),
a scoring dictionary must be defined before calling the scoring function::

    kind_dict = {"strong corroboration" : 2, 
                "weak corroboration1" : 1,
                "weak corroboration2" : 3,
                "weak corroboration3" : 5,
                "hanging extension" : 40, 
                "full extension" : 41, 
                "internal extension" : 42, 
                "specification" : 30, 
                "dir contradiction" : 10,
                "sign contradiction" : 11,
                "att contradiction" : 12,
                "flagged1" : 20,
                "flagged2" : 21,
                "flagged3" : 22}
    scored = score_reading(reading_df,model_df,graph,reading_cols,kind_values = kind_dict,attributes = attributes)

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
.. autofunction:: violin_tutorial.use_violin_script

To run `use_violin_script.py` at the command line: ::

    python use_violin_script.py test_input/ModelA.csv test_input/RA2_reading.xlsx output extend 50%




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
