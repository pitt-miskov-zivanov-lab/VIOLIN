Input and Output Files
======================

This page details the type and format of the files VIOLIN takes as
input, as well as the format of VIOLIN output

Input
-----
As VIOLIN is developed as part of a larger modeling framework, the default files
expected are model spreadsheets in the BioRECIPE format, and a machine reading
output spreadsheet:

.. list-table:: Model input spreadsheet
   :widths: 10 10 10 10 10 10
   :header-rows: 0   
   :class: tight-table   
 
   * - **Element Name**
     - **Element ID**
     - **Element Type**
     - **Variable**
     - **Positive Regulators**
     - **Negative Regulators**
   * - Name
     - ID
     - Type
     - variable name
     - pos regulators list
     - neg regulators list

The minimum required information from a model is: element/target name, element/target type, element/target ID, positive and negative
regulator/source names, types, and IDs. The nomenclature does not have to be specific, as VIOLIN has the capabilities
to recognize all of the commonly used vocabulary.

.. list-table:: LEE input spreadsheet
   :widths: 10 10 10 10 10 10 10 10 10 10 10 
   :header-rows: 0   
   :class: tight-table   
 
   * - **Element Name**
     - **Element ID**
     - **Element Type**
     - **Positive Reg Name**
     - **Positive Reg ID**
     - **Positive Reg Type**
     - **Negative Reg Name**
     - **Negative Reg ID**
     - **Negative Reg Type**
     - **PMCID**
     - **Evidence**
   * - Name
     - ID
     - Type
     - Pos regulator
     - reg ID
     - reg type
     - Neg regulator
     - reg ID
     - reg type
     - PMCID
     - evidence text

The minimum information required from reading output is: element/target name, element/target type, element/target ID, positive and negative
regulator/source names, types, and IDs. The nomenclature does not have to be specific, as VIOLIN has the capabilities
to recognize all of the commonly used vocabulary.


Additionally, VIOLIN has built-in functions to allow for other input file formats,
found on the :doc:`formatting` page.

Currently accepted file types are comma-separated files (**.csv**), 
tab-separated files (**.txt** or **.tsv**), and excel spreadsheets (**.xlsx**).

Output
------

VIOLIN output is sorted into five **.csv** output files:

- Total output - all interactions listed by total score
- Corroborations - all interactions classified as corroborations, listed by total score
- Extensions - all interactions classified as extensions, listed by total score
- Contradictions - all interactions classified as contradictions, listed by total score
- Questionable - all interactions classified as questionable, listed by total score