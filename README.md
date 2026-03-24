# VIOLIN
<!-- [![Documentation Status](https://readthedocs.org/projects/melody-violin/badge/?version=latest)](https://melody-violin.readthedocs.io/en/latest/?badge=latest) -->
[![API Status](https://img.shields.io/website?url=http%3A%2F%2Fboheme-alb-616936326.us-east-1.elb.amazonaws.com%2Fdocs&up_message=live&up_color=FFEDD4&down_message=inactive&down_color=red&label=web)](http://boheme-alb-616936326.us-east-1.elb.amazonaws.com/)
[![Documentation Status](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pitt-miskov-zivanov-lab/VIOLIN/HEAD?labpath=%2Fexamples%2Fuse_VIOLIN.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/pitt-miskov-zivanov-lab/VIOLIN/blob/master/LICENSE)


## Introduction to VIOLIN
VIOLIN (**V**ersatile **I**nteraction **O**rganizing to **L**everage **I**nformation in **N**etworks) is a configurable, attribute-aware reconciliation framework that formally compares new interaction lists against structured baseline graphs, enabling systematic knowledge integration. VIOLIN classifies each interaction as a corroboration, contradiction, flagged case, or extension, and supports configurable attribute inclusion strategies and mismatch semantics to adjust reconciliation strictness. 

As automated extraction technologies continue to expand the volume of structured interaction data, reconciliation frameworks such as VIOLIN will play a critical role in maintaining coherence between evolving literature and curated knowledge representations. The ability to quantify corroboration, contradiction, ambiguity, and expansion in a transparent and configurable manner provides a foundation for reproducible and scalable model updating in AI-assisted scientific discovery.

## Note
This repository includes the VIOLIN stand-alone package, a collection of files that can be used as inputs to VIOLIN, Jupyter notebooks, and VIOLIN documentation. The web-based user-interface is available [here](http://boheme-alb-616936326.us-east-1.elb.amazonaws.com/violin) and the latest versions can be found by switching branch to [webdev](https://github.com/pitt-miskov-zivanov-lab/VIOLIN/tree/webdev).

## Repo Contents

<!-- - [Functionality](#Functionality) -->
- [Installation](#Installation)
- [Tutorial](#Tutorial)
- [Repository structure](#Repository%20structure)
- [User Interface](#User-Interface)
- [Citation](#Citation)
- [Acknowledgments](#Acknowledgments)
- [Funding](#Funding)
- [Support](#Support)

<!-- ## Functionality
- Identification: identifying the importance of the interactions for a user-defined cellular model
- Error Detection: finding the contradictions in interactions list
- Scoring: showing the quality of the interactions for a model -->

## Installation

1. Clone the VIOLIN repository to your computer.
   ```
   git clone https://github.com/pitt-miskov-zivanov-lab/VIOLIN.git
   ```
2. Navigate into the directory, install VIOLIN and its python dependencies.
   ```
   cd VIOLIN
   pip install -e .
   ```
3. Run the provided notebook (Check [Jupyter notebook installation](https://jupyter.org/install) here).
   ```
   jupyter notebook examples/use_VIOLIN.ipynb
   ```
4. (Optional) You can also interact with VIOLIN with command line:
   ```bash
   cd examples
   # For arguments info:
   python run_VIOLIN.py -h
   # Example:
   python run_VIOLIN.py --model ./input/models/SkMel133_biorecipe.xlsx --reading ./input/interactions/RA2_reading_BioRECIPE.xlsx --output ./ex --summary pie_plots --scheme '1' --show_plot
   ```

## Tutorial
This section briefly describes the inputs and outputs of VIOLIN, as well as how to use it via a Jupyter Notebook. All tutorial examples are located in the `examples` folder. The Binder notebook also provides an online runtime environment to run VIOLIN. 


VIOLIN I/O follows [BioRECIPE](https://melody-biorecipe.readthedocs.io/en/latest/) format standards to structure the inputs and outputs:
### Input
- `Model` - A model for classification, in BioRECIPE model format. See the example file [SkMel133_biorecipe.xlsx](https://github.com/pitt-miskov-zivanov-lab/VIOLIN/blob/master/examples/input/SkMel133_biorecipe.xlsx)
- `interaction set` - A interaction set, in BioRECIPE interaction list format. An example file is provided in [RA2_reading_BioRECIPE.xlsx](https://github.com/pitt-miskov-zivanov-lab/VIOLIN/blob/master/examples/input/interactions/RA2_reading_BioRECIPE.xlsx).

### Output
The output for each pair of a interaction set and a model consists of 5 files, which are:

| Output Annotations    | Contents                                             |
|-----------------------|------------------------------------------------------|
| #_corroborations.csv | All interactions classified as collaborations        |
| #_contradictions.csv | All interactions classified as contradictions        |
| #_extensions.csv     | All interactions classified as extensions            |
| #_flagged.csv        | All interactions classified as flagged               |
| #_outputDF.csv       | All interactions                                     |

### Usage
The example usage for VIOLIN could be found in Jupyter notebook and the python script:
- `use_VIOLIN.ipynb` - [Jupyter notebook](https://github.com/pitt-miskov-zivanov-lab/VIOLIN/blob/master/examples/use_VIOLIN.ipynb) used to run VIOLIN
- `run_VIOLIN.py` - [python script](https://github.com/pitt-miskov-zivanov-lab/VIOLIN/blob/master/examples/run_VIOLIN.py) used to run VIOLIN at the command line


## Repository structure 

- [`setup.py`](setup.py): python file that help set up python dependencies installation
- [`src/violin/`](src/violin/): directory that includes core python VIOLIN files
  - [`src/violin/formatting.py`](src/violin/formatting.py): functions for formatting input spreadsheet values.
  - [`src/violin/in_out.py`](src/violin/in_out.py): functions for converting files into readable tables and saving results.
  - [`src/violin/network.py`](src/violin/network.py): functions for building networks and searching paths between nodes.
  - [`src/violin/numeric.py`](src/violin/numeric.py): functions for matching the mutual nodes and attributes in the interaction set and model.
  - [`src/violin/scoring.py`](src/violin/scoring.py): functions for classifying interactions using the decision tree.
  - [`src/violin/visualize_violin.py`](src/violin/visualize_violin.py): functions for visualizing VIOLIN results.

- [`examples/`](examples/): directory that includes tutorial notebook and example inputs and outputs.
<!-- - [`docs/`](docs/): containing files supporting the repo's host on [Read the Docs](https://melody-violin.readthedocs.io) -->
- [`eval/`](eval/): directory that includes all testing files and scripts. 
- [`LICENSE.txt`](LICENSE.txt): MIT License
- [`README.md`](README.md): this is me!

## Evaluation
### Evaluation on four machine readers and classification schemes
The evaluation and experiment could be found in `eval` folder. `eval/input` includes the results from four readers. `RA` in the filenames stands for the interactions to be compared with the model of Melanoma SkMel133 cell line, and filenames with `RB` prefix are for the CD4+ T cell differentiation model. The output from four readers and their filtered results are included in input folder.
To get the results, run `evaluation_schemes.py` by:
```bash
cd eval
python evaluation_schemes.py --scheme <SCHEME> --reader <READER_NAME> --output <OUTPUT_DIRECTORY> --attributes <ATTRIBUTES_LIST>
```

e.g.:
```bash
python evaluation_schemes.py --scheme v1 --reader REACH --output ./
```
For more details, check the helper function of argparse or in `evaluation_schemes.py`.
`manual_curation` folder contains manually curated results (referred to as Group 5 and Group 6 in the paper). To evaluate classification accuracy on this data, run `mc_eval.py` by: 
```bash
python mc_eval.py --model <baseline_graph_filename.xlsx> --reading <curated_interaction_list.xlsx> --output <evaluation_results.json>
```

## User Interface
Graph structure is complicated to view with a spreadsheet, VIOLIN includes a user interface that visualizes the classification results. The interface is built using React for the frontend and FastAPI for the backend. To access the user interface, first access the [VIOLIN UI](http://boheme-alb-616936326.us-east-1.elb.amazonaws.com/violin) page and upload the interaction list and model files in corresponding windows. 


## Citation

If you use this code or find this effort useful or relevant to your research, please cite the following paper:
```
@article{luo2024context,
  title={Context-driven interaction retrieval and classification for modeling, curation, and reuse},
  author={Luo, Haomiao and Hansen, Casey and Telmer, Cheryl A and Tang, Difei and Arazkhani, Niloofar and Zhou, Gaoxiang and Spirtes, Peter and Miskov-Zivanov, Natasa},
  journal={bioRxiv},
  pages={2024--07},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Acknowledgments
The original version of interaction classifier was developed by Dr. Peter Spirtes (CMU) as part of the AIMCancer project (DARPA's Big Mechanism program) and later revised and further expanded by Dr. Casey Hansen during her PhD studies in Melody Lab (Pitt). Dr. Cheryl Telmer (CMU) has provided invaluable feedback and biological, as well as manually curated files for testing VIOLIN throughout the project. Dr. Gaoxiang Zhou (Pitt) assisted in maintaining the code and developed the Binder notebook. The current version of VIOLIN has been developed primarily by Haomiao Luo at MeLoDy Lab. The front-end for the web-based interface was developed by Niloofar Arazkhani and the back-end was set up on Amazon AWS by Difei Tang at MeLoDy Lab. 

## Funding
This work was funded by DARPA Big Mechanism award W911NF-17-1-0135, NSF EAGER award CCF-2324742, and NIH award R01LM014673.

## Support
Feel free to reach out via email nmzivanov@pitt.edu for additional support if you run into any error.
