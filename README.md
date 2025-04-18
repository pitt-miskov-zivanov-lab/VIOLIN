# VIOLIN
[![Documentation Status](https://readthedocs.org/projects/melody-violin/badge/?version=latest)](https://melody-violin.readthedocs.io/en/latest/?badge=latest)

## What is VIOLIN?
VIOLIN, stands for **V**ersatile **I**nteraction **O**rganizing to **L**everage **I**nformation in **N**etworks, is a tool written with python used to automatically classify and judge literature-extracted interactions from machine readers or interactions retrieved from database by comparing them to existing curated models. This comparison can help identify key interactions in the context of model.


## Repo Contents

<!-- - [Functionality](#Functionality) -->
- [Installation](#Installation)
- [Tutorial](#Tutorial)
- [Package Structure](#Package-Structure)
- [User Interface](#User-Interface)
- [Citation](#Citation)
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


## Tutorial
This section briefly describes the inputs and outputs of VIOLIN, as well as how to use it via a Jupyter Notebook. All tutorial examples are located in the `examples` folder.
### Input
- `Model` - A model for classification, in BioRECIPE model format. See the example file [modelA_BioRECIPE.xlsx](https://github.com/pitt-miskov-zivanov-lab/VIOLIN/blob/master/examples/input/ModelA_biorecipe.xlsx)
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
- `use_violin_script.py` - [python script](https://github.com/pitt-miskov-zivanov-lab/VIOLIN/blob/master/examples/use_violin_script.py) used to run VIOLIN at the command line



### Customized input
To use your own curated model for VIOLIN, please follow the BioRECIPE model format. Please also note that VIOLIN use unique identifiers for every element in the model to construct the graph. The unique identifiers could distinguish the elements by `name`, `type`, `subtype`, and `compartment ID`, these identifiers are used in regulator list for regarding elements. For example, if IL2 positively influences IL2R, then in `Positive Regulator List` of IL2R, the IL2 should be written as 'IL2_pn_nan_GO0005576', where the `name`, `type`, `subtype`, and `compartment ID` of IL2 is 'IL2', 'protein', '', and 'GO0005576'. Before running the VIOLIN, we encourage you proofread the model file using the command below:
```bash
python preprocessing/normalization.py --format model --file [MODEL_FILENAME]
```

## Package Structure

- [`setup.py`](setup.py): python file that help set up python dependencies installation
- [`src/violin/`](src/violin/): directory that includes core python VIOLIN files
  - [`src/violin/formatting.py`](src/violin/formatting.py): functions for formatting input spreadsheet values.
  - [`src/violin/in_out.py`](src/violin/in_out.py): functions for converting files into readable tables and saving results.
  - [`src/violin/network.py`](src/violin/network.py): functions for building networks and searching paths between nodes.
  - [`src/violin/numeric.py`](src/violin/numeric.py): functions for matching the mutual nodes and attributes in the interaction set and model.
  - [`src/violin/scoring.py`](src/violin/scoring.py): functions for classifying interactions using the decision tree.
  - [`src/violin/visualize_violin.py`](src/violin/visualize_violin.py): functions for visualizing VIOLIN results.

- [`examples/`](examples/): directory that includes tutorial notebook and example inputs and outputs
- [`docs/`](docs/): containing files supporting the repo's host on [Read the Docs](https://melody-violin.readthedocs.io)
- [`LICENSE.txt`](LICENSE.txt): MIT License
- [`README.md`](README.md): this is me!

## Evaluation
### Evaluation on four machine readers and classification schemes
The evaluation and experiment could be found in `eval` folder. `eval/input` includes the results from four readers. `RA` in the filenames stands for the interactions to be compared with the model of Melanoma SkMel133 cell line, and filenames with `RB` prefix are for the model of CD4+ T cell differentiation model. The output from four readers are cached in `input.zip`.
To get the results, run `evaluation_schemes.py` by:
```bash
cd eval
unzip input.zip
python evaluation_schemes.py --scheme [SCHEME] --reader [READER_NAME] --output [OUTPUT_DIRECTORY] --attributes [ATTRIBUTES_LIST]
```

e.g.
```bash
python evaluation_schemes.py --scheme v1 --reader REACH --output ./
```
For more details, check the helper function of argparse or in `evaluation_schemes.py`.
## User Interface
For your convenience, VIOLIN includes a user interface that visualizes the classification results. The interface is built using Next.js for the frontend and FastAPI for the backend.

### Running the frontend
Install npm (if not already installed yet) and start the development server:
```bash
brew install npm
npm run dev
```

### Running the backend
Install FastAPI to start the server:
```bash
pip install fastapi[standard]
fastapi run main.py
```
You could find out the interface at http://localhost:3000, where you can upload your interaction set and model file for classification.

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

## Funding

This work was funded in part by DARPA Big Mechanism award, AIMCancer (W911NF-17-1-0135); and in part by the NSF EAGER award CCF-2324742.

## Support
Feel free to reach out via email nmzivanov@pitt.edu for additional support if you run into any error.
