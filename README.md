# VIOLIN
[![Documentation Status](https://readthedocs.org/projects/melody-violin/badge/?version=latest)](https://melody-violin.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pitt-miskov-zivanov-lab/VIOLIN/HEAD?labpath=%2Fexamples%2Fuse_VIOLIN.ipynb)

### VIOLIN (Validating Interactions Of Likely Importance to the Network)

VIOLIN (Validating Interactions Of Likely Importance to the Network) is a tool used to automatically classify and judge literature-extracted interactions curated from machine readers by comparing them to existing models. This comparison can help identify key interactions for model extension.

## Contents

- [Functionality](#Functionality)
- [I/O](#IO)
- [Online Tutorial](#Online-Tutorial)
- [Offline Installation](#Offline-Installation)
- [Package Structure](#Package-Structure)
- [Reproducibility](#Reproducibility)
- [Citation](#Citation)
- [Funding](#Funding)
- [Support](#Support)

## Functionality
- Identification: identifying the importance of the interactions for a user-defined cellular model
- Error Detection: finding the contradictions in interactions list
- Scoring: showing the quality of the interactions for a model

## I/O

### Input/Output
- `exmaple/input` - [folder](https://github.com/pitt-miskov-zivanov-lab/VIOLIN/blob/master/examples/input) containing model files and automatically retrieved interactions lists.
- `example/output` - [example/output/<interactions_filename>_<class>.csv](https://github.com/pitt-miskov-zivanov-lab/VIOLIN/tree/master/examples/output), in BioRECIPE format, containing classified interactions (corroborations, contradictions, extensions, and flagged) and [example/output/<interactions_filename>_<class>_score.csv](https://github.com/pitt-miskov-zivanov-lab/VIOLIN/tree/master/examples/output) (evidence score, match score, kind score, epistemic score and total score) 
- `use_VIOLIN.ipynb` - [Jupyter notebook](https://github.com/pitt-miskov-zivanov-lab/VIOLIN/blob/master/examples/use_VIOLIN.ipynb) used to run VIOLIN
- `use_violin_script.py` - [python script](https://github.com/pitt-miskov-zivanov-lab/VIOLIN/blob/master/examples/use_violin_script.py) used to run VIOLIN at the command line

### I/O Annotations

| Input Annotations      | Query                                                                    | Curation Method |
|------------------------|--------------------------------------------------------------------------|-----------------|
| RA1_reading.xlsx       | Melanoma                                                                 | REACH           |
| RA2_reading.xlsx       | MEK, ERK, AKT, GSK3, P70RSK, S6, CDK4, 4EBP1, YB1, SRC, CHK2, MTOR, PI3K | REACH explorer  |
| RA2_0_1_reading.xlsx   | RA2 without LEEs involving biological processes or chemicals             | Manual          |
| RA2_0_1_1_reading.xlsx | RA2_0_1 without LEEs which were redundant or irrelevant                  | Manual          |
| RA3_reading.xlsx       | MAPK/ERK pathway                                                         | REACH explorer  |
| RA4_reading.xlsx       | RPS6K1                                                                   | REACH explorer  |

For each reading set `R#` in `test_input`

| Output Annotations    | Contents                                     |
|-----------------------|----------------------------------------------|
| R#_corroborations.csv | All LEEs judged as collaborations            |
| R#_contradictions.csv | All LEEs judged as contradictions            |
| R#_extensions.csv     | All LEEs judged as extensions                |
| R#_flagged.csv        | All LEEs judged as flagged                   |
| R#_outputDF.csv       | All LEEs, in order of descending Total Score |

## Online Tutorial
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pitt-miskov-zivanov-lab/VIOLIN/HEAD?labpath=%2Fexamples%2Fuse_VIOLIN.ipynb)

Run the demonstrated example; or alternatively upload user-customized input files (see [I/O](#IO)) to the _input/_ directory on File Browser Tab (upper left corner) of Binder.

#### This interactive jupyter notebook walks you though all of the code and functions to:


## Offline Installation

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

## Package Structure

- [`setup.py`](setup.py): python file that help set up python dependencies installation
- [`src/violin/`](src/violin/): directory that includes core python VIOLIN files
  - [`src/violin/formatting.py`](src/violin/formatting.py): functions of preprocessing strings in tabular input;
  - [`src/violin/in_out.py`](src/violin/in_out.py): functions of reading interactions list and model file and writing VIOLIN outputs;
  - [`src/violin/network.py`](src/violin/in_out.py): functions of creating model network and finding paths between nodes;
  - [`src/violin/scoring.py`](src/violin/in_out.py): implementation of decision tree for classification;
  - [`src/violin/visualize_violin.py`](src/violin/in_out.py): functions of visualizing classifying results;
  - 
- [`examples/`](examples/): directory that includes tutorial notebook and example inputs and outputs
- [`environment.yml`](environment.yml): environment file, required by [Binder](https://mybinder.readthedocs.io/en/latest/using/config_files.html#environment-yml-install-a-conda-environment)
- [`docs/`](docs/): containing files supporting the repo's host on [Read the Docs](https://melody-violin.readthedocs.io)
- [`LICENSE.txt`](LICENSE.txt): MIT License
- [`README.md`](README.md): this is me!

## Reproducibility
The data and code for running experiments and creating figures in paper can be found in `data.zip`. The code is tested by `example/test` and `example/test_VIOLIN`.


## Citation

```bibtex
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

This work was funded in part by DARPA Big Mechanism award, AIMCancer (W911NF-17-1-0135) and in part by the NSF EAGER award CCF-2324742.

## Support
Feel free to reach out via email nmzivanov@pitt.edu for additional support if you run into any error.
