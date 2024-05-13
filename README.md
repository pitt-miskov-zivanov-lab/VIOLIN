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
- [Citation](#Citation)
- [Funding](#Funding)
- [Support](#Support)

## Functionality
- Function1: ???
- Function2: ???

## I/O

### Input/Output
- `test_input` - ??? input files used in VIOLIN publication
- `test_output` - ??? output files generated from files in test_input
- `use_VIOLIN.ipynb` - ??? ipy notebook used to run VIOLIN
- `use_violin_script.py` - ??? python script used to run VIOLIN at the command line

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

1. ???.
2. ???.
3. ???.

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
   jupyter notebook examples/use_VIOIN.ipynb
   ```

## Package Structure

- [`setup.py`](setup.py): python file that help set up python dependencies installation
- [`src/`](src/): directory that includes core python VIOLIN files
  - [`src/violin/formatting.py`](src/violin/formatting.py): functions of xxx;
  - [`src/violin/in_out.py`](src/violin/in_out.py): functions of xxx;
  - ???
- [`examples/`](examples/): directory that includes tutorial notebook and example inputs and outputs
- [`environment.yml`](environment.yml): environment file, required by [Binder](https://mybinder.readthedocs.io/en/latest/using/config_files.html#environment-yml-install-a-conda-environment)
- [`docs/`](docs/): containing files supporting the repo's host on [Read the Docs](https://melody-violin.readthedocs.io)
- [`LICENSE.txt`](LICENSE.txt): MIT License
- [`README.md`](README.md): this is me!

## Citation

_???_

## Funding

This work was funded in part by DARPA Big Mechanism award, AIMCancer (W911NF-17-1-0135); and in part by the University of Pittsburgh, Swanson School of Engineering.

## Support
Feel free to reach out via email nmzivanov@pitt.edu for additional support if you run into any error.
