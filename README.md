# VIOLIN
[![Documentation Status](https://readthedocs.org/projects/theviolin/badge/?version=latest)](https://theviolin.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pitt-miskov-zivanov-lab/VIOLIN/HEAD?labpath=%2Fexamples%2Fuse_VIOLIN.ipynb)

## Description of files and folders

- `test_input` - input files used in VIOLIN publication
- `test_output` - output files generated from files in test_input
- `use_VIOLIN.ipynb` - ipy notebook used to run VIOLIN
- `use_VIOLIN.py` - python script used to run VIOLIN at the command line


## test_input

| File                   | Query                                                                    | Curation Method |
|------------------------|--------------------------------------------------------------------------|-----------------|
| RA1_reading.xlsx       | Melanoma                                                                 | REACH           |
| RA2_reading.xlsx       | MEK, ERK, AKT, GSK3, P70RSK, S6, CDK4, 4EBP1, YB1, SRC, CHK2, MTOR, PI3K | REACH explorer  |
| RA2_0_1_reading.xlsx   | RA2 without LEEs involving biological processes or chemicals             | Manual          |
| RA2_0_1_1_reading.xlsx | RA2_0_1 without LEEs which were redundant or irrelevant                  | Manual          |
| RA3_reading.xlsx       | MAPK/ERK pathway                                                         | REACH explorer  |
| RA4_reading.xlsx       | RPS6K1                                                                   | REACH explorer  |

## test_output
- For each reading set `R#` in `test_input`

| File                  | Contents                                     |
|-----------------------|----------------------------------------------|
| R#_corroborations.csv | All LEEs judged as collaborations            |
| R#_contradictions.csv | All LEEs judged as contradictions            |
| R#_extensions.csv     | All LEEs judged as extensions                |
| R#_flagged.csv        | All LEEs judged as flagged                   |
| R#_outputDF.csv       | All LEEs, in order of descending Total Score |
