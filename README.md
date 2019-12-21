# OMG: Omics Mock Generator

Generates a mock dataset of _omics_ data (importable in EDD using the new import format): transcriptomics, proteomics, and metabolomics.

### Dependencies
- pip
- Python 3.6


### Installation

Please follow the next steps after cloning the repository:
```bash
$ pip install pipenv
$ pipenv install
```
You may want to use the `--user` flag in `pip` and prepend `~/.local/bin` to your path.

### Running locally

Before you run you need to have a model file in the data directory.

 * _R. opacus_ (**default**): There is a sample model file "reannotated_base_v3.sbml" in the test directory for you to get started with, that you can copy to the data directory. The only biomass reaction in this model is the default used for optimization.

 * _E. Coli_: There is a sample file "iECIAI39_1322.xml" in the test directory for this model organism. The default reaction name used from this test model file is "BIOMASS_Ec_iJO1366_core_53p95M"

 * Other: You can provide your own file in the data directory and modify the code to change the filename. Same edit is required for the reaction name that you want to generate the data from.

So, you can copy the desired model to the data directory

```bash
$ mkdir data
$ cp test/<selected.model> data/
```
and
```bash
$ pipenv run python omg
```
or
```bash
$ pipenv shell
(omg)$ ./omg

```

### Output

Example for _E. Coli_:

```sh
List of reactions related to BIOMASS production:
BIOMASS_Ec_iJO1366_WT_53p95M: E. coli biomass objective function (iJO1366) - WT - with 53.95 GAM estimate
BIOMASS_Ec_iJO1366_core_53p95M: E. coli biomass objective function (iJO1366) - core - with 53.95 GAM estimate
Condition parameter:  1
Displaying the reaction bounds after constraining them:
(-15, -15) 
```
The fake data gets written in the data directory. The condition parameter is responsible for controlling how many sets of fake or synthesized data you want to generate. Play with it to see multiple sets of data.

### Jupyter notebook (__legacy__)
- You can serve a local Jupyter server and run it locally
- Or you can copy the Jupyter notebook to Jupyter hub and run them there

### Contact
- For questions contact Somtirtha Roy at [somtirtharoy@lbl.gov](somtirtharoy@lbl.gov) or Jose M. Martí at [jmm@lbl.gov](jmm@lbl.gov).
