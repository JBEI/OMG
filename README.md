### Dependencies
- pip
- Python 3.x


### Installation

```bash
$ pip install pipenv
$ pipenv install
```

### Running locally

Before you run you need to have a model file in the data directory.
There is a sample file "iECIAI39_1322.xml" in the test directory for you to get started with that you can copy ito the data directory.

You can provide your own file in the data directory and modify the code to hange the file_name. Same edit is required for the reaction name that you want to generate the data from.
The default reaction name used from the test model file is "BIOMASS_Ec_iJO1366_core_53p95M"

```bash
$ cp test/iECIAI39.xml data/ 
$ pipenv shell
$ pipenv run python generate_test_omics_data.py
```

### Output

```sh
List of reactions related to BIOMASS production:
BIOMASS_Ec_iJO1366_WT_53p95M: E. coli biomass objective function (iJO1366) - WT - with 53.95 GAM estimate
BIOMASS_Ec_iJO1366_core_53p95M: E. coli biomass objective function (iJO1366) - core - with 53.95 GAM estimate
Condition parameter:  1
Displaying the reaction bounds after constraining them:
(-15, -15) 
```

### Jupyter notebook
- You can serve a local Jupyter server and run it locally
- Or you can copy the Jpyter notebook to Jupyter hub and run them there