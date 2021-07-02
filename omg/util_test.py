import json
import os
import pickle
import sys

import datatest as dt
import pandas as pd
import pytest
from pandas.util.testing import assert_frame_equal

from .utils import *

# sys.path.append('../src/')
# sys.path.append('/Users/somtirtharoy/workspace/Projects/OMG/omg/')
sys.path.append("../")


# FIXTURES


@pytest.fixture(scope="module")
def user_params_data():
    cwd = os.getcwd()

    with open(os.path.join(cwd, "omg/integration_tests/data/user_params.json")) as fh:
        user_param_dict = json.load(fh)
    return user_param_dict


@pytest.fixture(scope="module")
def get_inchikey_to_cid_data():

    return {
        "LTFMZDNNPPEQNG-KVQBGUIXSA-L": "CID:135398596",
        "YKBGVTZYEHREMT-KVQBGUIXSA-N": "CID:135398592",
        "RXKJFZQQPQGTFL-UHFFFAOYSA-N": "CID:670",
        "PHNGFPPXDJJADG-RRKCRQDMSA-L": "CID:135398613",
        "NBBJYMSMWIIQGU-UHFFFAOYSA-N": "CID:527",
        "VGONTNSXDCQUGY-RRKCRQDMSA-N": "CID:135398593",
        "GPRLSGONYQIRFK-UHFFFAOYSA-N": "CID:1038",
        "XMIIGOLPHOKFCH-UHFFFAOYSA-M": "CID:4740700",
        "UFHFLCQGNIYNRP-UHFFFAOYSA-N": "CID:783",
    }


@pytest.fixture(scope="module")
def get_experiment_description_file():

    return [
        "Line Name, Line Description, Part ID, Media, Shaking Speed, Starting OD, Culture Volume, Flask Volume, Growth Temperature, Replicate Count\n",
        "WT, Wild type E. coli, ABFPUB_000310, M9, 1, 0.1, 50, 200, 30, 1\n",
    ]


# TESTS
# =======================================================


def test_read_pubchem_id_file(get_inchikey_to_cid_data):
    fname = os.path.join(
        os.path.dirname(__file__), "integration_tests/data/inchikey_to_cid_test.txt"
    )
    inchikey_to_cid = read_pubchem_id_file(fname)
    assert inchikey_to_cid == get_inchikey_to_cid_data


def test_write_experiment_description_file(get_experiment_description_file):
    output_file_path = "data/"
    write_experiment_description_file(output_file_path, line_name="WT", label="")

    with open(f"{output_file_path}/EDD_experiment_description_file.csv") as fh:
        lines = fh.readlines()
        assert lines == get_experiment_description_file


# def test_write_omics_files(time_series_omics_data, omics_type, user_params, line_name='WT', al_format=False, label=''):

#     write_omics_files(proteomics_timeseries, 'proteomics', user_params, line_name='WT', label='_WT')
#     pass


def test_write_OD_data(user_params_data):
    proteomics_timeseries = {
        0.0: {
            "P77747": 0.0,
            "P02931": 0.0,
            "P02932": 0.0,
            "P06996": 0.0,
            "P68183": 0.0,
            "P68187": 0.0,
            "P0AEX9": 0.0,
            "P02916": 0.0,
            "P02943": 0.0,
        }
    }

    output_file_path = user_params_data["output_file_path"]
    # print(output_file_path)
    # print(os.path.exists(output_file_path))
    # print(os.path.isdir(output_file_path))
    if not os.path.exists(output_file_path):
        # print('===================================================> Here')
        os.makedirs(output_file_path)
    # write_in_edd_format(proteomics_timeseries, 'proteomics', user_params_data, 'WT', 'WT')


# def test_write_training_data_with_isopentenol(df, filename, output_file_path):
#     pass

# def test_write_external_metabolite(substrates, output_file_path, output_metabolites, line_name='WT', label=''):
#     pass
