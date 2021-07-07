import json
import os
import sys

import numpy as np
import pandas as pd
import pytest

from .utils import (
    read_pubchem_id_file,
    write_experiment_description_file,
    write_external_metabolite,
    write_OD_data,
)

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
def grid_data(user_params_data):
    t0 = user_params_data["timestart"]
    tf = user_params_data["timestop"]
    points = user_params_data["numtimepoints"]
    tspan, delt = np.linspace(t0, tf, points, dtype="float64", retstep=True)
    grid = (tspan, delt)
    return grid


@pytest.fixture(scope="module")
def cell_data(grid_data):
    cwd = os.getcwd()
    tspan = grid_data[0]
    with open(os.path.join(cwd, "omg/integration_tests/data/cell_data_test.txt")) as fh:
        cell_series = pd.Series([float(line.strip()) for line in fh.readlines()])
    cell = cell_series.rename({int(t): t for t in tspan})
    return cell


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
        "Line Name, Line Description, Part ID, Media, Shaking Speed, "
        "Starting OD, Culture Volume, Flask Volume, Growth Temperature, Replicate Count\n",
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


def test_write_OD_data(cell_data, user_params_data):

    cwd = os.getcwd()
    output_file_path = os.path.join(cwd, "omg/integration_tests/data/")
    write_OD_data(cell_data, output_file_path, line_name="WT", label="")

    # check files
    # read ground truth data
    with open(os.path.join(cwd, "omg/integration_tests/data/EDD_OD_ground.csv")) as fh:
        lines_expected = fh.readlines()

    written_file_path = os.path.join(cwd, "omg/integration_tests/data/EDD_OD.csv")
    with open(written_file_path) as fh:
        lines_actual = fh.readlines()

    try:
        assert lines_actual == lines_expected
    except Exception as ex:
        print(ex)
        print(
            f"The file EDD_OD.csv at {written_file_path} doesn't match the expected file EDD_OD_ground"
        )
    else:
        # delete written file if tests passes
        os.remove(written_file_path)


def test_write_external_metabolite(grid_data):
    Emets = {}
    cwd = os.getcwd()
    output_file_path = os.path.join(cwd, "omg/integration_tests/data/")

    # read Emets
    tspan = grid_data[0]
    index_map = {int(t): t for t in tspan}
    Emets = (
        pd.read_csv(os.path.join(cwd, "omg/integration_tests/data/Emets.csv"))
        .astype(float)
        .rename(index=index_map)
    )

    write_external_metabolite(Emets, output_file_path, line_name="WT", label="_WT")

    # check files
    # read ground truth data
    with open(
        os.path.join(
            cwd, "omg/integration_tests/data/EDD_external_metabolites_WT_ground.csv"
        )
    ) as fh:
        lines_expected = fh.readlines()

    written_file_path = os.path.join(
        cwd, "omg/integration_tests/data/EDD_external_metabolites_WT.csv"
    )
    with open(written_file_path) as fh:
        lines_actual = fh.readlines()

    try:
        assert lines_actual == lines_expected
    except Exception as ex:
        print(ex)
        print(
            f"The file EDD_external_metabolites_WT.csv at \
            {written_file_path} doesn't match the expected file EDD_external_metabolites_WT_ground"
        )
    else:
        # delete written file if tests passes
        os.remove(written_file_path)
