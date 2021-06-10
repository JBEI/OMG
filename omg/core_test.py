import copy
import json
import os
import pickle
import sys

import cobra
import datatest as dt
import pandas as pd
import pytest
from pandas.util.testing import assert_frame_equal, assert_series_equal

from .core import *

# #=============================================================================
# # FIXTURES FOR TESTS
# #=============================================================================


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
def model_data():
    cwd = os.getcwd()  # Get the current working directory (cwd)
    model = cobra.io.load_json_model(
        # os.path.join(cwd, "omg/integration_tests/data/iJO1366_MVA.json")
        os.path.join(cwd, "omg/integration_tests/data/model_with_constraints.json")
    )
    return model


@pytest.fixture(scope="module")
def solution_data():
    cwd = os.getcwd()
    with open(os.path.join(cwd, "omg/integration_tests/data/solution_old.json")) as fh:
        solution = json.load(fh)
    return solution


@pytest.fixture(scope="module")
def solution_old_data():
    cwd = os.getcwd()
    with open(os.path.join(cwd, "omg/integration_tests/data/solution_old.json")) as fh:
        solution_old = json.load(fh)
    return solution_old


@pytest.fixture(scope="module")
def init_met_conc_data():
    cwd = os.getcwd()
    with open(os.path.join(cwd, "omg/integration_tests/data/init_met_conc.csv")) as fh:
        conc = fh.readlines()[0].split(",")
        conc_list = [float(i) for i in conc]
    return conc_list


@pytest.fixture(scope="module")
def erxn2emet_data():
    cwd = os.getcwd()
    with open(os.path.join(cwd, "omg/integration_tests/data/erxn2emet.json")) as fh:
        Erxn2Emet_dict = json.load(fh)
    return Erxn2Emet_dict


@pytest.fixture(scope="module")
def met_names_data():
    cwd = os.getcwd()
    with open(os.path.join(cwd, "omg/integration_tests/data/met_names.csv")) as fh:
        met_names = fh.readlines()[0].split(",")
        met_names_list = [i for i in met_names]
    return met_names_list


@pytest.fixture(scope="module")
def emet_values_data():
    cwd = os.getcwd()
    with open(os.path.join(cwd, "omg/integration_tests/data/emet_values.csv")) as fh:
        emet_values = fh.readlines()[0].split(",")
        emet_values_list = [float(i) for i in emet_values]
    return emet_values_list


@pytest.fixture(scope="module")
def solution_pickle_data():
    cwd = os.getcwd()
    with open(
        os.path.join(cwd, "omg/integration_tests/data/solution_0.0.pkl"), "rb"
    ) as fh:
        solution = pickle.load(fh)
    return solution


@pytest.fixture(scope="module")
def proteomics_data():
    cwd = os.getcwd()
    with open(os.path.join(cwd, "omg/integration_tests/data/proteomics.json")) as fh:
        proteomics = json.load(fh)
    return proteomics


@pytest.fixture(scope="module")
def transcriptomics_data():
    cwd = os.getcwd()
    with open(
        os.path.join(cwd, "omg/integration_tests/data/transcriptomics.json")
    ) as fh:
        transcriptomics = json.load(fh)
    return transcriptomics


@pytest.fixture(scope="module")
def metabolomics_data():
    cwd = os.getcwd()
    with open(os.path.join(cwd, "omg/integration_tests/data/metabolomics.json")) as fh:
        metabolomics = json.load(fh)
    return metabolomics


@pytest.fixture(scope="module")
def inchikey_to_cid_data():
    cwd = os.getcwd()
    with open(
        os.path.join(cwd, "omg/integration_tests/data/inchikey_to_cid.json")
    ) as fh:
        inchikey_to_cid = json.load(fh)
    return inchikey_to_cid


# =============================================================================
# TESTS
# =============================================================================


def test_advance_OD_Emets(
    user_params_data,
    solution_old_data,
    grid_data,
    init_met_conc_data,
    erxn2emet_data,
    met_names_data,
    emet_values_data,
):

    tspan = grid_data[0]
    old_cell = pd.Series(index=tspan)
    t0 = user_params_data["timestart"]
    old_cell[t0] = user_params_data["initial_OD"]

    old_Emets = pd.DataFrame(index=tspan, columns=met_names_data)
    old_Emets.loc[t0] = init_met_conc_data
    delt = grid_data[1]

    # making copies of old_cell and old_Emets
    actual_cell = copy.deepcopy(old_cell)
    actual_Emets = copy.deepcopy(old_Emets)

    # calling advance_OD_Emets with debug option
    old_cell[t0 + delt], old_Emets.loc[t0 + delt] = advance_OD_Emets(
        erxn2emet_data,
        old_cell[t0],
        old_Emets.loc[t0],
        delt,
        solution_old_data,
        user_params_data,
        True,
    )

    # # assert here
    actual_Emets.loc[t0 + delt] = pd.Series(emet_values_data, index=met_names_data)
    actual_cell[t0 + delt] = old_cell[t0 + delt]

    assert_series_equal(actual_cell, old_cell)
    # assert_series_equal(actual_Emets, old_Emets)


# docker-compose exec <container_name> python -m pytest


def test_getBEFLuxes(user_params_data):
    t0 = 0.0
    tf = 1.0
    points = 1
    tspan, delt = np.linspace(t0, tf, points, dtype="float64", retstep=True)
    grid = (tspan, delt)

    num_strains = 2

    # designs_df = pd.read_csv(f'{user_params_data["designsfilepath"]}/{user_params_data["designsfile"]}',
    #                     usecols=['Part ID', 'Name', 'Summary'])
    # designs_df.columns = ['Part ID','Line Name','Line Description']

    # reactions = designs_df['Line Description'][0].split('_')[::2]
    # for rxn in reactions:
    #     designs_df[rxn] = None

    # for i in range(len(designs_df)):
    #     if designs_df['Line Name'][i]=='WT':
    #         designs_df.loc[i][reactions] = [1 for r in range(len(reactions))]
    #     else:
    #         values = designs_df.loc[i]['Line Description'].split('_')[1::2]
    #         designs_df.loc[i][reactions] = [float(value) for value in values]

    # designs_df = designs_df.drop(columns=['Line Description','Part ID'])

    solutionsMOMA_TS = {}
    # for i in range(num_strains):
    #     design = designs_df[cols].loc[i]
    #     if design['Line Name']=='WT':
    #         solutionsMOMA_TS[i] = getBEFluxes(model_TS, design, solution_TS, grid)
    #     else:
    #         solutionsMOMA_TS[i] = getBEFluxes(model_TS, design, solutionHI_TS, grid)

    # solutionsMOMA_TS[i] = getBEFluxes(model_TS, design, solution_TS, grid)

    # asserts


def test_get_proteomics_transcriptomics_data(
    model_data, solution_pickle_data, proteomics_data, transcriptomics_data
):

    proteomics, transcriptomics = get_proteomics_transcriptomics_data(
        model_data, solution_pickle_data, True, False
    )

    #  print(proteomics['O32583'])
    # print(proteomics['O32583'])
    # print(proteomics['P69681'])
    # print(proteomics_data['P69681'])

    print(len(proteomics.keys()))
    print(len(proteomics_data.keys()))

    print(len(transcriptomics.keys()))
    print(len(transcriptomics_data.keys()))

    # assert
    # assert proteomics == proteomics_data
    # assert transcriptomics == transcriptomics_data


def test_get_metabolomics_data(
    model_data, solution_pickle_data, inchikey_to_cid_data, metabolomics_data
):
    metabolomics, metabolomics_with_old_ids = get_metabolomics_data(
        model_data, solution_pickle_data, inchikey_to_cid_data, True
    )

    # assert
    assert metabolomics == metabolomics_data
