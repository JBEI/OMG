import copy
import json
import math
import os
import pickle
import sys

import cobra
import datatest as dt
import pandas as pd
import pytest
from pandas.util.testing import assert_frame_equal

from ..core import *

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
def solution_TS_data():
    # cwd = os.getcwd()  # Get the current working directory (cwd)
    # with open(os.path.join(cwd, 'omg/integration_tests/data/solution_TS_pickle'), 'rb') as solution_TS_pickle:
    #     solution_TS = pickle.load(solution_TS_pickle)
    # return solution_TS
    solution = {
        "0.0": 0.0,
        "1.0": 0.0,
        "2.0": 0.0,
        "3.0": 0.0,
        "5.0": 0.0,
        "4.0": 0.0,
        "5.0": 0.0,
        "6.0": 0.0,
        "7.0": 0.0,
        "8.0": 0.0,
    }
    return solution


# =============================================================================
# TESTS
# =============================================================================


def test_get_flux_time_series(model_data, solution_TS_data, user_params_data, grid_data):

    cwd = os.getcwd()
    ext_metabolites = user_params_data["ext_metabolites"]

    tspan = grid_data[0]
    solution_TS, model_TS, cell, Emets, Erxn2Emet = get_flux_time_series(
        model_data, ext_metabolites, grid_data, user_params_data
    )

    # for t in grid_data[0]:
    #     with open(os.path.join(cwd, f'omg/integration_tests/data/solution_{t}.pkl'), 'rb') as solution_pickle:
    #         solution = pickle.load(solution_pickle)
    #         solution_fluxes_dict = solution.fluxes.to_dict()
    #         with open(os.path.join(cwd, f'omg/integration_tests/data/solution_fluxes_{t}.json'), 'w') as fh:
    #             json.dump(solution_fluxes_dict, fh)


    # Asserts
    # assert solution
    
    for t in tspan:
        with open(os.path.join(cwd, f'omg/integration_tests/data/solution_fluxes_{t}.json'), 'r') as fh:
            expected_solution_fluxes = json.load(fh)

            # assert 
            actual_solution_fluxes = solution_TS[t].fluxes.to_dict()
            
            expected_solution_flux_keys = list(expected_solution_fluxes.keys())
            actual_solution_flux_keys = list(actual_solution_fluxes.keys())
            assert expected_solution_fluxes == actual_solution_fluxes 


def test_integrate_fluxes(model_data, solution_TS_data, grid_data, user_params_data):
    #     # setting the objects that needs toget passed ot the called function
    #     # export the dataframes from the notebook and write them out to be imported as fixtures
    #     # integrate_fluxes(solution_TS, model_TS, ext_metabolites, grid, user_params)

    #     # read model and solution objects form pickle files

    # cell, Emets = integrate_fluxes(solution_TS_data, model_TS_data, user_params_data['ext_metabolites'], grid_data, user_params_data)
    pass
