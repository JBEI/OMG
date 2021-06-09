import copy
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
    return {
        "host": "ecoli",  # ecoli or ropacus
        "modelfile": "../data/models/iJO1366_MVA.json",
        "cerevisiae_modelfile": "../data/models/iMM904.json",
        "timestart": 0.0,
        "timestop": 8.0,
        "numtimepoints": 9,
        "designsfile": "ice_mo_strains.csv",
        "designsfilepath": "../data/",
        "mapping_file": "../mapping/inchikey_to_cid.txt",
        "output_file_path": "../data/omg_output",
        "edd_omics_file_path": "../data/omg_output/edd/",
        "numreactions": 8,
        "numinstances": 96,
        "ext_metabolites": {
            "glc__D_e": 22.203,
            "nh4_e": 18.695,
            "pi_e": 69.454,
            "so4_e": 2.0,
            "mg2_e": 2.0,
            "k_e": 21.883,
            "na1_e": 103.7,
            "cl_e": 27.25,
            "isoprenol_e": 0.0,
            "ac_e": 0.0,
            "for_e": 0.0,
            "lac__D_e": 0.0,
            "etoh_e": 0.0,
        },
        "initial_OD": 0.01,
        "BIOMASS_REACTION_ID": "BIOMASS_Ec_iJO1366_core_53p95M",
    }


@pytest.fixture(scope="module")
def grid_data(user_params_data):
    t0 = user_params_data["timestart"]
    tf = user_params_data["timestop"]
    points = user_params_data["numtimepoints"]
    tspan, delt = np.linspace(t0, tf, points, dtype="float64", retstep=True)
    grid = (tspan, delt)
    return grid


@pytest.fixture(scope="module")
def model_TS_data():
    cwd = os.getcwd()  # Get the current working directory (cwd)
    # with open(os.path.join(cwd, 'omg/integration_tests/data/model_TS.pickle'), 'rb') as model_TS_pickle:
    #     model_TS = pickle.load(model_TS_pickle)
    # return model_TS
    model = cobra.io.load_json_model(
        os.path.join(cwd, "omg/integration_tests/data/iJO1366_MVA.json")
    )
    return model


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


def test_get_flux_time_series(model_TS_data, solution_TS_data, user_params_data):
    # model = cobra.io.load_json_model('omg/integration_tests/data/iJO1366_MVA.json')
    ext_metabolites = user_params_data["ext_metabolites"]

    t0 = user_params_data["timestart"]
    tf = user_params_data["timestop"]
    points = user_params_data["numtimepoints"]
    grid = np.linspace(t0, tf, points, dtype="float64", retstep=True)

    solution_TS, model_TS, cell, Emets, Erxn2Emet = get_flux_time_series(
        model_TS_data, ext_metabolites, grid, user_params_data
    )

    # Asserts
    # compare fluxes from the solution and compare that
    for i, sol in solution_TS.items():
        if not isinstance(sol, float):
            for index, flux in sol.fluxes.items():
                if (
                    "BIOMASS_Ec_iJO1366_WT_53p95M" in index
                ):  # change the hardcoded the biomass reaction name
                    assert flux == solution_TS_data[str(i)]


def test_integrate_fluxes(model_TS_data, solution_TS_data, grid_data, user_params_data):
    #     # setting the objects that needs toget passed ot the called function
    #     # export the dataframes from the notebook and write them out to be imported as fixtures
    #     # integrate_fluxes(solution_TS, model_TS, ext_metabolites, grid, user_params)

    #     # read model and solution objects form pickle files

    # cell, Emets = integrate_fluxes(solution_TS_data, model_TS_data, user_params_data['ext_metabolites'], grid_data, user_params_data)
    pass
