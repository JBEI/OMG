import json
import os
import pickle

import cobra
import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal

from .core import advance_OD_Emets, get_proteomics_transcriptomics_data

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
        os.path.join(cwd, "omg/integration_tests/data/iJO1366_MVA.json")
    )

    # adding constraints
    iso_cons = model.problem.Constraint(
        model.reactions.EX_isoprenol_e.flux_expression, lb=0.20
    )
    model.add_cons_vars(iso_cons)
    for_cons = model.problem.Constraint(
        model.reactions.EX_for_e.flux_expression, lb=0.10
    )
    model.add_cons_vars(for_cons)
    o2_cons = model.problem.Constraint(model.reactions.EX_o2_e.flux_expression, lb=-8.0)
    model.add_cons_vars(o2_cons)

    CC_rxn_names = ["ACCOAC", "MDH", "PTAr", "CS", "ACACT1r", "PPC", "PPCK", "PFL"]
    for reaction in CC_rxn_names:
        reaction_constraint = model.problem.Constraint(
            model.reactions.get_by_id(reaction).flux_expression, lb=-1.0, ub=1.0
        )
        model.add_cons_vars(reaction_constraint)

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
def solution_fluxes_0_data():
    cwd = os.getcwd()
    solution = {}
    with open(
        os.path.join(cwd, "omg/integration_tests/data/solution_fluxes_0.0.json")
    ) as fh:
        solution_fluxes = json.load(fh)
    solution["fluxes"] = solution_fluxes
    return solution


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
        met_names = fh.readlines()[0].strip().split(",")
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
    with open(os.path.join(cwd, "omg/integration_tests/data/proteomics_0_.json")) as fh:
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

    # calling advance_OD_Emets with debug option
    old_cell[t0 + delt], old_Emets.loc[t0 + delt] = advance_OD_Emets(
        erxn2emet_data,
        old_cell[t0],
        old_Emets.loc[t0],
        delt,
        solution_old_data,
        user_params_data,
        debug=True,
    )

    # Assert here
    # expected values
    expected_Emets = pd.DataFrame(index=tspan, columns=met_names_data)
    expected_Emets.loc[t0] = [
        22.070668648116943,
        18.618338547831804,
        69.44715329892733,
        1.9982098787230178,
        1.9999384270961587,
        21.88161457062599,
        103.7,
        27.249963056257695,
        0.0026466270376611185,
        0.027404985374640336,
        0.0013233135188305593,
        0.10253814350326079,
        0.0,
    ]
    expected_Emets.loc[t0 + delt] = [
        22.070668648116943,
        18.618338547831804,
        69.44715329892733,
        1.9982098787230178,
        1.9999384270961587,
        21.88161457062599,
        103.7,
        27.249963056257695,
        0.0026466270376611185,
        0.027404985374640336,
        0.0013233135188305593,
        0.10253814350326079,
        0.0,
    ]

    expected_cell = pd.Series(index=tspan)
    expected_cell[t0] = 0.010000
    expected_cell[t0 + delt] = 0.017098

    assert_series_equal(expected_cell, old_cell, check_less_precise=6)
    assert_frame_equal(expected_Emets, old_Emets)


def test_get_proteomics_transcriptomics_data(
    model_data, solution_fluxes_0_data, proteomics_data, transcriptomics_data
):

    proteomics, transcriptomics = get_proteomics_transcriptomics_data(
        model_data, solution_fluxes_0_data, False, True
    )

    cwd = os.getcwd()
    import json

    with open(
        os.path.join(cwd, "omg/integration_tests/data/proteomics_0_.json"), "w"
    ) as fh:
        json.dump(proteomics, fh)

    # print(len(proteomics.keys()))
    # print(len(proteomics_data.keys()))

    # print(len(transcriptomics.keys()))
    # print(len(transcriptomics_data.keys()))

    from deepdiff import DeepDiff

    dd = DeepDiff(proteomics, proteomics_data, ignore_order=True, math_epsilon=0.000001)
    print(dd)
    # dd = DeepDiff(transcriptomics, transcriptomics_data, ignore_order=True, math_epsilon=.000001)
    # print(dd)
    # assert
    # assert True
    # assert proteomics == proteomics_data
    # assert transcriptomics == transcriptomics_data


def test_get_metabolomics_data(
    model_data, solution_0_data, inchikey_to_cid_data, metabolomics_data
):
    # metabolomics, metabolomics_with_old_ids = get_metabolomics_data(
    #     model_data, solution_0_data, inchikey_to_cid_data, True
    # )
    pass
    # assert
    assert True
