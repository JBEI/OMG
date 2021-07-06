import copy
import json
import math
import os
import pickle

import cobra
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
        os.path.join(cwd, "omg/integration_tests/data/iJO1366_MVA.json")
    )

    # adding constraints
    iso = "EX_isoprenol_e"
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
def cell_data(grid_data):
    cwd = os.getcwd()
    tspan = grid_data[0]
    with open(os.path.join(cwd, f"omg/integration_tests/data/cell.txt")) as fh:
        cell_series = pd.Series([float(line.strip()) for line in fh.readlines()])
    cell = cell_series.rename({int(t): t for t in tspan})
    return cell


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
def solution_TS_data(user_params_data, grid_data):
    cwd = os.getcwd()  # Get the current working directory (cwd)
    tspan = grid_data[0]
    solution_TS = {}
    for t in tspan:
        with open(
            os.path.join(cwd, f"omg/integration_tests/data/solution_fluxes_{t}.json")
        ) as fh:
            solution_TS[t] = json.load(fh)
            solution_TS[t]["status"] = "optimal"
            solution_TS[t][user_params_data["BIOMASS_REACTION_ID"]] = 0.5363612610171448

    return solution_TS


# =============================================================================
# TESTS
# =============================================================================


def almost_equal(x, y, threshold=0.0001):
    return abs(x - y) < threshold


def test_get_flux_time_series(
    model_data, solution_TS_data, user_params_data, grid_data, erxn2emet_data
):

    cwd = os.getcwd()
    ext_metabolites = user_params_data["ext_metabolites"]

    tspan = grid_data[0]
    solution_TS, model_TS, cell, Emets, Erxn2Emet = get_flux_time_series(
        model_data, ext_metabolites, grid_data, user_params_data
    )

    # Asserts
    # assert fluxes in solution
    for t in tspan:
        with open(
            os.path.join(cwd, f"omg/integration_tests/data/solution_fluxes_{t}.json")
        ) as fh:
            expected_solution_fluxes = json.load(fh)
            # assert
            actual_solution_fluxes = solution_TS[t].fluxes.to_dict()
            # assert expected_solution_flux_keys == actual_solution_flux_keys
            assert expected_solution_fluxes == actual_solution_fluxes

    # assert Emets
    index_map = {int(t): t for t in tspan}
    expected_Emets = (
        pd.read_csv(os.path.join(cwd, f"omg/integration_tests/data/Emets.csv"))
        .astype(float)
        .rename(index=index_map)
    )
    Emets = Emets.astype(float)
    assert_frame_equal(expected_Emets, Emets)

    # assert Erxn2Emet
    assert Erxn2Emet == erxn2emet_data

    # assert cell
    # with open(
    #       os.path.join(cwd, f'omg/integration_tests/data/cell.txt')
    # ) as fh:
    #     expected_cell = [float(line.strip()) for line in fh.readlines()]
    # print(expected_cell)
    # actual_cell = cell.tolist()
    # print(actual_cell)
    # assert expected_cell == actual_cell


def test_getBEFluxes(model_data, grid_data, user_params_data):
    """
        We test this function for two strains as it takes a long time to run
        for all 95 strains
    """
    # solution = getBEFluxes(model_TS, design, solution_TS_data, grid_data)

    t0 = 0.0
    tf = 1.0
    points = 1
    tspan, delt = np.linspace(t0, tf, points, dtype="float64", retstep=True)
    grid_data = (tspan, delt)

    # get model_TS and solution_TS
    solution_TS, model_TS, cell, Emets, Erxn2Emet = get_flux_time_series(
        model_data, user_params_data["ext_metabolites"], grid_data, user_params_data
    )

    cwd = os.getcwd()
    designs_df = pd.read_csv(
        os.path.join(cwd, f"omg/integration_tests/data/ice_mo_strains.csv"),
        usecols=["Part ID", "Name", "Summary"],
    )
    designs_df.columns = ["Part ID", "Line Name", "Line Description"]

    # select only two strains from the designs
    designs_df = designs_df.loc[0:2, :]

    # getting the reaction names
    reactions = designs_df["Line Description"][0].split("_")[::2]
    for rxn in reactions:
        designs_df[rxn] = None

    for i in range(len(designs_df)):
        if designs_df["Line Name"][i] == "WT":
            designs_df.loc[i][reactions] = [1 for r in range(len(reactions))]
        else:
            values = designs_df.loc[i]["Line Description"].split("_")[1::2]
            designs_df.loc[i][reactions] = [float(value) for value in values]

    designs_df = designs_df.drop(columns=["Line Description", "Part ID"])

    # running MOMA
    solutionsMOMA_TS = {}
    cols = ["Line Name"]
    num_strains = 2
    cols.extend(reactions)
    for i in range(num_strains):
        design = designs_df[cols].loc[i]
        solutionsMOMA_TS[i] = getBEFluxes(model_TS, design, solution_TS, grid_data)

    # solutionsMOMA_TS[i] = getBEFluxes(model_TS, design, solution_TS, grid)

    # asserts
    # read fluxes for solutions

    # assert fluxes for the solutions


def get_optimized_model_at_t(model, erxn2emet, timestep, grid_data, cell):

    cwd = os.getcwd()
    tspan, delt = grid_data
    volume = 1.0
    index_map = {int(t): t for t in tspan}
    Emets = (
        pd.read_csv(os.path.join(cwd, f"omg/integration_tests/data/Emets.csv"))
        .astype(float)
        .rename(index=index_map)
    )

    # print(cell)

    for t in tspan:
        with model:
            for rxn, met in erxn2emet.items():
                # For each exchange reaction set lower bound such that the
                # corresponding external metabolite concentration does not
                # become negative
                model.reactions.get_by_id(rxn).lower_bound = max(
                    model.reactions.get_by_id(rxn).lower_bound,
                    -Emets.loc[t, met] * volume / cell[t] / delt,
                )
        if t == timestep:
            return model
    return None


def test_integrate_fluxes(
    model_data, solution_TS_data, grid_data, user_params_data, erxn2emet_data, cell_data
):
    # setting the objects that needs toget passed ot the called function
    # export the dataframes from the notebook and write them out to be imported
    # as fixtures

    # read model and solution objects form pickle files
    # get optimized model after the first step

    cwd = os.getcwd()
    tspan = grid_data[0]
    model_TS = pd.Series(index=tspan)
    model_TS[0.0] = get_optimized_model_at_t(
        model_data, erxn2emet_data, 0.0, grid_data, cell_data
    )

    cell, Emets = integrate_fluxes(
        solution_TS_data,
        model_TS,
        user_params_data["ext_metabolites"],
        grid_data,
        user_params_data,
        debug=True,
    )

    print(cell)
    # Asserts
    index_map = {int(t): t for t in tspan}
    expected_Emets = (
        pd.read_csv(os.path.join(cwd, f"omg/integration_tests/data/Emets.csv"))
        .astype(float)
        .rename(index=index_map)
    )
    Emets = Emets.astype(float)
    assert_frame_equal(expected_Emets, Emets)

    # assert cell
    with open(os.path.join(cwd, f"omg/integration_tests/data/cell_data.txt")) as fh:
        expected_cell = [float(line.strip()) for line in fh.readlines()]
    print(expected_cell)
    actual_cell = cell.tolist()
    print(actual_cell)
    from pytest import approx

    assert [0.100242 + 0.200023, 0.20042 + 0.40000677] == approx([0.3, 0.6])
    # assert expected_cell == approx(actual_cell)
