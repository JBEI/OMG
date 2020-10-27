import pytest
import pandas as pd
import datatest as dt
from pandas.util.testing import assert_frame_equal
import sys
import os
import omg
import pickle
# import capsys


def setup_function(function):
    print("setting up", function)
    user_params = {
        'host': 'ecoli', # ecoli or ropacus
        'modelfile': '../data/models/iJO1366_MVA.json',
        'cerevisiae_modelfile': '../data/models/iMM904.json',
        'timestart': 0.0,
        'timestop': 8.0,
        'numtimepoints': 9,
        'designsfile': 'ice_mo_strains.csv',
        'designsfilepath': '../data/',
        'mapping_file': '../mapping/inchikey_to_cid.txt',
        'output_file_path': '../data/omg_output',
        'edd_omics_file_path': '../data/omg_output/edd/',
        'numreactions': 8,
        'numinstances': 96,
        'ext_metabolites': {
            'glc__D_e': 22.203,
            'nh4_e': 18.695,
            'pi_e': 69.454,
            'so4_e': 2.0,
            'mg2_e': 2.0,
            'k_e': 21.883,
            'na1_e': 103.7,
            'cl_e': 27.25,
            'isoprenol_e': 0.0,
            'ac_e': 0.0,
            'for_e': 0.0,
            'lac__D_e': 0.0,
            'etoh_e': 0.0
        },
        'initial_OD': 0.01,
        'BIOMASS_REACTION_ID': 'BIOMASS_Ec_iJO1366_core_53p95M'
    } 

@pytest.fixture(scope="module")
def user_params_data():
    user_params = {
        'host': 'ecoli', # ecoli or ropacus
        'modelfile': '../data/models/iJO1366_MVA.json',
        'cerevisiae_modelfile': '../data/models/iMM904.json',
        'timestart': 0.0,
        'timestop': 8.0,
        'numtimepoints': 9,
        'designsfile': 'ice_mo_strains.csv',
        'designsfilepath': '../data/',
        'mapping_file': '../mapping/inchikey_to_cid.txt',
        'output_file_path': '../data/omg_output',
        'edd_omics_file_path': '../data/omg_output/edd/',
        'numreactions': 8,
        'numinstances': 96,
        'ext_metabolites': {
            'glc__D_e': 22.203,
            'nh4_e': 18.695,
            'pi_e': 69.454,
            'so4_e': 2.0,
            'mg2_e': 2.0,
            'k_e': 21.883,
            'na1_e': 103.7,
            'cl_e': 27.25,
            'isoprenol_e': 0.0,
            'ac_e': 0.0,
            'for_e': 0.0,
            'lac__D_e': 0.0,
            'etoh_e': 0.0
        },
        'initial_OD': 0.01,
        'BIOMASS_REACTION_ID': 'BIOMASS_Ec_iJO1366_core_53p95M'
    }

    return user_params

@pytest.fixture(scope="module")
def grid_data(user_params_data):
    t0 = user_params['timestart']
    tf = user_params['timestop']
    points = user_params['numtimepoints']
    tspan, delt = np.linspace(t0, tf, points, dtype='float64', retstep=True)
    grid = (tspan, delt)
    return grid

@pytest.fixture(scope="module")
def model_TS_data():
    with open('test/data/model_TS_pickle', 'rb') as model_TS_pickle:
        model_TS = pickle.load(model_TS_pickle)
    return model_TS

@pytest.fixture(scope="module")
def solution_TS_data():
    with open('test/data/solution_TS_pickle', 'rb') as solution_TS_pickle:
        solution_TS = pickle.load(solution_TS_pickle)
    return solution_TS

@pytest.fixture(scope='module')
@dt.working_directory(__file__)
def mydata():
    fname = os.path.join(os.path.dirname(__file__), 'data/Emets.csv')
    return pd.read_csv(fname)
    # return 0

def test_Emets(mydata):
    # setting the objects that needs toget passed ot the called function
    # export the dataframes from the notebook and write them out to be imported as fixtures
    # integrate_fluxes(solution_TS, model_TS, ext_metabolites, grid, user_params)
    
    # cell, Emets = omg.integrate_fluxes(solution_TS, model_TS, ext_metabolites, grid, user_params)
    pass

@pytest.fixture(scope='module')
def refdata():
    return pd.DataFrame(
        data=[
            ['x', 50],
            ['y', 30],
            ['z', 20],
        ],
        columns=['A', 'C'],
    )

def test_total(refdata):
    requirement = refdata['C'].sum()
    dt.validate(requirement, 100)

# NOTE: How to have an elegant solution so that we can add the tolerance criteria to the 
# test without a brute force approach
def test_dataframe(mydata):
    print("===================", os.getcwd())
    print("===================", os.path.dirname(__file__))
    print("===================", os.path.join(os.path.dirname(__file__), 'data/Emets_new.csv'))
    fname = os.path.join(os.path.dirname(__file__), 'data/Emets.csv')
    df = pd.read_csv(fname)
    dt.validate(df, mydata)
    # assert 23 == (21+2)






# def test_myoutput(capsys):  # or use "capfd" for fd-level
#     print("hello")
#     sys.stderr.write("world\n")
#     captured = capsys.readouterr()
#     assert captured.out == "hello\n"
#     assert captured.err == "world\n"
#     print("next")
#     captured = capsys.readouterr()
#     assert captured.out == "next\n"