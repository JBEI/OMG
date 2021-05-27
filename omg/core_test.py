import cobra
import pytest
import pandas as pd
import datatest as dt
from pandas.util.testing import assert_frame_equal
import pickle
import os
import sys
import copy

from .core import *

# #=============================================================================
# # FIXTURES FOR TESTS
# #=============================================================================

@pytest.fixture(scope="module")
def user_params_data():
    return {
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
def grid_data(user_params_data):
    t0 = user_params_data['timestart']
    tf = user_params_data['timestop']
    points = user_params_data['numtimepoints']
    tspan, delt = np.linspace(t0, tf, points, dtype='float64', retstep=True)
    grid = (tspan, delt)
    return grid

@pytest.fixture(scope="module")
def model_TS_data():
    cwd = os.getcwd()  # Get the current working directory (cwd)
    with open(os.path.join(cwd, 'omg/integration_tests/data/model_TS.pickle'), 'rb') as model_TS_pickle:
        model_TS = pickle.load(model_TS_pickle)
    return model_TS

@pytest.fixture(scope="module")
def solution_TS_data():
    cwd = os.getcwd()  # Get the current working directory (cwd)
    with open(os.path.join(cwd, 'omg/integration_tests/data/solution_TS_pickle'), 'rb') as solution_TS_pickle:
        solution_TS = pickle.load(solution_TS_pickle)
    return solution_TS

# @pytest.fixture(scope='module')
# @dt.working_directory(__file__)
# def mydata():
#     fname = os.path.join(os.path.dirname(__file__), 'data/Emets.csv')
#     return pd.read_csv(fname)
   
   
#=============================================================================
# TESTS 
#=============================================================================

def test_advance_OD_Emets(user_params_data, grid_data):
    Erxn2Emet =  {
        'EX_glc__D_e': 'glc__D_e',
        'EX_k_e': 'k_e',
        'EX_ac_e': 'ac_e',
        'EX_lac__D_e': 'lac__D_e',
        'EX_mg2_e': 'mg2_e',
        'EX_na1_e': 'na1_e',
        'EX_nh4_e': 'nh4_e',
        'EX_cl_e': 'cl_e',
        'EX_pi_e': 'pi_e',
        'EX_so4_e': 'so4_e',
        'EX_etoh_e': 'etoh_e',
        'EX_for_e': 'for_e',
        'EX_isoprenol_e': 'isoprenol_e'
        }

    tspan = grid_data[0] 
    old_cell = pd.Series(index=tspan)
    t0 = user_params_data['timestart']
    old_cell[t0] = user_params_data['initial_OD']

    met_names = ['glc__D_e', 'nh4_e', 'pi_e', 'so4_e', 'mg2_e', 'k_e', 'na1_e', 'cl_e', 'isoprenol_e', 'ac_e', 'for_e', 'lac__D_e', 'etoh_e']
    initial_concentrations = [22.203, 18.695, 69.454, 2.0, 2.0, 21.883, 103.7, 27.25, 0.0, 0.0, 0.0, 0.0, 0.0]
    old_Emets = pd.DataFrame(index=tspan, columns=met_names)
    old_Emets.loc[t0] = initial_concentrations    
    delt = grid_data[1]

    # solution = []
    # with open('omg/integration_tests/data/solution_0.0.pkl', 'rb') as solution_pickle:
    #     solution_old = pickle.load(solution_pickle)
    # print(solution_old)

    solution = {
        user_params_data['BIOMASS_REACTION_ID']:  
    }

    # making copies of old_cell and old_Emets
    actual_cell = copy.deepcopy(old_cell)
    actual_Emets = copy.deepcopy(old_Emets)

    old_cell[t0+delt], old_Emets.loc[t0+delt] = advance_OD_Emets(Erxn2Emet, old_cell[t0], old_Emets[t0], delt, solution_old, user_params_data)

    # # assert here
    actual_cell[t0+delt] = 0.017098
    Emet_values  = [21.8444, 18.4873, 69.4354, 1.99515, 1.99983, 21.8792, 103.7, 27.2499, 0.00717176, 0.0742613, 0.00358588, 0.277855, 0]
    actual_Emets.loc[t0+delt] = pd.Series(Emet_values, index=tspan)

    # assert actual_cell == old_cell 
    # assert actual_Emets == old_Emets 
























# # NOTE: How to have an elegant solution so that we can add the tolerance criteria to the 
# # test without a brute force approach
# def test_dataframe(mydata):
#     print("===================", os.getcwd())
#     print("===================", os.path.dirname(__file__))
#     print("===================", os.path.join(os.path.dirname(__file__), 'data/Emets_new.csv'))
#     fname = os.path.join(os.path.dirname(__file__), 'data/Emets.csv')
#     df = pd.read_csv(fname)
#     dt.validate(df, mydata)
#     # assert 23 == (21+2)




# # this uses the assert_frame_equal function from andas and uses the 
# # absolute tolerance parameter to account for tolerance
# def test_dataframe_using_pandas(mydata):
#     fname = os.path.join(os.path.dirname(__file__), 'data/Emets.csv')
#     df = pd.read_csv(fname)
#     assert_frame_equal(df, mydata, atol=0.0000001)

# # Custom method to include tolerance for comparing dataframes
# def test_dataframe_custom(mydata):
#    fname = os.path.join(os.path.dirname(__file__), 'data/Emets_new.csv')
#    df = pd.read_csv(fname) 
#    diff = abs(df - mydata)
#    print('=================Dataframe==================')
#    print(diff)