import pytest
import pandas as pd
import datatest as dt
from pandas.util.testing import assert_frame_equal
import pickle
import os
import sys

# sys.path.append('../src/')
# sys.path.append('/Users/somtirtharoy/workspace/Projects/OMG/omg/')

from .utils import *


# FIXTURES

@pytest.fixture(scope="module")
def get_inchikey_to_cid_data():
    
    return {
        'LTFMZDNNPPEQNG-KVQBGUIXSA-L': 'CID:135398596',
        'YKBGVTZYEHREMT-KVQBGUIXSA-N': 'CID:135398592',
        'RXKJFZQQPQGTFL-UHFFFAOYSA-N': 'CID:670',
        'PHNGFPPXDJJADG-RRKCRQDMSA-L': 'CID:135398613',
        'NBBJYMSMWIIQGU-UHFFFAOYSA-N': 'CID:527',
        'VGONTNSXDCQUGY-RRKCRQDMSA-N': 'CID:135398593',
        'GPRLSGONYQIRFK-UHFFFAOYSA-N': 'CID:1038',
        'XMIIGOLPHOKFCH-UHFFFAOYSA-M': 'CID:4740700',
        'UFHFLCQGNIYNRP-UHFFFAOYSA-N': 'CID:783'
    }

@pytest.fixture(scope="module")
def get_experiment_description_file():
    
    return [
        "Line Name, Line Description, Part ID, Media, Shaking Speed, Starting OD, Culture Volume, Flask Volume, Growth Temperature, Replicate Count\n",
        "WT, Wild type E. coli, ABFPUB_000310, M9, 1, 0.1, 50, 200, 30, 1\n"
    ]






# TESTS
# =======================================================

def test_read_pubchem_id_file(get_inchikey_to_cid_data):
    fname = os.path.join(os.path.dirname(__file__), 'integration_tests/data/inchikey_to_cid_test.txt')
    inchikey_to_cid = read_pubchem_id_file(fname) 
    assert inchikey_to_cid == get_inchikey_to_cid_data
    
def test_write_experiment_description_file(get_experiment_description_file):
    output_file_path = 'data/'
    write_experiment_description_file(output_file_path, line_name='WT', label='')

    with open(f'{output_file_path}/EDD_experiment_description_file.csv') as fh:
        lines = fh.readlines()
        assert lines == get_experiment_description_file 

# def test_write_in_al_format():
#     pass
        
# def test_write_in_edd_format(time_series_omics_data, omics_type, user_params, line_name, label=''):
#     pass
                    
# def test_write_omics_files(time_series_omics_data, omics_type, user_params, line_name='WT', al_format=False, label=''):
#     pass
    
# def test_write_OD_data(cell, output_file_path, line_name='WT', label=''):
#     pass
    
# def test_write_training_data_with_isopentenol(df, filename, output_file_path):
#     pass

# def test_write_external_metabolite(substrates, output_file_path, output_metabolites, line_name='WT', label=''):
#     pass   

# def test_get_random_number():
#     pass

# def test_add_random_noise():
#     pass

# def test_get_list_of_reactions(file_name):
#     pass

# def test_read_model(file_name):
#     pass