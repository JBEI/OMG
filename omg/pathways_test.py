import os
import pickle
import sys

import cobra
import datatest as dt
import pandas as pd
import pytest
from pandas.util.testing import assert_frame_equal

from .pathways import *

# sys.path.append('../src/')
# sys.path.append('/Users/somtirtharoy/workspace/Projects/OMG/omg/')


def test_model_has_IPP_pathway():
    fname = os.path.join(
        os.path.dirname(__file__), "integration_tests/data/iJO1366_test_with_IPP.json"
    )
    model = cobra.io.load_json_model(fname)
    assert model_has_IPP_pathway(model) == True


def test_add_isopentenol_pathway():
    fname = os.path.join(
        os.path.dirname(__file__), "integration_tests/data/iJO1366_test_wo_IPP.json"
    )
    model = cobra.io.load_json_model(fname)
    sce_fname = os.path.join(
        os.path.dirname(__file__), "integration_tests/data/iMM904.json"
    )
    sce = cobra.io.load_json_model(sce_fname)
    mod_model = add_isopentenol_pathway(model, sce, {}, write_file=False)
    assert model_has_IPP_pathway(mod_model) == True
