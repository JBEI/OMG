import os

import cobra

from .pathways import add_isopentenol_pathway, model_has_IPP_pathway

# sys.path.append('../src/')
# sys.path.append('/Users/somtirtharoy/workspace/Projects/OMG/omg/')


def test_model_has_IPP_pathway():
    cwd = os.getcwd()
    fname = os.path.join(cwd, "data/tests/iJO1366_test_with_IPP.json")
    model = cobra.io.load_json_model(fname)
    assert model_has_IPP_pathway(model) == True


def test_add_isopentenol_pathway():
    cwd = os.getcwd()
    fname = os.path.join(cwd, "data/tests/iJO1366_test_wo_IPP.json")
    model = cobra.io.load_json_model(fname)
    sce_fname = os.path.join(cwd, "data/tests/iMM904.json")
    sce = cobra.io.load_json_model(sce_fname)
    mod_model = add_isopentenol_pathway(model, sce, {}, write_file=False)
    assert model_has_IPP_pathway(mod_model) == True
