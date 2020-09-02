"""
omg: Omics Mock Generator

Generates a mock dataset of omics data (importable in EDD):
transcriptomics, proteomics, and metabolomics

Requirements: Python 3.7.2, cobra, numpy, pandas.
"""

__author__ = 'LBL-QMM'
__copyright__ = 'Copyright (C) 2019 Berkeley Lab'
__license__ = 'GNU Affero General Public License Version 3'
__status__ = 'Alpha'
__date__ = 'Dec 2019'
__version__ = '0.1.1'

import argparse
import collections as col
import os
import random
import re
import statistics
import sys
import urllib.parse
import urllib.request
import warnings
from shutil import copyfile
from enum import Enum
from typing import NewType, Dict, List, Any, OrderedDict, Counter

import cobra
from cobra.util.array import create_stoichiometric_matrix
import numpy as np
import pandas as pd
from cobra.exceptions import OptimizationError, Infeasible

# Type annotations
Filename = NewType('Filename', str)

# Enumerations
class Omics(Enum):
    """Enumeration with supported omics data types."""
    PROTEOMICS = 0
    TRANSCRIPTOMICS = 1
    METABOLOMICS = 2

    def __str__(self):
        return f'{str(self.name).lower()}'


# Constants
UNIPROT_URL = '''https://www.uniprot.org/uploadlists/'''
CTS_URL = '''https://cts.fiehnlab.ucdavis.edu/rest/convert/'''
# HOST NAME
HOST_NAME: str = 'ropacus' 
# TODO: Move some constants to variables by program arguments
DATA_FILE_PATH: Filename = Filename('data')
# Output file path
OUTPUT_FILE_PATH: Filename = Filename('data/output')
# INCHIKEY_TO_CID_MAP_FILE_PATH: mapping file path to map inchikey to cids
INCHIKEY_TO_CID_MAP_FILE_PATH: Filename = Filename('mapping') 
# MODEL_FILENAME: Filename = Filename('iECIAI39_1322.xml')  # E. coli
MODEL_FILENAME: Filename = Filename('reannotated_base_v3.sbml')  # R. opacus
MODEL_FILEPATH: Filename = Filename('')
# Training file name
TRAINING_FILE_NAME: Filename = Filename('')
TRAINING_FILE_PATH: Filename = Filename('')
# Start time and stop time
TIMESTART: float = 0.0
TIMESTOP: float = 8.0
NUMPOINTS: int = 9
    
# Initial OD value
INITIAL_OD = 0.01
# number of reactions and instances
NUM_REACTIONS: int = None
NUM_INSTANCES: int = None

# NOTE: user input to the program
REACTION_ID_ECOLI: str = 'BIOMASS_Ec_iJO1366_core_53p95M'  # E. coli
REACTION_ID: str = 'biomass_target'  # R. opacus
# REACTION_ID: str = 'SRC_C00185_e'  # R. opacus
GENE_IDS_DBS: List[str] = ['kegg.genes']  # R. opacus
# GENE_IDS_DBS: List[str] = ['uniprot', 'goa', 'ncbigi']  # E. coli
UNITS: Dict[Omics, str] = {
    Omics.PROTEOMICS: 'proteins/cell',
    Omics.TRANSCRIPTOMICS: "FPKM",
    Omics.METABOLOMICS: "mg/L"
}
# Fix the flux value to -15 as we have data for this constraint
LOWER_BOUND: int = -15
UPPER_BOUND: int = -15

# Internals
_EPS = np.finfo(np.double).eps


def ansi(num: int):
    """Return function that escapes text with ANSI color n."""
    return lambda txt: f'\033[{num}m{txt}\033[0m'


# pylint: disable=invalid-name
gray, red, green, yellow, blue, magenta, cyan, white = map(ansi,
                                                           range(90, 98))


# pylint: enable=invalid-name


#=============================================================================

def get_flux_time_series(model, ext_metabolites, grid, user_params):
    '''
    Generate fluxes and OD
    '''
    
    ## First unpack the time steps for the grid provided
    tspan, delt = grid

    ## Create a panda series containing the cell concentation for each time point
    cell = pd.Series(index=tspan)
    cell0 = user_params['initial_OD'] # in gDW/L
    t0 = user_params['timestart']
    cell[t0] = cell0

    ## Create a dataframe that constains external metabolite names and their concentrations
    # First organize external metabolites and their initial concentrations
    met_names = []
    initial_concentrations = []
    for met, init_conc in ext_metabolites.items():
        met_names.append(met)
        initial_concentrations.append(init_conc)
    # Create dataframe containing external metabolites
    Emets = pd.DataFrame(index=tspan, columns=met_names)
    # Add initial concentrations for external metabolites
    Emets.loc[t0] = initial_concentrations    
    # Create Dictionary mapping exchange reactions to the corresponding external metabolite 
    Erxn2Emet = {r.id: r.reactants[0].id for r in model.exchanges if r.reactants[0].id in met_names}
        
    
    ## Create storage for timeseries of models and solutions
    # Model time series
    model_TS = pd.Series(index=tspan)
    # Solution time series
    solution_TS = pd.Series(index=tspan)

    
    ## Main for loop solving the model for each time step and adding the corresponding OD and external metabolites created
    volume = 1.0  # volume set arbitrarily to one because the system is extensive
    for t in tspan:
        # Adding constraints for each time point without permanent changes to the model
        with model:
            for rxn, met in Erxn2Emet.items():
                # For each exchange reaction set lower bound such that the corresponding 
                # external metabolite concentration does not become negative 
                model.reactions.get_by_id(rxn).lower_bound = max(model.reactions.get_by_id(rxn).lower_bound, 
                                                                -Emets.loc[t,met]*volume/cell[t]/delt)
            
            # Calculate fluxes 
            solution_t = model.optimize()
            
            # Store the solution and model for each timepoint for future use (e.g. MOMA)
            solution_TS[t] = solution_t
            model_TS[t] = model.copy()
                        
            # Calculate OD and external metabolite concentrations for next time point t+delta
            cell[t+delt], Emets.loc[t+delt] = advance_OD_Emets(Erxn2Emet, cell[t], Emets.loc[t], delt, solution_t, user_params)
            print(t, solution_t.status, solution_t[user_params['BIOMASS_REACTION_ID']])     # Minimum output for testing
            
    return solution_TS, model_TS, cell, Emets, Erxn2Emet

def advance_OD_Emets(Erxn2Emet, old_cell, old_Emets, delt, solution, user_params):
    # Output is same as input if nothing happens in the if clause
    new_cell  = old_cell
    new_Emets = old_Emets
    
    # Obtain the value of mu (growth rate)
    mu = solution[user_params['BIOMASS_REACTION_ID']]
    
    # Calculate OD and external metabolite concentrations for next step
    if solution.status == 'optimal' and mu > 1e-6:  # Update only if solution is optimal and mu is not zero, otherwise do not update
        # Calculating next time point's OD
        new_cell = old_cell *np.exp(mu*delt) 
        # Calculating external external metabolite concentrations for next time point
        for rxn, met in Erxn2Emet.items():
            new_Emets[met] = max(old_Emets.loc[met]-solution[rxn]/mu*old_cell*(1-np.exp(mu*delt)),0.0)  
    
    return new_cell, new_Emets

def getBEFluxes(model_TS, design, solution_TS, grid):
    ## Unpacking time points grid
    tspan, delt = grid
    
    ## Parameters for flux constraints
    high = 1.1
    low  = 0.50

    ## Unpack information for desired flux changes
    # Get names for reaction targets
    reaction_names =list(design.index[1:])
    # Find number of target reactions and number of designs (or strains changed)
    #n_reactions = design.shape[1] - 1
    #n_instances = design.shape[0] - 1
        
    ## Time series containing the flux solution obtained through MOMA
    solutionsMOMA_TS = pd.Series(index=tspan)

    ## Main loop: for each strain and at each time point, find new flux profile through MOMA    
    #for i in range(0,n_instances):
    for t in tspan:
        model = model_TS[t]
        sol1 = solution_TS[t] # Reference solution calculated for each time point
        with model:
            # Adding the fluxed modifications for chosen reactions
            for reaction in reaction_names:
                flux = sol1.fluxes[reaction]
                lbcoeff =low
                ubcoeff =high
                if flux < 0:
                    lbcoeff = high
                    ubcoeff = low

                reaction_constraint = model.problem.Constraint(model.reactions.get_by_id(reaction).flux_expression, 
                                            lb = sol1.fluxes[reaction]*design[reaction]*lbcoeff,
                                            ub = sol1.fluxes[reaction]*design[reaction]*ubcoeff)
                                            #lb = model.reactions.get_by_id(reaction).lower_bound*design[reaction],
                                            #ub = model.reactions.get_by_id(reaction).upper_bound*design[reaction])
                model.add_cons_vars(reaction_constraint)

            # Reference solution calculated for each time point in above cell for wild type
            #sol1 = solution_TS[t]

            # Moma solution for each time point
            sol2 = cobra.flux_analysis.moma(model, solution=sol1, linear=False) 
                
            # saving the moma solutions across timepoints
            solutionsMOMA_TS[t] = sol2
                
    return solutionsMOMA_TS

def integrate_fluxes(solution_TS, model_TS, ext_metabolites, grid, user_params):
    
    ## First unpack the time steps for the grid provided
    tspan, delt = grid

    ## Create a panda series containing the cell concentation for each time point
    cell = pd.Series(index=tspan)
    cell0 = user_params['initial_OD'] # in gDW/L
    t0 = user_params['timestart']
    cell[t0] = cell0
    
    ## Create a dataframe that constains external metabolite names and their concentrations (DUPLICATED CODE)
    # First organize external metabolites and their initial concentrations
    model = model_TS[0]
    met_names = []
    initial_concentrations = []
    for met, init_conc in ext_metabolites.items():
        met_names.append(met)
        initial_concentrations.append(init_conc)
    # Create dataframe containing external metabolites
    Emets = pd.DataFrame(index=tspan, columns=met_names)
    # Add initial concentrations for external metabolites
    Emets.loc[t0] = initial_concentrations    
    # Create Dictionary mapping exchange reactions to the corresponding external metabolite 
    Erxn2Emet = {r.id: r.reactants[0].id for r in model.exchanges if r.reactants[0].id in met_names}
    
    ## Main loop adding contributions for each time step
    for t in tspan:     
        # Calculate OD and external metabolite concentrations for next time point t+delta
        cell[t+delt], Emets.loc[t+delt] = advance_OD_Emets(Erxn2Emet, cell[t], Emets.loc[t], delt, solution_TS[t], user_params)
    
    return cell, Emets

def get_proteomics_transcriptomics_data(model, solution):
    """

    :param model:
    :param solution:
    :param condition:
    :return:
    """

    # pre-determined linear constant (NOTE: Allow user to set this via parameter)
    # DISCUSS!!
    k = 0.8
    q = 0.06

    proteomics = {}
    transcriptomics = {}

    rxnIDs = solution.fluxes.keys()
    for rxnId in rxnIDs:
        reaction = model.reactions.get_by_id(rxnId)
        for gene in list(reaction.genes):

            # this will ignore all the reactions that does not have the gene.annotation property
            # DISCUSS!!
            if gene.annotation:
                if 'uniprot' not in gene.annotation:
                    if 'goa' in gene.annotation:
                        protein_id = gene.annotation['goa']
                    else:
                        break
                else:
                    protein_id = gene.annotation['uniprot'][0]
                
                # add random noise wjhich is 5 percent of the signal
                noiseSigma = 0.05 * solution.fluxes[rxnId]/k;
                noise = noiseSigma*np.random.randn();
                proteomics[protein_id] = (solution.fluxes[rxnId]/k) + noise

            # create transcriptomics dict
            noiseSigma = 0.05 * proteomics[protein_id]/q;
            noise = noiseSigma*np.random.randn();
            transcriptomics[gene.id] = (proteomics[protein_id]/q) + noise

    return proteomics, transcriptomics

def get_metabolomics_data(model, solution, mapping_file):
    """

    :param model:
    :param condition:
    :return:
    """
    metabolomics = {}
    # get metabolites

    # read the inchikey to pubchem ids mapping file
    inchikey_to_cid = {}
    inchikey_to_cid = read_pubchem_id_file(mapping_file)
    
    # create the stoichoimetry matrix fomr the model as a Dataframe and convert all the values to absolute values
    sm = create_stoichiometric_matrix(model, array_type='DataFrame')
    sm = sm.abs()

    # get all the fluxes across reactions from the solution
    fluxes = solution.fluxes
    
    # calculating the dot product of the stoichiometry matrix and the fluxes to calculate the net change 
    # in concentration of the metabolites across reactions
    net_change_in_concentrations = sm.dot(fluxes)  
    # converting all na values to zeroes and counting the total number of changes that happens for each metabolite
    num_changes_in_metabolites = sm.fillna(0).astype(bool).sum(axis=1)
    
    for met_id, conc in net_change_in_concentrations.items():
        metabolite = model.metabolites.get_by_id(met_id)
        
        # if there is an inchikey ID for the metabolite
        if 'inchi_key' in metabolite.annotation:
            # if it is a list get the first element
            if type(metabolite.annotation['inchi_key']) is list:
                inchi_key = metabolite.annotation['inchi_key'][0]
            else:
                inchi_key = metabolite.annotation['inchi_key']
        
            if inchi_key in inchikey_to_cid.keys():
                # if the CID is not in the metabolomics dict keys AND the mapped value is not None and the reactions flux is not 0   
                if (inchikey_to_cid[inchi_key] not in metabolomics.keys()) and (inchikey_to_cid[inchi_key] is not None):
                    metabolomics[inchikey_to_cid[inchi_key]] = conc/num_changes_in_metabolites.iloc[num_changes_in_metabolites.index.get_loc(met_id)]
                    
                elif (inchikey_to_cid[inchi_key] is not None):
                    metabolomics[inchikey_to_cid[inchi_key]] += conc/num_changes_in_metabolites.iloc[num_changes_in_metabolites.index.get_loc(met_id)]
                
    return metabolomics

def get_multiomics(model, solution, mapping_file):
    """

    :param model: cobra model object
    :param solution: solution for the model optimization using cobra
    :param data_type: defines the type of -omics data to generate (all by default)
    :return:
    """

    proteomics = {}
    transcriptomics = {}
    fluxomics = {}
    metabolomics = {}

    proteomics, transcriptomics = get_proteomics_transcriptomics_data(model, solution)

    metabolomics = get_metabolomics_data(model, solution, mapping_file)

    return (proteomics, transcriptomics, metabolomics)

def read_pubchem_id_file(mapping_file):
    inchikey_to_cid = {}
    with open(mapping_file, 'r') as fh:
            try:
                line = fh.readline()
                while line:
                    # checking to ignore inchikey records with no cid mappings
                    if (len(line.split()) > 1):
                        inchikey_to_cid[line.split()[0]] = 'CID:'+line.split()[1]
                    else:
                        inchikey_to_cid[line.strip()] = None

                    line = fh.readline()
            # NOTE: propagated exception, raise
            except Exception as ex:
                print("Error in reading file!")
                print(ex)
    # fh.close()

    return inchikey_to_cid

def write_experiment_description_file(condition=1, line_name='WT'):
    # create the filename
    experiment_description_file_name = f'{OUTPUT_FILE_PATH}/experiment_description_file.csv'

    #write experiment description file
    try:
        with open(experiment_description_file_name, 'w') as fh:
            fh.write(f'Line Name, Line Description, Part ID, Media, Shaking Speed, Starting OD, Culture Volume, Flask Volume, Growth Temperature, Replicate Count\n')
            fh.write(f"{line_name}, KEIO wild type, ABF_001327, M9, 1, 0.1, 50, 200, 30, 1\n")
    except Exception as ex:
        print("Error in writing file!")
        print(ex)

    fh.close()

def write_omics_files(time_series_omics_data, omics_type, output_file_path, line_name='WT'):
    """

    :param dataframe:
    :param data_type:
    :return:
    """

    # create file number two: omics file
    # TODO: Need to change the units to actual relevant units
    unit_dict = { "fluxomics": 'g/L',\
            "proteomics": 'proteins/cell',\
            "transcriptomics": "FPKM",\
            "metabolomics": "mg/L"
            }

    # create the filenames
    omics_file_name: str = f'{output_file_path}/{omics_type}_fakedata_sample.csv'
    if not os.path.isdir(output_file_path):
        os.mkdir(output_file_path)
    
    # open a file to write omics data for each type and for all timepoints and constraints
    try:
        with open(omics_file_name, 'w') as fh:
            fh.write(f'Line Name,Measurement Type,Time,Value,Units\n')
            for timepoint, omics_dict in time_series_omics_data.items():
                dataframe = pd.DataFrame.from_dict(omics_dict, orient='index', columns=[f'{omics_type}_value'])
                for index, series in dataframe.iteritems():
                    for id, value in series.iteritems():
                        fh.write((f'{line_name},{id},{timepoint},{value},{unit_dict[omics_type]}\n'))

    except Exception as ex:
        print("Error in writing file!")
        print(ex)
    

def write_OD_data(cell, output_file_path, line_name='WT'):
    # create the filename
    OD_data_file: str = f'{output_file_path}/OD_fakedata_sample.csv'
    if not os.path.isdir(output_file_path):
        os.mkdir(output_file_path)

    # write experiment description file
    try:
        with open(OD_data_file, 'w') as fh:
            fh.write(f'Line Name,Measurement Type,Concentration,Units,Time,Value\n')
            for index, value in cell.items():
                fh.write((f'{line_name},Optical Density,0.75,g/L,{index},{value}\n'))

    except Exception as ex:
        print("Error in writing OD file")
        print(ex)
    
def write_training_data_with_isopentenol(df, filename):
    filename = f'{OUTPUT_FILE_PATH}/{filename}'
    df.to_csv(filename, header=True, index=False)

def write_external_metabolite(substrates, output_file_path, filename='external_metabolites.csv', linename='WT'):
    # create the filename
    external_metabolites: str = f'{output_file_path}/{filename}'
    if not os.path.isdir(output_file_path):
        os.mkdir(output_file_path)
        
    # get ammonium and glucose from substrates
    glucose = substrates.loc[:, 'glc__D_e']
    ammonium = substrates.loc[:, 'nh4_e']
    isopentenol = substrates.loc[:, 'isoprenol_e']

    try:
        with open(external_metabolites, 'w') as fh:
            # get ammonium and glucose from substrates
            fh.write(f'Line Name,Measurement Type,Time,Value,Units\n')
            for index, value in glucose.items():
                fh.write((f'{linename},CID:5793,{index},{value},mg/L\n'))
                
            for index, value in ammonium.items():
                fh.write((f'{linename},CID:16741146,{index},{value},mg/L\n'))

            # write out isopentenol concentrations
            for index, value in isopentenol.items():
                fh.write((f'{linename},CID:15983957,{index},{value},mg/L\n'))
    
    except Exception as ex:
        print("Error in writing OD file")
        print(ex)

def get_random_number():
    """

    :return:
    """
    random.seed(12312)
    return random.random()

def add_random_noise():
    """

    :return:
    """
    pass

def get_list_of_reactions(file_name):
    """

    :param file_name: Name of the model file (has to be xml for now)
    :return: None (prints the list of reactions that has mass in them)
    """

    # Load model¶depending on the kind of file (the file has to be xml)
    if file_name.endswith(".xml"):
        model = cobra.io.read_sbml_model(file_name)

    # Print out the reaction name and reaction id for all reactions related to BIOMASS production:
    print("List of reactions related to BIOMASS production:")
    for rxn in model.reactions:
        if rxn.name is not None and 'BIOMASS' in rxn.id:
            print("{}: {}".format(rxn.id, rxn.name))

def get_optimized_solution(model, reaction_id):
    """

    :param model:
    :param reaction_id:
    :return solution:
    """

    # fix the flux value to -15 as we have data for this constraint
    model.reactions.get_by_id(reaction_id).lower_bound = self.LOWER_BOUND
    model.reactions.get_by_id(reaction_id).upper_bound = self.UPPER_BOUND
    # print(model.reactions.get_by_id(reaction_id))

    print("Displaying the reaction bounds after constraining them:")
    print(model.reactions.get_by_id(reaction_id).bounds)

    # optimizing the model for only the selected reaction   
    # model.slim_optimize()

    # optimizing model
    solution = model.optimize()

    return solution

def read_model(file_name):
    """

    :param file_name:
    :return model:
    """

    # Load model¶depending on the kind of file
    if file_name.endswith(".xml"):
        model = cobra.io.read_sbml_model(file_name)
    elif file_name.endswith(".json"):
        model = cobra.io.load_json_model(file_name)

    return model

def model_has_IPP_pathway(model):
    '''
    We check if the model has the following reactions if so then it has the isopentenol pathway
    ['HMGCOAS','HMGCOAR','MEVK1','PMD','IPMPP','IPtrpp','IPtex','EX_isoprenol_e']
    '''
    reaction_list = ['HMGCOAS','HMGCOAR','MEVK1','PMD','IPMPP','IPtrpp','IPtex','EX_isoprenol_e']
    model_reactions = [r.id for r in model.reactions]
    for reac in reaction_list:
        if reac not in model_reactions:
            return False
    return True

def add_isopentenol_pathway(model, sce):
    '''
    Add isopentenol pathway by taking it from the model instance of S. cerevisiae,
    we used the iMM904.json model
    '''
    # Load S. cerevisiae model
    #     sce = cobra.io.load_json_model(f'data/{cerevisiae_modelfile}')
    
    # Add mevalonate pathway reactions from S. cerevisiae model
    for x in ['HMGCOAS','HMGCOAR','MEVK1','DPMVD']:
        r = sce.reactions.get_by_id(x).copy()
        r.gene_reaction_rule = ''
        model.add_reaction(r)
    
    # Update gene names
    model.reactions.get_by_id('HMGCOAS').gene_reaction_rule = 'HMGS'
    model.reactions.get_by_id('HMGCOAR').gene_reaction_rule = 'HMGR'
    model.reactions.get_by_id('MEVK1').gene_reaction_rule = 'MK'
    model.reactions.get_by_id('DPMVD').gene_reaction_rule = 'PMD'
    
    # Add IP to model
    m = model.metabolites.ipdp_c.copy()
    m.id = 'ipmp_c'
    m.name = 'Isopentenyl monophosphate'
    m.formula = 'C5H9O4P'
    m.charge = -2
    model.add_metabolites([m])
    # Update PMD reaction to convert mev-5p to IP
    model.reactions.get_by_id('DPMVD').id = 'PMD'
    model.reactions.get_by_id('PMD').add_metabolites({'5dpmev_c': 1.0, '5pmev_c': -1.0,
                                                      'ipdp_c': -1.0, 'ipmp_c': 1.0})
    # Add isoprenol (isopentenol)
    m = model.metabolites.ipmp_c.copy()
    m.id = 'isoprenol_c'
    m.name = 'Isopentenol'
    m.formula = 'C5H10O'
    m.charge = 0
    model.add_metabolites([m])
    # Add phosphatase reaction by AphA
    r = model.reactions.CHLabcpp.copy()
    r.id = 'IPMPP'
    r.name = 'Isopentenyl monophosphate phosphatase'
    r.gene_reaction_rule = 'AphA'
    model.add_reactions([r])
    r.add_metabolites({'chol_p': 1.0, 'atp_c': 1.0, 'chol_c': -1.0, 'adp_c': -1.0, 'h_c': -1.0, 'ipmp_c': -1.0, 'isoprenol_c': 1.0})
    
    # Add periplasmic and extracellular isoprenol
    m = model.metabolites.isoprenol_c.copy()
    m.id = 'isoprenol_p'
    m.compartment = 'p'
    model.add_metabolites([m])
    m = model.metabolites.isoprenol_c.copy()
    m.id = 'isoprenol_e'
    m.compartment = 'e'
    model.add_metabolites([m])
    # Add periplasmic and extracellular transport reactions
    r = model.reactions.ETOHtrpp.copy()
    r.id = 'IPtrpp'
    r.name = 'Isopentenol reversible transport via diffusion (periplasm)'
    r.gene_reaction_rule = ''
    model.add_reactions([r])
    r.add_metabolites({'etoh_p': 1.0, 'etoh_c': -1.0, 'isoprenol_p': -1.0, 'isoprenol_c': 1.0})
    r = model.reactions.ETOHtex.copy()
    r.id = 'IPtex'
    r.name = 'Isopentenol transport via diffusion (extracellular to periplasm)'
    r.gene_reaction_rule = ''
    model.add_reactions([r])
    r.add_metabolites({'etoh_e': 1.0, 'etoh_p': -1.0, 'isoprenol_e': -1.0, 'isoprenol_p': 1.0})
    # Add a boundary reaction
    r = model.reactions.EX_etoh_e.copy()
    r.id = 'EX_isoprenol_e'
    r.name = 'Isopentenol exchange'
    r.gene_reaction_rule = ''
    model.add_reactions([r])
    r.add_metabolites({'etoh_e': 1.0, 'isoprenol_e': -1.0})
    
    # Write model to files
    outputfilename = user_params['modelfile'].split('.')[0] + '_IPP.json'
    cobra.io.save_json_model(model, f'data/{outputfilename}')
    
    return model

#=============================================================================

class Ropacus():

    def __init__(self):
        self.time_series_omics_data = {}
        self.LOWER_BOUND = -15
        self.UPPER_BOUND = -15

    def generate_time_series_data(self, model):

        # intiializing omics dictionaries to contain data across timepoints
        proteomics_list: List = []
        transcriptomics_list: List = []
        fluxomics_list: List = []
        metabolomics_list: List = []

        # generating time series data for the following flux constraints
        # 6, 9, 12, 15 corresponding to the times 0, 3, 6, 9 hours
        # NOTE: The constraints and the timepoints should be supplied as command line inputs
        
        time_series_omics_data = {}
        experiment_timepoints = [0, 3, 6, 9]
        flux_constraints = [6, 9, 12, 15]
        # NOTE; constraints in flux_constraints, think about it
        for i in range(len(flux_constraints)):
            # Set global reactions bounds (in addition to local)
            self.LOWER_BOUND = flux_constraints[i]
            self.UPPER_BOUND = flux_constraints[i]
            cobra_config = cobra.Configuration()
            cobra_config.bounds = self.LOWER_BOUND, self.UPPER_BOUND

            # Print the list of reaction names related to BIOMASS production
            self.print_reactions(model)

            # get fake proteomics data and write it to XLSX file
            condition = 1
            self.generate_mock_data(model, condition)


    def add_random_noise(self):
        # TODO
        """

        :return:
        """
        pass

    def chemical_translation(self, dict_in: Dict[str, Any],
                             fmt_from: str = 'KEGG',
                             fmt_to: str = 'PubChem CID') -> Dict[str, Any]:
        """
        Proxy to UCDavis Chemical Translation Service (CTS). Maps the keys of
        the input dictionary keeping intact the values.

        Default behaviour: map KEGG Compounds into PubChem CIDs

        For details, see https://cts.fiehnlab.ucdavis.edu/services
        """

        dict_out: Dict[str, float] = {}
        print(gray('Mapping metabolites ids using CTS'), end='', flush=True)
        ids_in: List[str] = list(dict_in.keys())
        pattern = re.compile(
            r"""(?:"searchTerm":")(\w+)(?:","results":\[")(\w+)(?:"])""")
        for id_in in ids_in:
            mapping_str: str = f'{fmt_from}/{fmt_to}/{id_in}'
            mapping_data = urllib.parse.quote(mapping_str)
            mapping_req = urllib.request.Request(CTS_URL + mapping_data)
            with urllib.request.urlopen(mapping_req) as map_file:
                mapping = map_file.read().strip().decode('utf-8')
            match: re.Match = pattern.search(mapping)
            if match:
                assert match.group(1) == id_in
                id_out: str = match.group(2)
                if fmt_to == 'PubChem CID':
                    id_out = 'CID:' + id_out
                dict_out[id_out] = dict_in[id_in]
                print(green('.'), end='', flush=True)
                dprint(f'Metabolite {id_in} mapped to {id_out}')
            else:
                print(red('.'), end='', flush=True)
                dprint(yellow(f'Metabolite {id_in} mapping failed!'))
        print(green('OK!'))
        self.vprint(gray('Number of unmapped genes from'), fmt_from, gray('to'),
               fmt_to, gray(':'), yellow(len(dict_in) - len(dict_out)))
        return dict_out

    def dict_to_edd(self, omics_dict: Dict[str, float],
                    omics: Omics) -> pd.DataFrame:
        """Get dataframe with EDD format from dictionary with omics values"""

        edd: List[OrderedDict[str, Any]] = []
        sample: OrderedDict[str, Any]

        for measurement, value in omics_dict.items():
            sample = col.OrderedDict([
                ('Line Name', 'WT'),
                ('Measurement Type', measurement),
                ('Time', 0),  # TODO: Generalize for time-series
                ('Value', value),
                ('Units', UNITS[omics])
            ])
            edd.append(sample)

        return pd.DataFrame(edd)

    def dprint(self, *a, **k):
        """Print only if debug mode is enabled"""
        if args.debug:
            print(*a, **k)

    def generate_mock_data(self, model, cond):
        """

        :param model: cobra model object
        :param solution: solution for the model optimization using cobra
        :param data_type: defines the type of -omics data to generate (all by default)
        :return:
        """

        while cond:
            print(gray('Condition parameter:'), magenta(cond))
            cond -= 1
            self.optimize_solution(model, REACTION_ID)
            solution: cobra.Solution = cobra.core.solution.get_solution(
                model, raise_error=False)
            self.vprint(gray('Solution objective value:'), solution.objective_value)
            self.vprint(gray('Model summary after optimization:'))
            try:
                self.vprint(model.summary())
            #   self.vprint(model.metabolites.C00185_e.summary())
            except Infeasible:
                self.vprint(yellow(
                    'Model summary unavailable as solution was unfeasible!'))
                    # exit code here

            self.write_experiment_description(cond)
            self.get_omics_data(model, solution, cond)

    def gene_to_protein(self, dict_in: Dict[str, Any],
                        fmt_gene: str = 'KEGG_ID',
                        fmt_prot: str = 'ID') -> Dict[str, Any]:
        """
        From any dict whose keys are gene IDs, maps them to protein IDs and
        keeps the value intact


        Default behaviour: map KEGG IDs into UNIPROT IDs

        For details, see https://www.uniprot.org/help/api_idmapping
        """

        dict_out: Dict[str, float] = {}
        print(gray('Mapping genes into proteins using UNIPROT... '), end='')
        gene_ids: List[str] = list(dict_in.keys())
        mapping_params: Dict[str, str] = {
            'from': fmt_gene,
            'to': fmt_prot,
            'format': 'tab',
            'query': '\t'.join(gene_ids)
        }
        mapping_data = urllib.parse.urlencode(mapping_params)
        mapping_data = mapping_data.encode('utf-8')
        mapping_req = urllib.request.Request(UNIPROT_URL, mapping_data)
        with urllib.request.urlopen(mapping_req) as map_file:
            mapping = map_file.read().strip().decode('utf-8').split('\n')
        for gene2prot in mapping[1:]:
            gene, prot = gene2prot.split('\t', 1)
            dict_out[prot] = dict_in[gene]
            dprint('Gene', gene, 'mapped to protein', prot)
        if dict_out:
            print(green('OK!'))
            self.vprint(gray('Number of unmapped genes from'), fmt_gene, gray('to'),
                   fmt_prot, gray(':'), yellow(len(dict_in) - len(dict_out)))
        else:
            print(yellow('PROBLEM!'))
        return dict_out

    # NOTE: Name it consistently , generate_omics_data
    def get_omics_data(self, model: cobra.Model,
                       solution: cobra.Solution,
                       cond: int):
        """
        Core method that generates all omics data.

        :param model:
        :param solution:
        :param cond:
        :return:

        """
        # Pre-determined linear constants
        PROTE_SCALING: float = 10  # Scaling factor for fluxes to proteomics
        TRANS_SCALING: float = 1.2  # S.F. for proteomics to transcriptomics
        # TODO: Allow user to set those constants via parameters

        # The omics variable name should coincide with those elements of Omics
        proteomics: Dict[str, float] = {}
        transcriptomics: Dict[str, float] = {}
        metabolomics: Dict[str, float] = {}

        # Get values and statistics for proteomics and transcriptomics
        proteo_stats: Dict[str, Counter[str]] = {
            db + status: col.Counter() for db in GENE_IDS_DBS
            for status in ['_missing', '_success', '_zero']}
        metabolite_awflux: Dict[str, List[float]] = {}  # abs weighted fluxes

        rxn_ids: pd.Index = solution.fluxes.index
        # Cobra docs: Accessing reaction fluxes through a Solution object
        #  is the safer, preferred, and only guaranteed to be correct way.
        # NOTE: Put the operations in fucntions , more modular
        for rxn_id in rxn_ids:
            reaction: cobra.Reaction = model.reactions.get_by_id(rxn_id)
            flux: float = solution.fluxes[rxn_id]
            gene: cobra.Gene

            # Subloop 1/2: proteomics and transcriptomics
            for gene in reaction.genes:
                gene_id: str = ''
                # WARNING! Based on gene.annotation property populated
                gene_id_db: str = ''
                for gene_id_db in GENE_IDS_DBS:
                    try:
                        gene_id = gene.annotation[gene_id_db]
                    except KeyError:
                        proteo_stats[gene_id_db + '_missing'][gene_id] += 1
                    else:
                        # Populates proteomics and transcriptomics dicts if
                        #  related flux has a positive value
                        proteo: int = np.ceil(flux * PROTE_SCALING)
                        if proteo > _EPS:
                            # Accumulate in case of multiple genes
                            try:
                                proteomics[gene_id] += proteo
                            except KeyError:
                                proteomics[gene_id] = proteo
                            proteo_stats[gene_id_db + '_success'][gene_id] += 1
                        else:
                            proteo_stats[gene_id_db + '_zero'][gene_id] += 1
                        transc: float = proteo * TRANS_SCALING
                        if transc > _EPS * 1e+3:
                            transcriptomics[gene.id] = transc
                            break
                        else:
                            self.dprint(yellow('WARNING!'), gray('Gene'), gene.id,
                                gray('in reaction'), rxn_id,
                                gray('has no useful annotation. Skipping...'))

            # Subloop 2/2: metabolomics (partial)
            for metabolite, coeff in reaction.metabolites.items():
                awflux: float = abs(coeff * flux)  # absolute weighted flux
                if awflux < _EPS:
                    continue
                metabolite_id: str = metabolite.id.rsplit(
                    sep='_', maxsplit=1)[0]  # Remove suffixes _c, _e, etc
                try:
                    metabolite_awflux[metabolite_id].append(awflux)
                except KeyError:
                    metabolite_awflux[metabolite_id] = [awflux]

        # Metabolomics (final)
        # Alt: to avoid this loop use a moving average in the subloop above
        for metabolite, awfluxes in metabolite_awflux.items():
            metabolomics[metabolite] = statistics.mean(awfluxes)
        self.vprint(gray('Number of active metabolites:'), len(metabolomics))

        dprint(gray('Number of fluxes related to each gene (top 10)'))
        for gene_id_db in GENE_IDS_DBS:
            for status in ['_missing', '_success', '_zero']:
                self.dprint(gene_id_db + status, proteo_stats[
                    gene_id_db + status].most_common(10))

        # Map genes ids into protein ids accepted by EDD
        proteomics = self.gene_to_protein(proteomics)

        # Map metabolites ids into those accepted by EDD
        metabolomics = self.chemical_translation(metabolomics)

        # Write omics files
        for omic in Omics:  # NOTE: omics variable names are elements of Omics
            omics_df: pd.DataFrame = self.dict_to_edd(eval(str(omic)), omic)
            self.write_data_files(omics_df, omic, cond)

    def get_random_number(self):
        """

        :return:
        """
        random.seed(12312)
        return random.random()

    def optimize_solution(self, model: cobra.Model, reaction_id: str) -> None:
        """

        :param model:
        :param reaction_id:
        :return solution:
        """

        reaction: cobra.Reaction = model.reactions.get_by_id(reaction_id)
        self.vprint(gray('Reaction:'), reaction)
        if args.debug:
            print(blue('List of reactants:'))
            for reactant in reaction.reactants:
                print(reactant, reactant.name)
            print(blue('List of products:'))
            for product in reaction.products:
                print(product, product.name)


        # Set local reaction bounds
        model.reactions.get_by_id(reaction_id).lower_bound = LOWER_BOUND
        model.reactions.get_by_id(reaction_id).upper_bound = UPPER_BOUND
        self.vprint(gray('Displaying the reaction bounds after constraining them:'),
               blue(model.reactions.get_by_id(reaction_id).bounds))
        # Optimize the model using FBA
        print(gray('Optimizing the model using FBA... '), end='')
        model.slim_optimize()
        try:
            cobra.util.assert_optimal(model)
        except OptimizationError as error:
            print(yellow('PROBLEM!'), error)
        else:
            print(green('OK!'))

    def read_model(self, file_name):
        """

        :param file_name:
        :return model:
        """

        # Check presence of model file
        if not os.path.isfile(file_name):
            # NOTE: The error handling not consistent and will be dominated by the stack trace
            print(red('ERROR!'),
                  f'File {file_name} missing from the data dir!')
            raise IOError('Missing file')

        # Load model depending on the kind of file
        self.vprint(gray(f'Loading model in {file_name}... '), end='')
        if file_name.endswith('.xml') or file_name.endswith('.sbml'):
            model = cobra.io.read_sbml_model(file_name)
        elif file_name.endswith('.json'):
            model = cobra.io.load_json_model(file_name)
        else:
            # NOTE: stacktrace issue
            print(red('ERROR!'),
                  f'File {file_name} type not supported!')
            raise TypeError('Unsupported file format')
        self.vprint(green('OK!'))

        return model

    def print_reactions(self, model):
        """

        :param model:
        :return: None (prints the list of reactions that have BIOMASS in them)
        """

        # Print out the reaction name and reaction id for all reactions
        #   related to BIOMASS production:
        self.vprint(gray('List of reactions related to BIOMASS production:'))
        for rxn in model.reactions:
            if rxn.name is not None and 'biomass' in rxn.id.lower():
                self.vprint(f"{rxn.id} : {rxn.name}")


    # NOTE: pass everything to asingle print function and add the verbosity arg layer there
    def vprint(self, *a, **k):
        """Print only if verbose mode is enabled"""
        if args.verbose:
            print(*a, **k)

    def write_data_files(self, edd: pd.DataFrame, omics: Omics = None,
                         cond: int = 1) -> None:
        """
        Write the EDD dataframe into a xlsx file

        :param edd:
        :param omics:
        :param cond:
        :return:
        """
        omics_fname: Filename = Filename(
            os.path.join(DATA_FILE_PATH,
                         f'{omics}_mock{cond}.xlsx'))
        print(gray('Saving file'), magenta(omics_fname), gray('... '), end='')
        try:
            # NOTE: Both excel and CSV for both classes and make this method a part of the core class IMPORTANT!!!
            edd.to_excel(omics_fname,
                         sheet_name=f'{omics}',
                         index=False)
        
        # NOTE: Handle this error better. Handle errors so that you can make this into a library and propagate the errors for better handling
        except IOError as ex:
            print(red('ERROR!'))
            self.vprint(ex)
        else:
            print(green('OK!'))

    def write_experiment_description(self, cond=1):
        """

        :param cond:
        :return:
        """

        exp_desc_fname: Filename = Filename(
            os.path.join(
                DATA_FILE_PATH,
                f'EDD_Omics_Experiment_Description_mock{cond}.xlsx'))
        index_label = 'Line Name'
        exp_desc_cols = pd.Index([
            'Line Description',
            'Media',
            'Shaking speed',
            'Starting OD',
            'Culture Volume',
            'Flask Volume',
            'Growth Temperature',
            'Replicate Count',
        ], name=index_label)
        metadata_wt: Dict[str, Dict[str, Any]] = {'WT': {
            'Line Description': 'R. Opacus PD630 wild type (mock)',
            'Media': 'Mock media',
            'Shaking speed': 1.0,
            'Starting OD': 0.1,
            'Culture Volume': 50.0,
            'Flask Volume': 200.0,
            'Growth Temperature': 30.0,
            'Replicate Count': 1,
        }}
        exp_desc_df = pd.DataFrame.from_dict(metadata_wt,
                                             orient='index',
                                             columns=exp_desc_cols)
        print(gray('Saving file'), magenta(exp_desc_fname),
              gray('... '), end='')
        try:
            exp_desc_df.to_excel(exp_desc_fname,
                                 sheet_name='EXP_DESC',
                                 index_label=index_label)
        except IOError as ex:
            print(red('ERROR!'))
            self.vprint(ex)
        else:
            print(green('OK!'))

#======================================
# MAIN FUNCTION
#======================================

def check_debug(args):
    """Check debugging mode"""
    if args.debug:
        print(blue('INFO:'), gray('Debugging mode activated'))
        print(blue('INFO:'), gray('Active parameters:'))
        for key, val in vars(args).items():
            if val is not None and val is not False and val != []:
                print(gray(f'\t{key} ='), f'{val}')
        args.verbose = True  # Unconditional verbose mode activation
    elif not sys.warnoptions:
        warnings.simplefilter("ignore")

def generate_data_for_host(filename):
    global HOST_NAME 
    global DATA_FILE_PATH
    global OUTPUT_FILE_PATH
    
    # if data folder doesn't exist create it
    if not os.path.isdir(DATA_FILE_PATH):
        os.mkdir(DATA_FILE_PATH)
    if not os.path.isdir(OUTPUT_FILE_PATH):
        os.mkdir(OUTPUT_FILE_PATH)
        
    # copy the training file to the data folder
    src_file = f'{TRAINING_FILE_PATH}/{TRAINING_FILE_NAME}'
    dest_file = f'{DATA_FILE_PATH}/{TRAINING_FILE_NAME}'
    dest = copyfile(src_file, dest_file)
    
    MODEL_FILEPATH
    src_file = f'{MODEL_FILEPATH}/{MODEL_FILENAME}'
    dest_file = f'{DATA_FILE_PATH}/{MODEL_FILENAME}'
    dest = copyfile(src_file, dest_file)
    
    """
        Generate omics data for host and model name
    """
    if HOST_NAME == 'ecoli':
        # create instance of the E. Coli class
        ecoli = Ecoli()

        # read model file
        model = ecoli.read_model(filename)

        # generate ecoli synthetic data for model and condition
        condition = 1
        ecoli.generate_time_series_data(model, condition)

    elif HOST_NAME == 'ropacus':
        # create instance of the E. Coli class
        rop = Ropacus()

        # read model file
        model = rop.read_model(filename)

        # generate time series mock data for host
        rop.generate_time_series_data(model)

def main():
    """Main entry point to the script."""
    global REACTION_ID_ECOLI
    global DATA_FILE_PATH
    global HOST_NAME 
    global MODEL_FILENAME
    global MODEL_FILEPATH
    global TIMESTART
    global TIMESTOP
    global NUMPOINTS
    global TRAINING_FILE_NAME
    global TRAINING_FILE_PATH
    global INITIAL_OD
    
    # Argument Parser Configuration
    parser = argparse.ArgumentParser(
        description='Omics Mock Generator',
        epilog='%(prog)s -- {}'.format(__date__),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        help='enable debug mode (implies verbose mode)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='enable verbose mode'
    )
    parser.add_argument(
        '-V', '--version',
        action='version',
        version='%(prog)s release {} ({})'.format(__version__, __date__)
    )
    
    # parser.add_argument(
    #     '-ho', '--host',
    #     default='ropacus',
    #     help='specify host organism'
    # )
    # parser.add_argument(
    #     '-mf', '--modelfile',
    #     default='reannotated_base_v3.sbml',
    #     help='specify model file to use, should be in data folder'
    # )
    parser.add_argument(
        '-ho', '--host',
        default='ecoli',
        help='specify host organism'
    )
    parser.add_argument(
        '-mf', '--modelfile',
        # default='iJO1366_MVA.json',
        default='iJO1366_MVA.json',
        help='specify model file to use, should be in data folder'
    )
    parser.add_argument(
        '-mfp', '--modelfilepath',
        # default='iJO1366_MVA.json',
        default='sample_files',
        help='specify model file path to use'
    )
    parser.add_argument(
        '-tstart', '--timestart',
        default=0.0,
        help='specify the start time for generating the time series data'
    )
    parser.add_argument(
        '-tstop', '--timestop',
        default=9.0,
        help='specify the stop time for generating the time series data'
    )
    parser.add_argument(
        '-np', '--numpoints',
        default=9,
        help='specify the number of points between timestart and timestop for which to generate the time serTRAINING_FILE_PATHies data'
    )
    parser.add_argument(
        '-tf', '--trainingfile',
        default='training_data_8genes.csv',
        help='specify the training file name'
    )
    parser.add_argument(
        '-tfp', '--trainingfilepath',
        default='sample_files',
        help='specify the training file path name'
    )
    parser.add_argument(
        '-nr', '--numreactions',
        default=1,
        help='specify the number of reactions in the training file'
    )
    parser.add_argument(
        '-ni', '--numinstances',
        default=1,
        help='specify the number of instances/strains in the training file'
    )

    # user_params = {
    # 'host': 'ecoli', # ecoli or ropacus
    # 'modelfile': 'iJO1366_MVA.json',
    # 'timestart': 0.0,
    # 'timestop': 8.0,
    # 'numpoints': 9,
    # 'reactants': ['glc__D_e', 'nh4_e', 'pi_e', 'so4_e', 'mg2_e', 'k_e', 'na1_e', 'cl_e'],
    # 'initial_substrates': [22.203, 18.695, 69.454, 2.0, 2.0, 21.883, 103.7, 27.25],
    
    # }

    # Parse arguments
    args = parser.parse_args()

    # Program header
    print('\n=-= {} =-= v{} - {} =-= by {} =-=\n'.format(
        sys.argv[0], __version__, __date__, __author__))

    # Select cases depending on the debug flag
    check_debug(args)

    # check if host and model file has been mentioned
    HOST_NAME = args.host
    MODEL_FILEPATH = args.modelfilepath
    MODEL_FILENAME = args.modelfile
    TIMESTART = args.timestart
    TIMESTOP = args.timestop
    NUMPOINTS = args.numpoints 
    TRAINING_FILE_NAME = args.trainingfile
    TRAINING_FILE_PATH = args.trainingfilepath
    NUM_REACTIONS = args.numreactions
    NUM_INSTANCES = args.numinstances
    INITIAL_OD = args.initialod

    filename: Filename = Filename(f'{MODEL_FILEPATH}/{MODEL_FILENAME}')
    
    # get time series omics data for specified host and model
    generate_data_for_host(filename)


if __name__ == "__main__":
    # TODO: Ask for filename and reaction name and then generate the mock data
    main()
