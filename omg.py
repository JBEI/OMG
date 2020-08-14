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
from enum import Enum
from typing import NewType, Dict, List, Any, OrderedDict, Counter

import cobra
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
# Training file name
TRAINING_FILE_NAME: Filename = Filename('')
# Start time and stop time
TIMESTART: float = 0.0
TIMESTOP: float = 8.0
NUMPOINTS: int = 9

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
# Host class: contaiing all core methods used for 
# generating omics data for any host organism. All hosts should inherit from 
# this class and add methods specific to them.

class Host():

    def __init__(self):
        pass

#=============================================================================

class Ecoli(Host):

    def __init__(self):
        self.time_series_omics_data = {}
        self.LOWER_BOUND = -15
        self.UPPER_BOUND = 1000

    def generate_time_series_data(self, model, condition):
        global TRAINING_FILE_NAME 
        global REACTION_ID_ECOLI
        global DATA_FILE_PATH 

        # intiializing omics dictionaries to contain data across timepoints
        proteomics_list: List = []
        transcriptomics_list: List = []
        fluxomics_list: List = []
        metabolomics_list: List = []
       
        # The whole idea of using "batch simulation concepts from Joonhoon" is to
        # estimate glucose consumption values (written as flux constraints below)
        # is a more realistic way
        
        # In order to update flux values and estimate concentration of glucose, we 
        # assume concentrations at t = 0h to be "subs0" (refer to batch simulation notebook).
        # We also assume "volume" = 1.0 and OD at t = 0h (cell0)
        # Assume time points, you can keep this time points or change them to what Joonhoon has
        # Step1: Evaluate flux of glucose at current time point by this equation:
        # model.reactions.get_by_id(k).lower_bound = max(model.reactions.get_by_id(k).lower_bound,
        #                                                  -subs.loc[t,v]*volume/cell[t]/delt)
        # Step2: Solve the model using the function "get_optimized_solution" present below and generate
        # data for that time point
        # Step 3: Calculate mu where mu = solution(biomass)
        # Step 4: Calculate OD for next time point: cell[t+delt] = cell[t]*np.exp(mu*delt)
        # Step 5: Calculate glucose for next time point t+deltat: 
        # subs.loc[t+delt,v] = max(subs.loc[t,v]-sol[k]/mu*cell[t]*(1-np.exp(mu*delt)),0.0)
        # Go back to step 1
        
        time_series_omics_data = {}

        # time steps to calculate the biomass production for
        t0 = TIMESTART
        tf = TIMESTOP
        points = NUMPOINTS
        tspan, delt = np.linspace(t0, tf, points, dtype='float64', retstep=True)

        # step interval
        # delt = tspan[1] - tspan[0]

        # panda series containing OD values for the timepionts
        cell = pd.Series(index=tspan)
        cell0 = 0.01 # in gDW/L
        cell[t0] = cell0
        
        #======
        # to debug precision error in when creating irrational step timepoints 
        # and they dont match the indexes of the timepoints in the cell Series
        debug = False
        if debug:
            print(tspan)
            print(cell)
            print(type(tspan))
            for t in tspan:
                print(t, )
            sys.exit()
            for t in tspan:
                print(f'{t} ---> {cell[t]}')
                if t < 8.0:
                    cell[t+delt] = cell[t]*np.random.random_sample()
                
            print(cell)
            sys.exit()
        #======

        # reactants
        comp = ['glc__D_e', 'nh4_e', 'pi_e', 'so4_e', 'mg2_e', 'k_e', 'na1_e', 'cl_e']

        # Dataframe containing substrates
        subs = pd.DataFrame(index=tspan, columns=comp)
        # initial substrate vlaues
        subs0 = [22.203, 18.695, 69.454, 2.0, 2.0, 21.883, 103.7, 27.25] # in mM
        subs.loc[t0] = subs0

        # Panda series containing the isopentanol concentrations
        # and solutions to the wild type after adding isopentenol
        # pathway and introducing isopentenol fluzes and solving for the 
        # optimum solution for maimum biomass
        conc_iso = pd.Series(index=tspan)
        conc_iso[tspan[0]] = 0.0

        sol_time_wild = pd.Series(index=tspan)
        
        # exterior substrates
        subs_ext = {r.id: r.reactants[0].id for r in model.exchanges if r.reactants[0].id in comp}

        # NOTE: put the body of the for loop inside a function
        for t in tspan:
            # Not changing the model but adding constraints for each time point
            with model:
                for k, v in subs_ext.items():
                    # why do we set volume to one? Is it arbitrary?
                    volume = 1.0
                    # print(model.reactions.get_by_id(k).lower_bound)
                    # print(-subs.loc[t,v]*volume/cell[t]/delt)

                    # Set global reactions bounds (in addition to local)
                    # set the lower bound to the maximum of the lower bound in the model and the change of glucose 
                    model.reactions.get_by_id(k).lower_bound = max(model.reactions.get_by_id(k).lower_bound, -subs.loc[t,v]*volume/cell[t]/delt)
                    self.LOWER_BOUND = model.reactions.get_by_id(k).lower_bound
                    # print(self.LOWER_BOUND)
                    # self.UPPER_BOUND = -15
                    cobra_config = cobra.Configuration()
                    cobra_config.bounds = self.LOWER_BOUND, self.UPPER_BOUND

                    # get fake proteomics data and write it to XLSX file
                    condition = 1
                    proteomics, transcriptomics, fluxomics, metabolomics, solution = self.generate_fake_data(model, condition)

                    # Step 3: Calculate mu where mu = solution(biomass)
                    mu = solution[REACTION_ID_ECOLI]

                    # Step 4: Calculate OD for next time point: cell[t+delt] = cell[t]*np.exp(mu*delt)
                    # print("===================")
                    # print("t, t+delt, cell[t]: ", t, t+delt, cell[t])
                    # print("cell[t]*np.exp(mu*delt): ", cell[t]*np.exp(mu*delt) )
                    # print("cell[t], mu, delt: ", cell[t], mu, delt)
                    cell[t+delt] = cell[t]*np.exp(mu*delt)
                    # print("t+delt, cell[t+delt]: ", t+delt, cell[t+delt])
                    # if np.isnan(cell[t+delt]):
                    #     print("I am here")
                    #     print("cell: ", cell)
                    #     return

                    # Step 5: Calculate glucose for next time point t+deltat:
                    subs.loc[t+delt,v] = max(subs.loc[t,v]-solution[k]/mu*cell[t]*(1-np.exp(mu*delt)),0.0) 

                # appending the dictionaries to a master list that keeps track of the timepoints associated with the data generated
                proteomics_list.append((proteomics, t))
                transcriptomics_list.append((transcriptomics, t))
                fluxomics_list.append((fluxomics, t))
                metabolomics_list.append((metabolomics, t))
 
                # optimize model using pFBA after inducing isopentenol and formate formation 
                # and get isopentenol concentrations
                # NOTE: pass the training file as an argument from the cli
                training_data_file = f'{DATA_FILE_PATH}/{TRAINING_FILE_NAME}'
                sol_time_wild = self.generate_isopentenol_concentrations(model, sol_time_wild, training_data_file, t, tspan, delt, cell, subs, subs_ext, conc_iso)
                # print(sol_time_wild)

        # generate training data for reactions with isopentenol production after optimizing model using MOMA
        # NOTE: This is not working as I do not have the 'cplex' solver installed and it is looking for a 
        # qp-solver that is not there to solve the MOMA optimization
        # Have to run this in the jprime server
        
        # self.generate_isopentenol_and_solution_for_biomass_using_moma(model, sol_time_wild, training_data_file, tspan, delt, cell, subs, subs_ext)
        # print(type(subs))
        # print(subs)
        # sys.exit()

        time_series_omics_data = {'proteomics': proteomics_list, 'transcriptomics': transcriptomics_list, 'fluxomics': fluxomics_list, 'metabolomics': metabolomics_list}
        
        # write all the data generated
        self.write_experiment_description_file(condition)
        self.write_omics_files(time_series_omics_data)
        self.write_OD_data(cell)
        # write external metabolites in subs: Ammonia and glucose and isoprenol concentrations
        self.write_external_metabolite(subs, conc_iso)

    # This uses the modified E. Coli model that has the added isopentenol pathway
    # QUESTION: Do we need to add the isopentenol pathway to it, if not provided?DO we need to check it?
    def generate_isopentenol_concentrations(self, model, sol_time_wild, training_data_file, timepoint, tspan, delt, cell, subs, subs_ext, conc_iso):
        iso = 'EX_isoprenol_e'
        df = pd.read_csv(training_data_file)

        # Calculating the number of reactions that should be modified (n_genes) and 
        # number of strains for which isoprenol concentration should be estimated 
        n_reactions = df.shape[1] - 1
        n_instances = df.shape[0] - 1
        # print(n_reactions,n_instances)

        # Inserting the isoprenol concentration as the last column in the dataframe
        df.insert(loc=n_reactions+1, column='Isoprenol Concentration (mM)', value=None)

        iso_cons = model.problem.Constraint(model.reactions.EX_isoprenol_e.flux_expression,
                                lb = 0.20)
        model.add_cons_vars(iso_cons)
        for_cons = model.problem.Constraint(model.reactions.EX_for_e.flux_expression,
                                lb = 0.10)
        model.add_cons_vars(for_cons)
        # display(model.summary())
        sol_t = model.optimize()
        # storing the solution for each timepoint which are going to be reference solutions for moma (see below)
        sol_time_wild[timepoint] = sol_t
        mu = sol_t[REACTION_ID_ECOLI]

        if sol_t.status == 'optimal' and mu > 1e-6:
            # Calculating next time point's OD
            cell[timepoint+delt] = cell[timepoint]*np.exp(mu*delt)
            for k, v in subs_ext.items():
                # Calculating substrate's concentration for next time point
                subs.loc[timepoint+delt,v] = max(subs.loc[timepoint,v]-sol_t[k]/mu*cell[timepoint]*(1-np.exp(mu*delt)),0.0)
            if sol_t[iso] > 0:
                # Calculating isoprenol concentration for next time point
                conc_iso.loc[timepoint+delt] = conc_iso.loc[timepoint]-sol_t[iso]/mu*cell[timepoint]*(1-np.exp(mu*delt))
            else:
                conc_iso.loc[0:t] = 0
                conc_iso.loc[timepoint+delt] = conc_iso.loc[timepoint]-sol_t[iso]/mu*cell[timepoint]*(1-np.exp(mu*delt))
        else:
            cell[timepoint+delt] = cell[timepoint]
            for k, v in subs_ext.items():
                subs.loc[timepoint+delt,v] = subs.loc[timepoint,v]
            conc_iso.loc[timepoint+delt] = conc_iso.loc[timepoint]

        return sol_time_wild

    def generate_isopentenol_and_solution_for_biomass_using_moma(self, model, sol_time_wild, training_data_file, tspan, delt, cell, subs, subs_ext):
        model.solver = 'cplex'
        iso = 'EX_isoprenol_e'
        df = pd.read_csv(training_data_file)

        # The original e.coli iJO1366 model does not have the isoprenol pathway. 
        # Thus, performing simple flux balance analysis on the model will not allocate 
        # any flux for isoprenol production reaction. So, we modify the model so that 
        # it produces a small amount of isoprenol. In addition, we force a small amount 
        # of formate production which forces the model to activate the 'PFL' reaction.
        with model:
            # display(model.summary())
            # Constraint to force a small amount of isoprenol production
            iso_cons = model.problem.Constraint(model.reactions.EX_isoprenol_e.flux_expression,
                                        lb = 0.20)
            # Adding the constraint to the model
            model.add_cons_vars(iso_cons)
            # Constraint to force a small amount of formate production which would activate the "PFL" reaction
            for_cons = model.problem.Constraint(model.reactions.EX_for_e.flux_expression,
                                        lb = 0.10)
            # Adding the constraint to the model
            model.add_cons_vars(for_cons)
            WT_FBA_sol = cobra.flux_analysis.pfba(model)
            print(WT_FBA_sol.status, WT_FBA_sol[REACTION_ID_ECOLI], WT_FBA_sol[iso])

        # Calculating the number of reactions that should be modified (n_genes) and 
        # number of strains for which isoprenol concentration should be estimated 
        n_reactions = df.shape[1] - 1
        n_instances = df.shape[0] - 1

        # Inserting the isoprenol concentration as the last column in the dataframe
        df.insert(loc=n_reactions+1, column='Isoprenol Concentration (mM)', value=None)

        # Panda series containing the isopentanol concentrations
        # and solutions to the wild type after adding isopentenol
        # pathway and introducing isopentenol fluzes and solving for the 
        # optimum solution for maimum biomass
        conc_iso = pd.Series(index=tspan)
        conc_iso[tspan[0]] = 0.0

        volume = 1.0
        # For each strain
        for i in range(0,n_instances):
            # At each time point
            for t in tspan:
                # Adding constraints to the model at each time point for each strain without globally changing the model
                with model:
                    for k, v in subs_ext.items():
                        model.reactions.get_by_id(k).lower_bound = max(model.reactions.get_by_id(k).lower_bound,
                                                                -subs.loc[t,v]*volume/cell[t]/delt)
                    # Adding the fluxed modifications for chosen reactions
                    cons1 = model.problem.Constraint(model.reactions.ACCOAC.flux_expression, 
                                                    lb = WT_FBA_sol['ACCOAC']*df.iloc[i,1],
                                                    ub = WT_FBA_sol['ACCOAC']*df.iloc[i,1])
                    model.add_cons_vars(cons1)
                
                    cons2 = model.problem.Constraint(model.reactions.MDH.flux_expression,
                                                    lb = WT_FBA_sol['MDH']*df.iloc[i,2],
                                                    ub = WT_FBA_sol['MDH']*df.iloc[i,2])
                    model.add_cons_vars(cons2)
                
                    cons3 = model.problem.Constraint(model.reactions.PTAr.flux_expression,
                                                    lb = WT_FBA_sol['PTAr']*df.iloc[i,3],
                                                    ub = WT_FBA_sol['PTAr']*df.iloc[i,3])
                    model.add_cons_vars(cons3)
                
                    cons4 = model.problem.Constraint(model.reactions.CS.flux_expression,
                                                    lb = WT_FBA_sol['CS']*df.iloc[i,4],
                                                    ub = WT_FBA_sol['CS']*df.iloc[i,4])
                    model.add_cons_vars(cons4)
                
                    cons5 = model.problem.Constraint(model.reactions.ACACT1r.flux_expression,
                                                    lb = WT_FBA_sol['ACACT1r']*df.iloc[i,5],
                                                    ub = WT_FBA_sol['ACACT1r']*df.iloc[i,5])
                    model.add_cons_vars(cons5)
                
                    cons6 = model.problem.Constraint(model.reactions.PPC.flux_expression,
                                                    lb = WT_FBA_sol['PPC']*df.iloc[i,6],
                                                    ub = WT_FBA_sol['PPC']*df.iloc[i,6])
                    model.add_cons_vars(cons6)
                
                    cons7 = model.problem.Constraint(model.reactions.PPCK.flux_expression,
                                                    lb = WT_FBA_sol['PPCK']*df.iloc[i,7],
                                                    ub = WT_FBA_sol['PPCK']*df.iloc[i,7])
                
                    model.add_cons_vars(cons7)
                
                    cons8 = model.problem.Constraint(model.reactions.PFL.flux_expression,
                                                    lb = WT_FBA_sol['PFL']*df.iloc[i,8],
                                                    ub = WT_FBA_sol['PFL']*df.iloc[i,8])
                
                    model.add_cons_vars(cons8)
                    
                    # Reference solution calculated for each time point in above cell for wild type
                    sol1 = sol_time_wild[t]
                    # print(sol_time_wild)
                    # print(sol1)

                    # Moma solution for each time point
                    sol2 = cobra.flux_analysis.moma(model, solution=sol1, linear=False)
                    mu = sol2[REACTION_ID_ECOLI]
                    print(i,t, sol2.status, mu)
                    if sol2.status == 'optimal' and mu > 1e-6:
                        cell[t+delt] = cell[t]*np.exp(mu*delt)
                        for k, v in subs_ext.items():
                            subs.loc[t+delt,v] = max(subs.loc[t,v]-sol2[k]/mu*cell[t]*(1-np.exp(mu*delt)),0.0)
                        if sol2[iso] > 0:
                            conc_iso.loc[t+delt] = conc_iso.loc[t]-sol2[iso]/mu*cell[t]*(1-np.exp(mu*delt))
                        else:
                            conc_iso.loc[0:t] = 0
                            conc_iso.loc[t+delt] = conc_iso.loc[t]-sol2[iso]/mu*cell[t]*(1-np.exp(mu*delt))
                    else:
                        cell[t+delt] = cell[t]
                        for k, v in subs_ext.items():
                            subs.loc[t+delt,v] = subs.loc[t,v]
                        conc_iso.loc[t+delt] = conc_iso.loc[t]
            
            
            # Storing the final concentration for all strains
            df.iloc[i,9] = conc_iso.iloc[-1]
            print(conc_iso)
            print(i,sol2[iso],conc_iso.iloc[-1])

            # write out the training dataset with isopentenol production concentrations
            # filename = 'training_data_8genes_withiso.csv'
            # self.write_training_data_with_isopentenol(df, filename)

    def generate_fake_data(self, model, condition):
        """

        :param model: cobra model object
        :param solution: solution for the model optimization using cobra
        :param data_type: defines the type of -omics data to generate (all by default)
        :return:
        """

        self.proteomics = {}
        self.transcriptomics = {}
        self.fluxomics = {}
        self.metabolomics = {}

        # reaction_id of choice passed to the function# hardcoded here for this particular file (Need to convert this to an interactive cli program)
        reaction_id = REACTION_ID_ECOLI

        # while condition:
            # print("Condition parameter: ", condition)
        condition-=1
        solution = self.get_optimized_solution(model, reaction_id)
        # solution: cobra.Solution = cobra.core.solution.get_solution(
        #     model, raise_error=False)

        proteomics, transcriptomics, fluxomics = self.get_proteomics_transcriptomics_fluxomics_data(model, solution, condition)
        
        metabolomics = self.get_metabolomics_data(model, condition)
        
        return (proteomics, transcriptomics, fluxomics, metabolomics, solution)

    # NOTE: 
    def read_pubchem_id_file(self):
        inchikey_to_cid = {}
        filename = f'{INCHIKEY_TO_CID_MAP_FILE_PATH}/inchikey_to_cid.txt'
        with open(filename, 'r') as fh:
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

    def get_metabolomics_data(self, model, condition):
        """

        :param model:
        :param condition:
        :return:
        """
        metabolomics = {}
        # get metabolites
        # NOTE: Need to find a better algorithm. This is O(n^3)

        # read the inchikey to pubchem ids mapping file
        inchikey_to_cid = {}
        inchikey_to_cid = self.read_pubchem_id_file()

        for met in model.metabolites:
            # get associated reactions
            for reaction in list(met.reactions):
                # get dictionary of associated metabolites and their concentrations
                for metabolite, conc in reaction._metabolites.items():
                    if metabolite.id == met.id:
                        # map the BIGG ids to CIDs using the inchikeys in the metabolites and the ampping file
                        # that we have generated from Pubchem
                        # remember that not all Inchikeys dont have a mappping to a CIDs and there are
                        # multiple mappings for some Inchikeys
                        if 'inchi_key' in met.annotation:
                            if type(met.annotation['inchi_key']) is list:
                                inchi_key = met.annotation['inchi_key'][0]
                            else:
                                inchi_key = met.annotation['inchi_key']
                            
                            if inchi_key in inchikey_to_cid.keys():
                                if inchikey_to_cid[inchi_key] not in metabolomics.keys():
                                    if inchikey_to_cid[inchi_key] is not None:
                                        metabolomics[inchikey_to_cid[inchi_key]] = abs(conc)
                                else:
                                    if inchikey_to_cid[inchi_key] is not None:
                                        metabolomics[inchikey_to_cid[inchi_key]] += abs(conc)
            # getting number of associated reactions and averaging the metabolic concentration value
            num_reactions = len(list(met.reactions))

            # check if inchi_key attribite present else ignore metabolite
            if 'inchi_key' in met.annotation.keys() and inchi_key in inchikey_to_cid.keys():
                if inchikey_to_cid[inchi_key] is not None:
                    metabolomics[inchikey_to_cid[inchi_key]]/=num_reactions

        return metabolomics

    def get_proteomics_transcriptomics_fluxomics_data(self, model, solution, condition):
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
        fluxomics = {}

        # print(solution.fluxes['EX_cm_e'])
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
                            # print("HERE")
                        else:
                            break
                    else:
                        protein_id = gene.annotation['uniprot'][0]
                        # print("HERERERERERERERER")

                    # create proteomics dict
                    # Adding noise which is 5% of the signal data. signal + signal*0.05 = signal*1.05
                    # print(rxnId)
                    # print(solution.fluxes)
                    # print(type(solution.fluxes))
                    # print(solution.fluxes[rxnId])
                    # print(protein_id)

                    proteomics[protein_id] = (solution.fluxes[rxnId]/k)*1.05
                    fluxomics[rxnId] = solution.fluxes[rxnId]

                # create transcriptomics dict
                transcriptomics[gene.id] = (proteomics[protein_id]/q)*1.05

        return proteomics, transcriptomics, fluxomics

    def write_experiment_description_file(self, condition=1, line_name='WT'):
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

    def write_omics_files(self, time_series_omics_data, condition=1, line_name='WT'):
        """

        :param dataframe:
        :param data_type:
        :param condition:
        :return:
        """

        # create file number two: omics file
        # TODO: Need to change the units to actual relevant units
        unit_dict = { "fluxomics": 'g/L',\
                "proteomics": 'proteins/cell',\
                "transcriptomics": "FPKM",\
                "metabolomics": "mg/L"
                }

        # for each omics type data
        for omics_type, omics_list in time_series_omics_data.items():
            # create the filenames
            omics_file_name: str = f'{OUTPUT_FILE_PATH}/{omics_type}_fakedata_sample_{condition}.csv'
            
            # open a file to write omics data for each type and for all timepoints and constraints
            try:
                with open(omics_file_name, 'w') as fh:
                    fh.write(f'Line Name,Measurement Type,Time,Value,Units\n')
                    for omics_dict, timepoint in omics_list:
                        dataframe = pd.DataFrame.from_dict(omics_dict, orient='index', columns=[f'{omics_type}_value'])
                        for index, series in dataframe.iteritems():
                            for id, value in series.iteritems():
                                fh.write((f'{line_name},{id},{timepoint},{value},{unit_dict[omics_type]}\n'))

            except Exception as ex:
                print("Error in writing file!")
                print(ex)
        
            fh.close()

    def write_OD_data(self, cell, line_name='WT'):
        # create the filename
        OD_data_file: str = f'{OUTPUT_FILE_PATH}/OD_fakedata_sample.csv'

        # write experiment description file
        try:
            with open(OD_data_file, 'w') as fh:
                fh.write(f'Line Name,Measurement Type,Concentration,Units,Time,Value\n')
                for index, value in cell.items():
                    # print(index, value)
                    fh.write((f'{line_name},Optical Density,0.75,g/L,{index},{value}\n'))

        except Exception as ex:
            print("Error in writing OD file")
            print(ex)
        
    def write_training_data_with_isopentenol(self, df, filename):
        filename = f'{DATA_FILE_PATH}/{filename}'
        df.to_csv(filename, header=True, index=False)

    def write_external_metabolite(self, substrates, isopentenol_conc, filename='external_metabolites.csv', linename='WT'):
        # create the filename
        external_metabolites: str = f'{OUTPUT_FILE_PATH}/{filename}'
        # get ammonium and glucose from substrates
        glucose = substrates.loc[:, 'glc__D_e']
        ammonium = substrates.loc[:, 'nh4_e']

        try:
            with open(external_metabolites, 'w') as fh:
                # get ammonium and glucose from substrates
                fh.write(f'Line Name,Measurement Type,Time,Value,Units\n')
                for index, value in glucose.items():
                    fh.write((f'{linename},CID:5793,{index},{value},mg/L\n'))
                    
                for index, value in ammonium.items():
                    fh.write((f'{linename},CID:16741146,{index},{value},mg/L\n'))

                # write out isopentenol concentrations
                for index, value in isopentenol_conc.items():
                    fh.write((f'{linename},CID:15983957,{index},{value},mg/L\n'))
        
        except Exception as ex:
            print("Error in writing OD file")
            print(ex)

    def get_random_number(self):
        """

        :return:
        """
        random.seed(12312)
        return random.random()

    def add_random_noise(self):
        """

        :return:
        """
        pass

    def get_list_of_reactions(self, file_name):
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

    def get_optimized_solution(self, model, reaction_id):
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

    def read_model(self, file_name):
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
    global TRAINING_FILE_NAME 
    global REACTION_ID_ECOLI
    global DATA_FILE_PATH
    global HOST_NAME 
    global MODEL_FILENAME
    global TIMESTART
    global TIMESTOP
    global NUMPOINTS
    global TRAINING_FILE_NAME
    
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
        help='specify the number of points between timestart and timestop for which to generate the time series data'
    )
    parser.add_argument(
        '-tf', '--trainingfile',
        default='training_data_8genes.csv',
        help='specify the training file name placed in the data directory in the OMG library'
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

    # if data folder doesn't exist create it
    if not os.path.isdir(DATA_FILE_PATH):
        os.mkdir(DATA_FILE_PATH)
    if not os.path.isdir(OUTPUT_FILE_PATH):
        os.mkdir(OUTPUT_FILE_PATH)

    # check if host and model file has been mentioned
    HOST_NAME = args.host
    MODEL_FILENAME = args.modelfile
    TIMESTART = args.timestart
    TIMESTOP = args.timestop
    NUMPOINTS = args.numpoints 
    TRAINING_FILE_NAME = args.trainingfile

    filename: Filename = Filename(os.path.join(DATA_FILE_PATH, MODEL_FILENAME))
    # reaction_id = 'EX_glc__D_e'
    
    # get time series omics data for specified host and model
    generate_data_for_host(filename)


if __name__ == "__main__":
    # TODO: Ask for filename and reaction name and then generate the mock data
    main()
