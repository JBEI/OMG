"""
omg: Omics Mock Generator

Generates a mock dataset of omics data (importable in EDD):
transcriptomics, proteomics, and metabolomics

Requirements: Python 3.7.2, cobra, numpy, pandas.
"""

__author__ = "LBL-QMM"
__copyright__ = "Copyright (C) 2019 Berkeley Lab"
__license__ = ""
__status__ = "Alpha"
__date__ = "Dec 2019"
__version__ = "0.1.1"


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
from shutil import copyfile
from typing import Any, Counter, Dict, List, NewType, OrderedDict

import cobra
import numpy as np
import pandas as pd
from cobra.exceptions import Infeasible, OptimizationError
from cobra.util.array import create_stoichiometric_matrix

from .utils import *


def get_flux_time_series(model, ext_metabolites, grid, user_params):
    """
    Generate fluxes for all timepoints and corresponding OD data

    :param model: a host model on which to run FBA on and calculate the fluxes
    :param ext_metabolites: external metabolites we are interested in
    :grid: tuple having the time span over which we calculate the time series
            data and the step in betweent he timepoints
    :user_params: user params that have all the user parameters supplied by the
        user to customize the solution
    :return solution:
    :model_TS:
    :cell:
    :Emets:
    :Erxn2Emet:
    """
    ## First unpack the time steps for the grid provided
    tspan, delt = grid

    # Create a panda series containing the cell concentation for each time point
    cell = pd.Series(index=tspan)
    cell0 = user_params["initial_OD"]  # in gDW/L
    t0 = user_params["timestart"]
    cell[t0] = cell0

    # Create a dataframe that constains external metabolite names
    # and their concentrations
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
    # Create Dictionary mapping exchange reactions to
    # the corresponding external metabolite
    Erxn2Emet = {
        r.id: r.reactants[0].id
        for r in model.exchanges
        if r.reactants[0].id in met_names
    }

    ## Create storage for timeseries of models and solutions
    # Model time series
    model_TS = pd.Series(index=tspan)
    # Solution time series
    solution_TS = pd.Series(index=tspan)

    ## Main for loop solving the model for each time step and adding
    # the corresponding OD and external metabolites created
    # volume set arbitrarily to one because the system is extensive
    volume = 1.0
    for t in tspan:
        # for t in [0.0, 1.0, 2.0]:

        # Adding constraints for each time point without
        # permanent changes to the model
        with model:
            for rxn, met in Erxn2Emet.items():
                # For each exchange reaction set lower bound such
                # that the corresponding
                # external metabolite concentration does not become negative
                model.reactions.get_by_id(rxn).lower_bound = max(
                    model.reactions.get_by_id(rxn).lower_bound,
                    -Emets.loc[t, met] * volume / cell[t] / delt,
                )

            # Calculate fluxes
            solution_t = model.optimize()

            # Store the solution and model for each timepoint
            # for future use (e.g. MOMA)
            solution_TS[t] = solution_t
            model_TS[t] = model.copy()

            # Calculate OD and external metabolite concentrations
            # for next time point t+delta
            cell[t + delt], Emets.loc[t + delt] = advance_OD_Emets(
                Erxn2Emet, cell[t], Emets.loc[t], delt, solution_t, user_params
            )

            print(
                t, solution_t.status, solution_t[user_params["BIOMASS_REACTION_ID"]]
            )  # Minimum output for testing

    return solution_TS, model_TS, cell, Emets, Erxn2Emet


def advance_OD_Emets(
    Erxn2Emet,
    old_cell,
    old_Emets,
    delt,
    solution,
    user_params,
    timestep=None,
    debug=False,
):
    """
    Get the concentration of the external metabolites over time

    :param Erxn2Emet:
    :param old_cell:
    :param old_Emets:
    :param delt:
    :param solution:
    :param user_params:

    :return
    :new_cell:
    :new_Emets:
    """

    # Output is same as input if nothing happens
    # in the if solution_status == "optimal" and mu > 1e-6 clause
    new_cell = old_cell
    new_Emets = old_Emets

    # Obtain the value of mu (growth rate)
    mu = solution[user_params["BIOMASS_REACTION_ID"]]

    # Calculate OD and external metabolite concentrations for next step
    # Update only if solution is optimal and mu is not zero,
    # otherwise do not update
    if debug:
        solution_status = solution["status"]
    else:
        solution_status = solution.status

    if solution_status == "optimal" and mu > 1e-6:
        # Calculating next time point's OD
        new_cell = old_cell * np.exp(mu * delt)
        # Calculating external external metabolite concentrations
        # for next time point
        for rxn, met in Erxn2Emet.items():
            new_Emets[met] = max(
                old_Emets.loc[met]
                - solution[rxn] / mu * old_cell * (1 - np.exp(mu * delt)),
                0.0,
            )

    return new_cell, new_Emets


def getBEFluxes(model_TS, design, solution_TS, grid):
    """
    Get the fluxes for the bio-engineered strains

    :param model_TS:
    :param design:
    :param solution_TS:
    :param grid:

    :return
    :solutionsMOMA_TS
    """

    ## Unpacking time points grid
    tspan, delt = grid

    ## Parameters for flux constraints
    high = 1.1
    low = 0.50

    # Unpack information for desired flux changes
    # Get names for reaction targets
    reaction_names = list(design.index[1:])
    # Find number of target reactions and number of designs (or strains changed)
    # n_reactions = design.shape[1] - 1
    # n_instances = design.shape[0] - 1

    # Time series containing the flux solution obtained through MOMA
    solutionsMOMA_TS = pd.Series(index=tspan)

    # Main loop: for each strain and at each time point,
    # find new flux profile through MOMA
    # for i in range(0,n_instances):
    for t in tspan:
        model = model_TS[t]
        # Reference solution calculated for each time point
        sol1 = solution_TS[t]
        with model:
            # Adding the fluxed modifications for chosen reactions
            for reaction in reaction_names:
                flux = sol1.fluxes[reaction]
                lbcoeff = low
                ubcoeff = high
                if flux < 0:
                    lbcoeff = high
                    ubcoeff = low

                reaction_constraint = model.problem.Constraint(
                    model.reactions.get_by_id(reaction).flux_expression,
                    lb=sol1.fluxes[reaction] * design[reaction] * lbcoeff,
                    ub=sol1.fluxes[reaction] * design[reaction] * ubcoeff,
                )

                model.add_cons_vars(reaction_constraint)

            # Reference solution calculated for each time point
            # in above cell for wild type
            # sol1 = solution_TS[t]

            # Moma solution for each time point
            sol2 = cobra.flux_analysis.moma(model, solution=sol1, linear=False)

            # saving the moma solutions across timepoints
            solutionsMOMA_TS[t] = sol2

    return solutionsMOMA_TS


def integrate_fluxes(
    solution_TS, model_TS, ext_metabolites, grid, user_params, debug=False
):
    """
    Get the concentration of the external metabolites over time

    :param solution_TS:
    :param model_TS:
    :param ext_metabolites:
    :param grid:
    :param user_params:

    :return
    :cell:
    :Emets:
    """

    ## First unpack the time steps for the grid provided
    tspan, delt = grid

    ## Create a panda series containing the cell concentation
    # for each time point
    cell = pd.Series(index=tspan)
    cell0 = user_params["initial_OD"]  # in gDW/L
    t0 = user_params["timestart"]
    cell[t0] = cell0

    ## Create a dataframe that constains external metabolite names
    # and their concentrations (DUPLICATED CODE)
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
    # Create Dictionary mapping exchange reactions to the
    # corresponding external metabolite
    Erxn2Emet = {
        r.id: r.reactants[0].id
        for r in model.exchanges
        if r.reactants[0].id in met_names
    }

    ## Main loop adding contributions for each time step
    for t in tspan:
        # Calculate OD and external metabolite concentrations
        # for next time point t+delta

        cell[t + delt], Emets.loc[t + delt] = advance_OD_Emets(
            Erxn2Emet,
            cell[t],
            Emets.loc[t],
            delt,
            solution_TS[t],
            user_params,
            t,
            debug,
        )

    return cell, Emets


def get_optimized_solution(model, reaction_id, lower_bound, upper_bound):
    """
    Get the FBA solution of the model by maximizing biomass production

    :param model:
    :param reaction_id:
    :lower_bound:
    :upper_bound:

    :return solution:
    """

    # fix the flux value to -15 as we have data for this constraint
    model.reactions.get_by_id(reaction_id).lower_bound = lower_bound
    model.reactions.get_by_id(reaction_id).upper_bound = upper_bound
    # print(model.reactions.get_by_id(reaction_id))

    print("Displaying the reaction bounds after constraining them:")
    print(model.reactions.get_by_id(reaction_id).bounds)

    # optimizing the model for only the selected reaction
    # model.slim_optimize()

    # optimizing model
    solution = model.optimize()

    return solution


def get_proteomics_transcriptomics_data(model, solution, add_noise=True, debug=True):
    """
    Get transcriptomics and proteomics data

    :param model:
    :param solution:
    :param noise_zero:

    :return:
    :proteomics:
    :transcriptomics:
    """

    # pre-determined linear constant
    # (NOTE: Allow user to set this via parameter)
    # DISCUSS!!
    k = 0.8
    q = 0.06

    proteomics = {}
    transcriptomics = {}

    if debug:
        rxnIDs = solution["fluxes"].keys()
    else:
        rxnIDs = solution.fluxes.keys()

    for rxnId in rxnIDs:
        reaction = model.reactions.get_by_id(rxnId)
        for gene in list(reaction.genes):

            # this will ignore all the reactions that does not have
            # the gene.annotation property
            # DISCUSS!!
            if gene.annotation:
                if "uniprot" not in gene.annotation:
                    if "goa" in gene.annotation:
                        protein_id = gene.annotation["goa"]
                    else:
                        break
                else:
                    protein_id = gene.annotation["uniprot"][0]

                # add random noise which is 5 percent of the signal
                if debug:
                    reaction_flux = solution["fluxes"][rxnId]
                else:
                    reaction_flux = solution.fluxes[rxnId]

                noise = 0
                if add_noise:
                    noiseSigma = 0.05 * (reaction_flux / k)
                    noise = noiseSigma * np.random.randn()

                proteomics[protein_id] = abs((reaction_flux / k) + noise)

                # create transcriptomics dict
                noise = 0
                if add_noise:
                    noiseSigma = 0.05 * proteomics[protein_id] / q
                    noise = noiseSigma * np.random.randn()
                transcriptomics[gene.id] = abs((proteomics[protein_id] / q) + noise)

    return proteomics, transcriptomics


def get_metabolomics_data(model, solution, mapping_file, file_mapped=False):
    """
    Get metabolomics data

    :param model:
    :param solution:
    :mapping_file:

    :return:
    :metabolomics:
    :metabolomics_with_old_ids:
    """

    metabolomics = {}
    metabolomics_with_old_ids = {}
    # get metabolites

    # read the inchikey to pubchem ids mapping file
    if file_mapped:
        inchikey_to_cid = mapping_file
    else:
        inchikey_to_cid = {}
        inchikey_to_cid = read_pubchem_id_file(mapping_file)

    # create the stoichoimetry matrix fomr the model as a Dataframe and
    # onvert all the values to absolute values
    sm = create_stoichiometric_matrix(model, array_type="DataFrame")

    # get all the fluxes across reactions from the solution
    fluxes = solution.fluxes

    # calculating the dot product of the stoichiometry matrix and the fluxes
    # to calculate the net change
    # in concentration of the metabolites across reactions
    net_change_in_concentrations = sm.abs().dot(fluxes.abs())
    # net_change_in_concentrations = net_change_in_concentrations.abs()

    # converting all na values to zeroes and counting the total number of
    # changes that happens for each metabolite
    num_changes_in_metabolites = sm.fillna(0).astype(bool).sum(axis=1)

    for met_id, conc in net_change_in_concentrations.items():
        metabolite = model.metabolites.get_by_id(met_id)

        # if there is an inchikey ID for the metabolite
        if "inchi_key" in metabolite.annotation:
            # if it is a list get the first element
            if type(metabolite.annotation["inchi_key"]) is list:
                inchi_key = metabolite.annotation["inchi_key"][0]
            else:
                inchi_key = metabolite.annotation["inchi_key"]

            if inchi_key in inchikey_to_cid.keys():
                # if the CID is not in the metabolomics dict keys AND the
                # mapped value is not None and the reactions flux is not 0
                if (inchikey_to_cid[inchi_key] not in metabolomics.keys()) and (
                    inchikey_to_cid[inchi_key] is not None
                ):
                    metabolomics[inchikey_to_cid[inchi_key]] = (
                        conc
                        / num_changes_in_metabolites.iloc[
                            num_changes_in_metabolites.index.get_loc(met_id)
                        ]
                    )
                    metabolomics_with_old_ids[met_id] = (
                        conc
                        / num_changes_in_metabolites.iloc[
                            num_changes_in_metabolites.index.get_loc(met_id)
                        ]
                    )

                elif inchikey_to_cid[inchi_key] is not None:
                    metabolomics[inchikey_to_cid[inchi_key]] += (
                        conc
                        / num_changes_in_metabolites.iloc[
                            num_changes_in_metabolites.index.get_loc(met_id)
                        ]
                    )
                    metabolomics_with_old_ids[met_id] = (
                        conc
                        / num_changes_in_metabolites.iloc[
                            num_changes_in_metabolites.index.get_loc(met_id)
                        ]
                    )

    return metabolomics, metabolomics_with_old_ids


def get_multiomics(model, solution, mapping_file, old_ids=False):
    """
    Get mulitomics data

    :param model: cobra model object
    :param solution: solution for the model optimization using cobra
    :param mapping_file:
    :old_ids:

    :return:
    :proteomics:
    :transcriptomics:
    :metabolomics:
    """

    proteomics = {}
    transcriptomics = {}
    fluxomics = {}
    metabolomics = {}

    proteomics, transcriptomics = get_proteomics_transcriptomics_data(model, solution)

    metabolomics, metabolomics_with_old_ids = get_metabolomics_data(
        model, solution, mapping_file
    )

    if old_ids:
        return (proteomics, transcriptomics, metabolomics, metabolomics_with_old_ids)
    else:
        return (proteomics, transcriptomics, metabolomics)
