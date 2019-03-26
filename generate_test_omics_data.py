
import cobra
import pandas as pd


def generate_fake_data(model, condition):
    """

    :param model: cobra model object
    :param solution: solution for the model optimization using cobra
    :param data_type: defines the type of -omics data to generate (all by default)
    :return:
    """

    # reaction_id of choice passed to the function# hardcoded here for this particular file (Need to convert this to an interactive cli program)
    reaction_id = 'BIOMASS_Ec_iJO1366_core_53p95M'
    # solution = None
    solution = get_optimized_solution(model, reaction_id)

    get_protemics_transcriptomics_data(model, solution)

    get_metabolomics_data(model, solution, condition)




def get_metabolomics_data(model, solution, condition):
    metabolomics = {}
    # get metabolites
    for met in model.metabolites:
        # get associated reactions
        for reaction in list(met.reactions):
            # get dictionary of associated metaboolites and their concentrations
            for metabolite, conc in reaction._metabolites.items():
                if metabolite.id == met.id:
                    if met.id not in metabolomics.keys():
                        metabolomics[met.id] = abs(conc)
                    else:
                        metabolomics[met.id] += abs(conc)
        # getting number of associated reactions and averaging the metabolic concentration value
        num_reactions = len(list(met.reactions))
        metabolomics[met.id]/=num_reactions

    metabolomics_dataframe = pd.DataFrame.from_dict(metabolomics, orient='index', columns=['metabolomics_value'])
    # Write the dataframe into a csv file
    file_name = 'metabolomics_fakedata.csv'
    metabolomics_dataframe.to_csv(file_name, sep=',', encoding='utf-8')



def get_protemics_transcriptomics_data(model, solution):

    # pre-determined linear constant (NOTE: Allow user to set this via parameter)
    # DISCUSS!!
    k = 0.8
    q = 0.6

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
                    protein_id = gene.annotation['goa']
                else:
                    protein_id = gene.annotation['uniprot']

                # create proteomics dict
                proteomics[protein_id] = solution.fluxes[rxnId]/k

            # create transcriptomics dict
            transcriptomics[gene.id] = proteomics[protein_id]/q


    fake_file_name = 'proteomics_fakedata.csv'
    proteomics_dataframe = pd.DataFrame.from_dict(proteomics, orient='index', columns=['proteomics_value'])
    # Write the dataframe into a csv file
    proteomics_dataframe.to_csv(fake_file_name, sep=',', encoding='utf-8')

    fake_file_name = 'transcriptomics_fakedata.csv'
    transcriptomics_dataframe = pd.DataFrame.from_dict(transcriptomics, orient='index', columns=['transcriptomics_value'])
    # Write the dataframe into a csv file
    transcriptomics_dataframe.to_csv(fake_file_name, sep=',', encoding='utf-8')





def get_list_of_reactions(file_name):
    """

    :param file_name: Name of the model file (has to be xml for now)
    :return: None (prints the list of reactions that has mass in them)
    """

    # Load model¶depending on the kind of file (the file has to be xml)
    if file_name.endswith(".xml"):
        model = cobra.io.read_sbml_model(file_name)

    # Print out the reaction name and reaction id for all reactions related to BIOMASS production:
    for rxn in model.reactions:
        if rxn.name is not None and 'BIOMASS' in rxn.id:
            print("{}: {}".format(rxn.id, rxn.name))



def get_optimized_solution(model, reaction_id):
    """

    :param model:
    :param reaction_id:
    :return:
    """

    # check to see if the reaction is bound (need to add a check here)
    # DISCUSS!!
    print(model.reactions.get_by_id(reaction_id).bounds)

    # fix the flux value to -15 as we have data for this constraint
    model.reactions.get_by_id(reaction_id).lower_bound = -15
    model.reactions.get_by_id(reaction_id).upper_bound = -15
    model.reactions.get_by_id(reaction_id)

    solution = model.optimize()

    return solution


def read_model(file_name):
    """

    :param file_name:
    :param reaction_id:
    :return:
    """

    # Load model¶depending on the kind of file
    if file_name.endswith(".xml"):
        model = cobra.io.read_sbml_model(file_name)
    elif file_name.endswith(".json"):
        model = cobra.io.load_json_model(file_name)

    return model




if __name__ == "__main__":
    # file_name = 'iJR904.json'
    file_name = 'data/iECIAI39_1322.xml'
    # reaction_id = 'EX_glc__D_e'

    # spits out the list of reaction names related to BIOMASS production
    get_list_of_reactions(file_name)

    # read model
    model = read_model(file_name)

    # get fake proteomics data and write it to file
    # NOTE: (need to add type of file, only CSV for now)
    condition = 1
    generate_fake_data(model, condition)


