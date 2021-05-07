
import cobra
import pandas as pd
import random 
import os


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

    return inchikey_to_cid

def write_experiment_description_file(output_file_path, line_name='WT', label=''):
    
    # HARD CODED ONLY FOR WILD TYPE!
    
    if not os.path.isdir(output_file_path):
        os.mkdir(output_file_path)
            
    # create the filename
    experiment_description_file_name = f'{output_file_path}/EDD_experiment_description_file{label}.csv'

    #write experiment description file
    try:
        with open(experiment_description_file_name, 'w') as fh:
            fh.write(f'Line Name, Line Description, Part ID, Media, Shaking Speed, Starting OD, Culture Volume, Flask Volume, Growth Temperature, Replicate Count\n')
            if line_name == 'WT':
                line_descr = 'Wild type E. coli'
                part_id = 'ABFPUB_000310'
            else:
                line_descr = ''
                part_id = 'ABFPUB_000310'  #THIS SHOULD BE CHANGED!
            fh.write(f"{line_name}, {line_descr}, {part_id}, M9, 1, 0.1, 50, 200, 30, 1\n")
    except Exception as ex:
        print("Error in writing file!")
        print(ex)

def write_in_al_format(time_series_omics_data, omics_type, user_params, label=''):
    
    try:
        output_file_path = user_params['al_omics_file_path']
        if not os.path.isdir(output_file_path):
            os.mkdir(output_file_path)
            
        for timepoint, omics_dict in time_series_omics_data.items():
            al_file_name = f'{output_file_path}/AL_{omics_type}_{timepoint}_hrs{label}.csv'

            with open(al_file_name, 'w') as ofh:
                dataframe = pd.DataFrame.from_dict(omics_dict, orient='index', columns=[f'{omics_type}_value'])
                for index, series in dataframe.iteritems():
                    for id, value in series.iteritems():
                        ofh.write(f'{id},{value}\n')
    except:
        print('Error in writing in Arrowland format')
        
def write_in_edd_format(time_series_omics_data, omics_type, user_params, line_name, label=''):
    
    # Dictionary to map omics type with the units of measurement
    unit_dict = { "fluxomics": 'mmol/gdwh',\
        "proteomics": 'proteins/cell',\
        "transcriptomics": "FPKM",\
        "metabolomics": "mM"
    }
    
    # write in EDD format
    output_file_path = user_params['edd_omics_file_path']
    # create the filenames
    omics_file_name: str = f'{output_file_path}/EDD_{omics_type}{label}.csv'
            
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
                    
def write_omics_files(time_series_omics_data, omics_type, user_params, line_name='WT', al_format=False, label=''):
    """

    :param dataframe:
    :param data_type:
    :return:
    """
    
    # check which format we have to create the data in
    if not al_format:
        # write the omics files in EDD format by separating in terms of the timepoints
        write_in_edd_format(time_series_omics_data, omics_type, user_params, line_name, label=label)
    
    else:
        # write the omics files in ARROWLAND format by separating in terms of the timepoints
        write_in_al_format(time_series_omics_data, omics_type, user_params, label=label)
    
def write_OD_data(cell, output_file_path, line_name='WT', label=''):
    # create the filename
    OD_data_file: str = f'{output_file_path}/EDD_OD{label}.csv'
    if not os.path.isdir(output_file_path):
        os.mkdir(output_file_path)

    # write experiment description file
    try:
        with open(OD_data_file, 'w') as fh:
            fh.write(f'Line Name,Measurement Type,Time,Value,Units\n')
            for index, value in cell.items():
                fh.write((f'{line_name},Optical Density,{index},{value},n/a\n'))

    except Exception as ex:
        print("Error in writing OD file")
        print(ex)
    
def write_training_data_with_isopentenol(df, filename, output_file_path):
    filename = f'{output_file_path}/{filename}'
    df.to_csv(filename, header=True, index=False)

def write_external_metabolite(substrates, output_file_path, output_metabolites, line_name='WT', label=''):
    # create the filename      
    external_metabolites: str = f'{output_file_path}/EDD_external_metabolites{label}.csv'
    if not os.path.isdir(output_file_path):
        os.mkdir(output_file_path)
        
    # Table for metabolites to be exported
    glucose     = substrates.loc[:, 'glc__D_e']
    ammonium    = substrates.loc[:, 'nh4_e']
    isopentenol = substrates.loc[:, 'isoprenol_e']
    acetate     = substrates.loc[:, 'ac_e']
    formate     = substrates.loc[:, 'for_e']
    lactate     = substrates.loc[:, 'lac__D_e']
    ethanol     = substrates.loc[:, 'etoh_e']
     
    output_metabolites = {
        "5793": glucose, "12988": isopentenol, "175": acetate, "283": formate, "612": lactate, "702": ethanol}
    
    # Write file lines
    try:
        with open(external_metabolites,'w') as fh:
            # Top header
            fh.write(f'Line Name,Measurement Type,Time,Value,Units\n')
            # Metabolite lines
            for cid in output_metabolites:
                met = output_metabolites[cid]
                for index,value in met.items():
                    fh.write((f'{line_name},CID:{cid},{index},{value},mM\n'))
   
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
