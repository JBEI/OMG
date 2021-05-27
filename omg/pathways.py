import cobra

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

def add_isopentenol_pathway(model, sce, user_params, write_file=True):
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
        model.add_reactions([r])
    
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
    if write_file:
        outputfilename = user_params['modelfile'].split('.')[0] + '_IPP.json'
        cobra.io.save_json_model(model, f'data/{outputfilename}')
    
    return model

