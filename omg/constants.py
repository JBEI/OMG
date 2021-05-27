# Constants
UNIPROT_URL = """https://www.uniprot.org/uploadlists/"""
CTS_URL = """https://cts.fiehnlab.ucdavis.edu/rest/convert/"""
# HOST NAME
HOST_NAME: str = "ropacus"
# TODO: Move some constants to variables by program arguments
DATA_FILE_PATH: Filename = Filename("data")
# Output file path
OUTPUT_FILE_PATH: Filename = Filename("data/output")
# INCHIKEY_TO_CID_MAP_FILE_PATH: mapping file path to map inchikey to cids
INCHIKEY_TO_CID_MAP_FILE_PATH: Filename = Filename("mapping")
# MODEL_FILENAME: Filename = Filename('iECIAI39_1322.xml')  # E. coli
MODEL_FILENAME: Filename = Filename("reannotated_base_v3.sbml")  # R. opacus
MODEL_FILEPATH: Filename = Filename("")
# Training file name
TRAINING_FILE_NAME: Filename = Filename("")
TRAINING_FILE_PATH: Filename = Filename("")
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
REACTION_ID_ECOLI: str = "BIOMASS_Ec_iJO1366_core_53p95M"  # E. coli
REACTION_ID: str = "biomass_target"  # R. opacus
# REACTION_ID: str = 'SRC_C00185_e'  # R. opacus
GENE_IDS_DBS: List[str] = ["kegg.genes"]  # R. opacus
# GENE_IDS_DBS: List[str] = ['uniprot', 'goa', 'ncbigi']  # E. coli
UNITS: Dict[Omics, str] = {
    Omics.PROTEOMICS: "proteins/cell",
    Omics.TRANSCRIPTOMICS: "FPKM",
    Omics.METABOLOMICS: "mM",
}
# Fix the flux value to -15 as we have data for this constraint
LOWER_BOUND: int = -15
UPPER_BOUND: int = -15

# Internals
_EPS = np.finfo(np.double).eps
