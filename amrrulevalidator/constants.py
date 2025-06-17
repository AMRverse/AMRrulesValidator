# Current version of the AMR rules specification
SPEC_VERSION = "v0.6"

# for the current spec version, this are the columns that are expected in the rules file
CANONICAL_COLUMNS = ["ruleID", "txid", "organism", "gene", "nodeID", "protein accession", "HMM accession", "nucleotide accession", "ARO accession", "mutation", "variation type", "gene context", "drug", "drug class", "phenotype", "clinical category", "breakpoint", "breakpoint standard", "breakpoint condition", "PMID", "evidence code", "evidence grade", "evidence limitations", "rule curation note"]

# the following variables are the allowed values for columns in the rules file, where these are not drawn from outside sources like NCBI or CARD
# NCBI/CARD values are included in the ResourceManager class, inside resources.py
PHENOTYPE = ['wildtype', 'nonwildtype']

GENE_CONTEXT = ['core', 'acquired']

CLINICAL_CAT = ["S", "I", "R"]
