# Current version of the AMR rules specification
SPEC_VERSION = "v0.6"

# for the current spec version, this are the columns that are expected in the rules file
CANONICAL_COLUMNS = ["ruleID", "txid", "organism", "gene", "nodeID", "protein accession", "HMM accession", "nucleotide accession", "ARO accession", "mutation", "variation type", "gene context", "drug", "drug class", "phenotype", "clinical category", "breakpoint", "breakpoint standard", "breakpoint condition", "PMID", "evidence code", "evidence grade", "evidence limitations", "rule curation note"]

# the following variables are the allowed values for columns in the rules file, where these are not drawn from outside sources like NCBI or CARD
# NCBI/CARD values are included in the ResourceManager class, inside resources.py
PHENOTYPE = ['wildtype', 'nonwildtype']

GENE_CONTEXT = ['core', 'acquired']

CLINICAL_CAT = ["S", "I", "R"]

BREAKPOINT_CONDITIONS = [
            "-", "Endocarditis", "Endocarditis with combination treatment",
            "Intravenous", "Meningitis", "Meningitis, Endocarditis",
            "Non-endocarditis", "Non-meningitis", "Non-meningitis, Non-endocarditis",
            "Non-pneumonia", "Oral", "Oral, Infections originating from the urinary tract",
            "Oral, Other indications", "Oral, Uncomplicated urinary tract infection",
            "Pneumonia", "Prophylaxis", "Respiratory", "Screen", "Skin",
            "Uncomplicated urinary tract infection"
        ]

EVIDENCE_CODES = ["ECO:0001091 knockout phenotypic evidence", "ECO:0000012 functional complementation evidence", "ECO:0001113 point mutation phenotypic evidence", "ECO:0000024 protein-binding evidence", "ECO:0001034 crystallography evidence", "ECO:0000005 enzymatic activity assay evidence", "ECO:0000042 gain-of-function mutant phenotypic evidence", "ECO:0007000 high throughput mutant phenotypic evidence", "ECO:0001103 natural variation mutant evidence", "ECO:0005027 genetic transformation evidence", "ECO:0000020 protein inhibition evidence", "ECO:0006404 experimentally evolved mutant phenotypic evidence", "ECO:0000054 double mutant phenotype evidence"]

EVIDENCE_GRADES = ["high", "moderate", "low", "very low"]

EVIDENCE_LIMITATIONS = ["lacks evidence for this species", "lacks evidence for this genus", "lacks evidence for this allele", "lacks evidence of the degree to which MIC is affected", "low clinical relevance", "unknown clinical relevance", "statistical geno/pheno evidence but no experimental evidence", "conflicting evidence"]