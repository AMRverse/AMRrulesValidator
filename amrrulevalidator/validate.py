"""AMR rules validator module."""

from pathlib import Path

from amrrulevalidator.constants import CANONICAL_COLUMNS, SPEC_VERSION
from amrrulevalidator.utils.io import read_tsv, write_tsv
from amrrulevalidator.utils.resources import ResourceManager
from amrrulevalidator.checks import *


def get_column(column_name, rows):
    """Extract values for a specific column from a list of row dictionaries."""
    return [row[column_name] for row in rows if column_name in row]


def _normalize_blanks(rows):
    """
    Normalize blank/empty values in rows.
    
    Converts empty strings and "NA" to "-" in all cells.
    
    Args:
        rows: List of dictionaries containing row data
        
    Returns:
        List of dictionaries with normalized blank values
    """
    for row in rows:
        for key in row:
            # Convert empty strings and "NA" to "-"
            if row[key] == "" or row[key].upper() == "NA":
                row[key] = "-"
    
    return rows


def _flag_missing(rows):
    """
    Flag missing values in rows.
    
    Replaces "-" with "ENTRY MISSING" to highlight required fields
    that are missing data.
    
    Args:
        rows: List of dictionaries containing row data
        
    Returns:
        List of dictionaries with missing values flagged
    """
    for row in rows:
        for key in row:
            # Replace "-" with "ENTRY MISSING"
            if row[key] == "-":
                row[key] = "ENTRY MISSING"
    
    return rows


def run_validate(input_p: Path, output_p: Path, rm: ResourceManager) -> bool:
    """
    Validate an AMRrules file.
    
    Args:
        input_p: Path to the input TSV file
        output_p: Path to the output TSV file
        rm: ResourceManager instance for accessing reference data
        
    Returns:
        bool: True if validation was successful, False otherwise
    """

    # set up dict to capture which checks pass or fail
    summary_checks = {}
    
    print(f"\nValidating rules file: {input_p}")

    # Read rows from input file
    rows = read_tsv(input_p)
    
    # Ensure every row has the required column (fill with empty string if absent)
    print(f"\nChecking that all required columns for spec {SPEC_VERSION} are present...")
    found_columns = list(rows[0].keys()) if rows else []
    for column in CANONICAL_COLUMNS:
        if column not in found_columns:
            print(f"❌ {column} column not found in file.")
            print(f"Adding {column} column and filling with empty values, to enable validation to continue.")
    for row in rows:
        for col in CANONICAL_COLUMNS:
            if col not in row:
                row[col] = ""

    # Check ruleID
    print("\nChecking ruleID column...")
    rule_ids = get_column("ruleID", rows)
    rule_ids, summary_checks["ruleID"], rows = check_ruleIDs(rule_ids, rows)

    # Check txid and organism
    print("\nChecking txid column...")
    # parse the taxonomy file to get valid txids and organisms
    taxonomy_file_path = rm.dir / "ncbi_taxonomy.tsv"
    if not taxonomy_file_path.exists():
        print("❌ Cannot find NCBI taxonomy file. Run 'amrrule update-resources' to download it.")
        ncbi_organism_dict = None

    # Read the taxonomy file and store all the organisms and their txids in a dictionary
    with open(taxonomy_file_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        ncbi_organism_dict = {row['Accession']: row['Name'] for row in reader}
    
    txid_list = get_column("txid", rows)
    txid_not_empty, summary_checks["txid"], rows = check_txid(txid_list, rows, ncbi_organism_dict)
    
    print("\nChecking organism column...")
    organism_list = get_column("organism", rows)
    org_not_empty, summary_checks["organism"], rows = check_organism(organism_list, rows, ncbi_organism_dict)
    
    #Only check txid and organism together if both columns are not empty
    if txid_not_empty and org_not_empty and ncbi_organism_dict:
        print("\nChecking txid and organism are valid together...")
        summary_checks["txid and organism"], rows = check_txid_organism(txid_list, organism_list, rows, ncbi_organism_dict)
  
    # Check gene
    print("\nChecking gene column...")
    # we need ruleIDs to validate gene, so if this is None we can't fully validate gene
    summary_checks["gene"], rows = check_gene(get_column("gene", rows), rule_ids, rows)
    
    # Check gene accessions
    # okay so for these columns, at least one of them must have a value
    print("\nChecking nodeID, refseq accession, GenBank accession and HMM accession columns...")
        
    # Print placeholder message about database version
    print(f"\nChecking against AMRFinderPlus database version {rm.get_amrfp_db_version()}...")

    refseq_prot_accessions, refseq_nucl_accessions = rm.refseq_accessions()
    
    # Check accessions
    summary_checks["gene accessions"], rows = check_id_accessions(
        get_column("nodeID", rows), 
        get_column("protein accession", rows), 
        get_column("nucleotide accession", rows), 
        get_column("HMM accession", rows), 
        get_column("variation type", rows), 
        refseq_prot_accessions, refseq_nucl_accessions, rm.refseq_nodes(), rm.hmm_accessions(), rows
    )

    # Check ARO accession
    print("\nChecking ARO accession column...")
    aro_terms = rm.aro_terms()  # Get ARO terms from ResourceManager
    summary_checks["ARO accession"], rows = check_aro(get_column("ARO accession", rows), aro_terms, rows)

    # Check variation type
    print("\nChecking variation type column...")
    summary_checks["variation type"], rows = check_variation(get_column("variation type", rows), rows)



    # Check mutation and variation type compatibility
    print("\nChecking mutation and variation type columns are compatible...")
    summary_checks["variation type mutation concordance"], rows = check_mutation_variation(
        get_column("mutation", rows), 
        get_column("variation type", rows), 
        rows
    )

    print("\nChecking gene context column...")
    summary_checks["gene context"], rows = check_context(
        get_column("gene context", rows), 
        get_column("variation type", rows), 
        rows
    )

    print("\nChecking drug and drug class columns...")
    summary_checks["drug and drug class"], rows = check_drug_drugclass(
            get_column("drug", rows), 
            get_column("drug class", rows),
            rows,
            rm)

    # Check phenotype
    print("\nChecking phenotype column...")
    summary_checks["phenotype"], rows = check_phenotype(
        get_column("phenotype", rows),
        rows)

    # Check clinical category
    print("\nChecking clinical category column...")
    summary_checks["clinical category"], rows = check_clinical(get_column("clinical category", rows), rows)

    print("Writing output file...")
    # Write the processed rows to the output file
    write_tsv(rows, output_p, CANONICAL_COLUMNS)

    # Check breakpoint
    if "breakpoint" in found_columns:
        summary_checks["breakpoint"] = check_if_not_missing(
            get_column("breakpoint", rows), 
            "breakpoint"
        )
    else:
        print(f"\n❌ No breakpoint column found in file. Spec {SPEC_VERSION} requires this column to be present. Continuing to validate other columns...")
        summary_checks["breakpoint"] = False

    # Check clinical category and breakpoint concordance
    if "clinical category" in found_columns and "breakpoint" in found_columns:
        summary_checks["clinical category and breakpoint concordance"] = check_sir_breakpoint(
            get_column("clinical category", rows), 
            get_column("breakpoint", rows)
        )
    else:
        summary_checks["clinical category and breakpoint concordance"] = False

    # Check breakpoint standard
    if "breakpoint standard" in found_columns:
        summary_checks["breakpoint standard"] = check_bp_standard(
            get_column("breakpoint standard", rows)
        )
    else:
        print(f"\n❌ No breakpoint standard column found in file. Spec {SPEC_VERSION} requires this column to be present. Continuing to validate other columns...")
        summary_checks["breakpoint standard"] = False
    
    # Check breakpoint condition
    if "breakpoint condition" in found_columns:
        breakpoint_condition_list = [
            "-", "Endocarditis", "Endocarditis with combination treatment",
            "Intravenous", "Meningitis", "Meningitis, Endocarditis",
            "Non-endocarditis", "Non-meningitis", "Non-meningitis, Non-endocarditis",
            "Non-pneumonia", "Oral", "Oral, Infections originating from the urinary tract",
            "Oral, Other indications", "Oral, Uncomplicated urinary tract infection",
            "Pneumonia", "Prophylaxis", "Respiratory", "Screen", "Skin",
            "Uncomplicated urinary tract infection"
        ]
        summary_checks["breakpoint condition"] = check_if_allowed_value(
            get_column("breakpoint condition", rows), 
            "breakpoint condition", 
            breakpoint_condition_list, 
            missing_allowed=True
        )
    else:
        print(f"\n❌ No breakpoint condition column found in file. Spec {SPEC_VERSION} requires this column to be present. Continuing to validate other columns...")
        summary_checks["breakpoint condition"] = False

    # Check PMID
    if "PMID" in found_columns:
        summary_checks["PMID"] = check_if_not_missing(
            get_column("PMID", rows), 
            "PMID"
        )
    else:
        print(f"\n❌ No PMID column found in file. Spec {SPEC_VERSION} requires this column to be present. Continuing to validate other columns...")
        summary_checks["PMID"] = False

    # Check evidence code
    if "evidence code" in found_columns:
        summary_checks["evidence code"] = check_evidence_code(
            get_column("evidence code", rows)
        )
    else:
        print(f"\n❌ No evidence code column found in file. Spec {SPEC_VERSION} requires this column to be present. Continuing to validate other columns...")
        summary_checks["evidence code"] = False

    # Check evidence grade and limitations
    if "evidence grade" in found_columns and "evidence limitations" in found_columns:
        summary_checks["evidence grade and limitations"] = check_evidence_grade_limitations(
            get_column("evidence grade", rows), 
            get_column("evidence limitations", rows)
        )
    else:
        for column in ["evidence grade", "evidence limitations"]:
            if column not in found_columns:
                print(f"\n❌ {column} column not found in file.")
        print(f"\n❌ Both evidence grade and limitations columns required for spec {SPEC_VERSION}.")
        summary_checks["evidence grade and limitations"] = False

    # Print summary of checks
    passed_checks = [check for check, status in summary_checks.items() if status]
    failed_checks = [check for check, status in summary_checks.items() if not status]

    print("\nSummary of checks:")
    print(f"✅ Passed: {len(passed_checks)}")
    for check in passed_checks:
        print(f" - {check}")
    print(f"❌ Failed: {len(failed_checks)}")
    for check in failed_checks:
        print(f"  - {check}")
    print(f"Checked against AMRFinderPlus database version: {rm.get_amrfp_db_version()}")
    
    # Process data for output
    print("\nProcessing data for output...")
    
    # Normalize blank values (convert empty strings and "NA" to "-")
    rows = _normalize_blanks(rows)
    
    # Flag missing values (replace "-" with "ENTRY MISSING")
    rows = _flag_missing(rows)
    
    print("Writing output file...")
    # Write the processed rows to the output file
    write_tsv(rows, output_p, CANONICAL_COLUMNS)
    
    return True
