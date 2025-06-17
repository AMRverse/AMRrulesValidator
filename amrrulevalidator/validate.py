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

    print("\nChecking breakpoint column...")
    summary_checks["breakpoint"], rows = check_breakpoint(
        get_column("breakpoint", rows), 
        get_column("clinical category", rows), 
        rows
    )

    # Check breakpoint standard
    print("\nChecking breakpoint standard column...")
    summary_checks["breakpoint standard"], rows = check_bp_standard(get_column("breakpoint standard", rows), rows)
    
    # Check breakpoint condition
    print("\nChecking breakpoint condition column...")
    summary_checks["breakpoint condition"], rows = check_bp_condition(get_column("breakpoint condition", rows), rows)

    # Check PMID
    print("\nChecking PMID column...")
    summary_checks["PMID"], rows = check_PMID(get_column("PMID", rows), rows)

    # Check evidence code
    print("\nChecking evidecnce code column...")
    summary_checks["evidence code"], rows = check_evidence_code(get_column("evidence code", rows), rows)
 

    # Check evidence grade and limitations
    print("\nChecking evidence grade and limitations columns...")
    summary_checks["evidence grade and limitations"], rows = check_evidence_grade_limitations(get_column("evidence grade", rows), get_column("evidence limitations", rows), rows)

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

    print("Writing output file...")
    # Write the processed rows to the output file
    write_tsv(rows, output_p, CANONICAL_COLUMNS)

    
    return True
