"""AMR rules validator module."""

from pathlib import Path

from amrrulesvalidator.constants import CANONICAL_COLUMNS
from amrrulesvalidator.utils.io import read_tsv, write_tsv
from amrrulesvalidator.utils.resources import ResourceManager
from amrrulesvalidator.checks import *


def get_column(column_name, rows):
    """Extract values for a specific column from a list of row dictionaries."""
    return [row[column_name] for row in rows if column_name in row]


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
    # Read rows from input file
    rows = read_tsv(input_p)
    
    # Ensure every row has every canonical column (fill with empty string)
    for row in rows:
        for col in CANONICAL_COLUMNS:
            if col not in row:
                row[col] = ""
    
    # set up dict to capture which checks pass or fail
    summary_checks = {}
    
    print(f"\nValidating rules file: {input_p}")

    # Print column headers
    print("\nChecking that all required columns for spec v0.6 are present...")
    found_columns = list(rows[0].keys()) if rows else []
    for column in CANONICAL_COLUMNS:
        if column not in found_columns:
            print(f"❌ {column} column not found in file.")
    print("\nContinuing to validate values in found columns...")

    # Check ruleID
    if "ruleID" in found_columns:
        rule_ids = get_column("ruleID", rows)
        summary_checks["ruleID"] = check_ruleIDs(rule_ids)
    else:
        print("\n❌ No ruleID column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
        rule_ids = None
        summary_checks["ruleID"] = False

    # Check txid and organism
    if "txid" in found_columns and "organism" in found_columns:
        txid_list = get_column("txid", rows)
        organism_list = get_column("organism", rows)
        summary_checks["txid and organism"] = check_organism(txid_list, organism_list, rm)
    else:
        print("\n❌ No organism or txid column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
        summary_checks["txid and organism"] = False

    # Check gene
    if "gene" in found_columns and "ruleID" in found_columns:
        summary_checks["gene"] = check_gene(get_column("gene", rows), rule_ids)
    elif "gene" in found_columns and not "ruleID" in found_columns:
        print("\n❌ No ruleID column found in file. Spec v0.6 requires this column to be present, and cannot validate gene without it. Continuing to validate other columns...")
        summary_checks["gene"] = False
    else:
        print("\n❌ No gene column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
        summary_checks["gene"] = False
    
    # Check gene accessions
    accession_columns = ["nodeID", "protein accession", "nucleotide accession", "HMM accession"]
    if all(col in found_columns for col in accession_columns):
        refseq_url = None  # Replace with ResourceManager call when implemented
        hmm_url = None  # Replace with ResourceManager call when implemented
        amrfp_nodes_url = None  # Replace with ResourceManager call when implemented
        
        # Print placeholder message about database version
        print(f"\nChecking against AMRFinderPlus database version (placeholder)...")
        
        # Check accessions
        if "variation type" in found_columns:
            summary_checks["gene accessions"] = check_id_accessions(
                get_column("nodeID", rows), 
                get_column("protein accession", rows), 
                get_column("nucleotide accession", rows), 
                get_column("HMM accession", rows), 
                get_column("variation type", rows), 
                refseq_url, amrfp_nodes_url, hmm_url
            )
        else:
            summary_checks["gene accessions"] = check_id_accessions(
                get_column("nodeID", rows), 
                get_column("protein accession", rows), 
                get_column("nucleotide accession", rows), 
                get_column("HMM accession", rows), 
                None, refseq_url, amrfp_nodes_url, hmm_url
            )
    else:
        for column in accession_columns:
            if column not in found_columns:
                print(f"\n❌ {column} column not found in file.")
        print("\n❌ Spec v0.6 requires all of nodeID, protein accession, nucleotide accession, HMM accession and variation type columns to be present in order to validate. Continuing to validate other columns...")
        summary_checks["gene accessions"] = False

    # Check ARO accession
    if "ARO accession" in found_columns:
        aro_terms = None  # Replace with ResourceManager call when implemented
        summary_checks["ARO accession"] = check_aro(get_column("ARO accession", rows), aro_terms)
    else:
        print("\n❌ No ARO accession column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
        summary_checks["ARO accession"] = False
    
    # Check mutation
    if "mutation" in found_columns:
        summary_checks["mutation"] = check_mutation(get_column("mutation", rows))
    else:
        print("\n❌ No mutation column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
        summary_checks["mutation"] = False
    
    # Check variation type
    if "variation type" in found_columns:
        variation_allowed_types = [
            "Gene presence detected", "Protein variant detected", 
            "Nucleotide variant detected", "Promoter variant detected", 
            "Inactivating mutation detected", "Gene copy number variant detected", 
            "Nucleotide variant detected in multi-copy gene", 
            "Low frequency variant detected", "Combination"
        ]
        summary_checks["variation type"] = check_if_allowed_value(
            get_column("variation type", rows), 
            "variation type", 
            variation_allowed_types
        )
    else:
        print("\n❌ No variation type column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
        summary_checks["variation type"] = False

    # Check mutation and variation type compatibility
    if "variation type" in found_columns and "mutation" in found_columns:
        summary_checks["variation type mutation concordance"] = check_mutation_variation(
            get_column("mutation", rows), 
            get_column("variation type", rows)
        )
    else:
        summary_checks["variation type mutation concordance"] = False

    # Check gene context
    if "gene context" in found_columns:
        if "variation type" not in found_columns:
            summary_checks["gene context"] = check_context(
                get_column("gene context", rows), 
                None
            )
        else:
            summary_checks["gene context"] = check_context(
                get_column("gene context", rows), 
                get_column("variation type", rows)
            )
        
        # Check context and mutation concordance
        if "mutation" in found_columns:
            summary_checks["context and mutation concordance"] = check_context_mutation(
                get_column("mutation", rows), 
                get_column("gene context", rows)
            )
    else:
        print("\n❌ No gene context column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
        summary_checks["gene context"] = False

    # Check drug and drug class
    if "drug" in found_columns and "drug class" in found_columns:
        summary_checks["drug and drug class"] = check_drug_drugclass(
            get_column("drug", rows), 
            get_column("drug class", rows)
        )
    else:
        for column in ["drug", "drug class"]:
            if column not in found_columns:
                print(f"\n❌ {column} column not found in file.")
        print("\n❌ Spec v0.6 requires at least both drug and drug class columns to be present in order to validate. Continuing to validate other columns...")
        summary_checks["drug and drug class"] = False

    # Check phenotype
    if "phenotype" in found_columns:
        summary_checks["phenotype"] = check_if_allowed_value(
            get_column("phenotype", rows), 
            "phenotype", 
            ["wildtype", "nonwildtype"]
        )
        
        # Check phenotype and context concordance
        if "gene context" in found_columns:
            summary_checks["phenotype and context concordance"] = check_phenotype_context(
                get_column("phenotype", rows), 
                get_column("gene context", rows)
            )
    else:
        print("\n❌ No phenotype column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
        summary_checks["phenotype"] = False

    # Check clinical category
    if "clinical category" in found_columns:
        summary_checks["clinical category"] = check_if_allowed_value(
            get_column("clinical category", rows), 
            "clinical category", 
            ["S", "I", "R"]
        )
    else:
        print("\n❌ No clinical category column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
        summary_checks["clinical category"] = False

    # Check breakpoint
    if "breakpoint" in found_columns:
        summary_checks["breakpoint"] = check_if_not_missing(
            get_column("breakpoint", rows), 
            "breakpoint"
        )
    else:
        print("\n❌ No breakpoint column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
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
        print("\n❌ No breakpoint standard column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
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
        print("\n❌ No breakpoint condition column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
        summary_checks["breakpoint condition"] = False

    # Check PMID
    if "PMID" in found_columns:
        summary_checks["PMID"] = check_if_not_missing(
            get_column("PMID", rows), 
            "PMID"
        )
    else:
        print("\n❌ No PMID column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
        summary_checks["PMID"] = False

    # Check evidence code
    if "evidence code" in found_columns:
        summary_checks["evidence code"] = check_evidence_code(
            get_column("evidence code", rows)
        )
    else:
        print("\n❌ No evidence code column found in file. Spec v0.6 requires this column to be present. Continuing to validate other columns...")
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
        print("\n❌ Both evidence grade and limitations columns required for spec v0.6.")
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
    print("Checked against AMRFinderPlus database version: [placeholder]")

    # Write the processed rows to the output file
    write_tsv(rows, output_p, CANONICAL_COLUMNS)
    
    return True
