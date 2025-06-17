"""Validation check functions extracted from legacy script."""

import csv
import re
from pathlib import Path
from amrrulevalidator.utils.check_helpers import report_check_results, validate_pattern, check_values_in_list, check_if_col_empty
from amrrulevalidator.constants import *



def check_if_not_missing(value_list, col_name, list_unique=False):

    invalid_indices = [index for index, value in enumerate(value_list) if value.strip() == '' or value.strip() in ['NA', '-']]

    if not invalid_indices:
        print("✅ All " + col_name + " values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print(col_name + " column must contain a value that is not NA or '-'")
        for index in invalid_indices:
            print(f"Row {index + 2}: {value_list[index]}")
    if list_unique:
        unique_values = set(value_list)
        print(f"\nUnique {col_name} values: {', '.join(map(str, unique_values))}")
    # now return check value
    if not invalid_indices:
        return True
    else:
        return False


def check_ruleIDs(id_list, rows):
    """
    Checks that rule IDs are unique and have the same prefix.
    Args:
        id_list: List of rule IDs to check.
        rows: list of row dictionaries to flag missing values in.
    Returns:
        tuple: A tuple containing:
            - A set of unique rule IDs.
            - A boolean indicating if the check passed.
            - The modified rows with any invalid rule IDs flagged.
    """

    invalid_dict = {}
    prefix_options = []
    valid_rule_ids = []

    # if all values are empty, fail the check
    ruleID_missing, rows = check_if_col_empty(id_list, 'ruleID', rows=rows)
    if ruleID_missing:
        return None, False, rows

    # Check for consistent prefix
    prefix_options = []
    for index, rule_id in enumerate(id_list):
        if rule_id.strip() == '' or rule_id.strip() in ['NA', '-']:
            # If the rule ID is empty, NA, or '-', flag it as missing
            rows[index]['ruleID'] = 'ENTRY MISSING'
            invalid_dict[index] = "Rule ID is empty, 'NA', or '-'"
            continue
        # otherwise, we have something in this cell, so first let's check that it starts with three capital letters
        # add a flag to this cell if that's the case
        if not re.match(r'^[A-Z]{3}', rule_id.strip()):
            invalid_dict[index] = f"Rule ID '{rule_id}' does not start with an expected prefix (should be three capital letters matching the organism)."
            rows[index]['ruleID'] = 'CHECK VALUE: ' + rule_id
            continue
        # if it does, we can check the prefix
        rule_prefix = rule_id.strip()[:3]  # Get the first three characters as prefix
        # if prefix_options is empty, we can set it as the expected prefix
        if not prefix_options:
            prefix_options.append(rule_prefix)
            expected_prefix = rule_prefix
        else:
           # otherwise, we should check that the prefix is consistent with the previous one
           if rule_prefix != expected_prefix:
               invalid_dict[index] = f"Inconsistent prefix: {rule_prefix} (expected {expected_prefix})"
               rows[index]['ruleID'] = 'CHECK VALUE: ' + rule_id
               # add the different prefix to the list
               prefix_options.append(rule_prefix)
               continue
        # if we get here, the ruleID is valid, so we can add it to our set of ruleIDs that we will output later
        # make a note though if the ruleID is a duplicate that we've already seen
        if rule_id.strip() in valid_rule_ids:
            invalid_dict[index] = f"Duplicate rule ID: {rule_id.strip()}"
            rows[index]['ruleID'] = 'CHECK VALUE: ' + rule_id
        else:
            # if it's not a duplicate, we can add it to the valid rule IDs
            valid_rule_ids.append(rule_id.strip())
    
    # Generate success/failure message
    success_message = f"All ruleIDs are valid (prefix: {expected_prefix})"
    failure_message = "Rule IDs must be unique and have the same prefix."
    
    # Report results
    check_result = report_check_results(
        check_name="ruleID",
        invalid_dict=invalid_dict,
        success_message=success_message,
        failure_message=failure_message
    )

    if len(prefix_options) > 1:
        print(f"\nMultiple rule ID prefixes found: {', '.join(prefix_options)}")

    return valid_rule_ids, check_result, rows


def check_txid(txid_list, rows, ncbi_organism_dict):
    """ Check that the taxonomic IDs are valid."""

    txid_missing, rows = check_if_col_empty(txid_list, 'txid', rows=rows)

    if txid_missing:
        print("❌ txid column is empty. Please provide values in the column to validate.")
        return False, False, rows

    if not ncbi_organism_dict:
        print("❌ NCBI organism dictionary is not provided. Cannot validate txids without it.")
        return True, False, rows
    
    # for each txid, check that its a value (so not empty, NA, or '-'), and that it's in the ncbi_organism_dict keys
    invalid_dict, rows = check_values_in_list(
        value_list=txid_list,
        col_name='txid',
        allowed_values=set(ncbi_organism_dict.keys()),
        missing_allowed=False,
        rows=rows,
        fail_reason = "is not a valid NCBI taxonomic ID"
    )

    check_result = report_check_results(
        check_name="txid",
        invalid_dict=invalid_dict,
        success_message="All txids are valid",
        failure_message="Txids must be present, not 'NA' or '-'. Txids should be in the NCBI taxonomy list, as per file resources/ncbi_taxonomy.tsv."
    )

    return True, check_result, rows


def check_organism(organism_list, rows, ncbi_organism_dict):

    org_missing, rows = check_if_col_empty(organism_list, 'organism', rows=rows)

    if org_missing:
        print("❌ Organism column is empty. Please provide values in the column to validate.")
        return False, False, rows
    
    if not ncbi_organism_dict:
        print("❌ NCBI organism dictionary is not provided. Cannot validate organism names without it.")
        return True, False, rows
    
    invalid_dict = {}

    # Check if the organism names are valid
    for index, organism in enumerate(organism_list):
        # first check that the organism name is not empty, NA, or '-'
        organism = organism.strip()
        if organism in ['NA', '-', '']:
            rows[index]['organism'] = 'ENTRY MISSING'
            organism_list[index] = 'ENTRY MISSING'
            continue
        # now check that the organism name starts with 's__'
        if not organism.startswith('s__'):
            rows[index]['organism'] = 'CHECK VALUE: ' + organism
            invalid_dict = {index: f"Organism name {organism} does not start with 's__'"}
            continue
        # if all those pass, now check if it's in the NCBI organism dictionary
        organism_name = organism.replace('s__', '', 1)  # Remove 's__' prefix for comparison
        if organism_name not in ncbi_organism_dict.values():
            rows[index]['organism'] = 'CHECK VALUE: ' + organism
            invalid_dict = {index: f"Organism name {organism} is not in the NCBI taxonomy list"}
            continue

    check_result = report_check_results(
        check_name="organism",
        invalid_dict=invalid_dict,
        success_message="All organisms are valid",
        failure_message="Organisms must be present, not 'NA' or '-'. Organisms should be in the NCBI taxonomy list, as per file resources/ncbi_taxonomy.tsv. Organisms should start with the prefix 's__'.",
        unique_values = set(organism_list)
    )

    return True, check_result, rows


def check_txid_organism(txid_list, organism_list, rows, ncbi_organism_dict):
  
    # this check will only occur if both columns are present
    invalid_dict = {}
    # for all valid taxids, check that the associated organism name matches the one in the list
    for index, txid in enumerate(txid_list):
        txid = txid.strip()
        if txid in ncbi_organism_dict:
            expected_organism = ncbi_organism_dict[txid]
            # we need to split the 's__' prefix from the organism name before checking
            current_organism = organism_list[index].strip().replace('s__', '', 1)
            if current_organism != expected_organism:
                invalid_dict[index] = f"Organism name {current_organism} does not match expected name {expected_organism} for taxid {txid}"
                rows[index]['organism'] = f'CHECK VALUE: {organism_list[index]}'
                rows[index]['txid'] = f'CHECK VALUE: {organism_list[index]}'
  
    success_message = "All txid-organism pairs are valid"
    failure_message = "Txid-organism pairs must match the what is in resources/ncbi_taxonomy.tsv."

    check_result = report_check_results(
        check_name="txid-organism",
        invalid_dict=invalid_dict,
        success_message=success_message,
        failure_message=failure_message
    )
    
    return check_result, rows


def check_gene(gene_list, rule_list, rows):
    
    # first check if the gene column is empty
    gene_missing, rows = check_if_col_empty(gene_list, 'gene', rows=rows)

    if gene_missing:
        print("❌ Gene column is empty. Please provide values in this column to validate.")
        return False, rows
    
    invalid_dict = {}
    
    # if gene isn't empty, first check if there are any empty values
    for index, gene in enumerate(gene_list):
        gene = gene.strip()
        if gene in ['NA', '-', '']:
            rows[index]['gene'] = 'ENTRY MISSING'
            invalid_dict[index] = "Gene is empty, 'NA', or '-'"

    # now we want to check for gene names that are actually combo rules - if there are any, we want to check that any rule IDs mentioned here are present in the file already
    # if there is a value in gene list that follows the format of three capital letters followed by a string of four numbers, this is one to compare against rule ids
    if rule_list:
        # we will use a regex to find the ruleIDs in the gene list
        pattern = re.compile(r'[A-Z]{3}\d{4}')
        for index, gene in enumerate(gene_list):
            if index not in invalid_dict.keys():
                matches = pattern.findall(gene)
                for match in matches:
                    if match not in rule_list:
                        invalid_dict[index] = f"ruleID {match} is not present in the list of rules"
                        rows[index]['gene'] = 'CHECK VALUE: ' + gene
                        break
    else:
        print("\nNo rule IDs available, skipping combinatorial rule check in gene column.")
    
    check_result = report_check_results(
        check_name="gene",
        invalid_dict=invalid_dict,
        success_message="All gene values are valid",
        failure_message="Gene column must contain a value that is not empty, NA or '-'. If the gene is a combinatorial rule, it must match an existing ruleID in the ruleID column.",
    )

    return check_result, rows


def check_id_accessions(nodeID_list, protein_list, nucleotide_list, hmm_list, variation_type_list, refseq_prot_accessions, refseq_nucl_accessions, refseq_node_ids, hmm_accessions, rows):
    
    nodeID_missing, rows = check_if_col_empty(nodeID_list, 'nodeID', rows=rows)
    protein_missing, rows = check_if_col_empty(protein_list, 'protein accession', rows=rows)
    nucleotide_missing, rows = check_if_col_empty(nucleotide_list, 'nucleotide accession', rows=rows)
    hmm_missing, rows = check_if_col_empty(hmm_list, 'HMM accession', rows=rows)

    # now check individual columns for allowable values
    # this function is actually multiple smaller checks
    # for each list, if any value is NA or empty, we should replace with ENTRY MISSING, as there should be a dash
    # secondly, if any value is not empty, a dash, or ENTRY MISSING, we should check if it's in the relevant accession list
    # finally, any rows that have all values empty, NA, '-' or ENTRY MISSING should be checked to see if they have variation type 'Combination'
    # this would make that row valid. Otherwise, the row is invalid

    invalid_node_dict = {}
    invalid_prot_dict = {}
    invalid_nucl_dict = {}
    invalid_hmm_dict = {}

    if not nodeID_missing:
        invalid_node_dict, rows = check_values_in_list(nodeID_list, refseq_node_ids, 'nodeID', rows, missing_allowed=True, fail_reason="is not a valid NCBI Reference Gene Hierarchy node ID")

    if not protein_missing:
        invalid_prot_dict, rows = check_values_in_list(protein_list, refseq_prot_accessions, 'protein accession', rows, missing_allowed=True, fail_reason="is not an NCBI Reference Gene Catalog protein accession")

    if not nucleotide_missing:
        invalid_nucl_dict, rows = check_values_in_list(nucleotide_list, refseq_nucl_accessions, 'nucleotide accession', rows, missing_allowed=True, fail_reason="is not an NCBI Reference Gene Catalog nucleotide accession")

    if not hmm_missing:
        invalid_hmm_dict, rows = check_values_in_list(hmm_list, hmm_accessions, 'HMM accession', rows, missing_allowed=True, fail_reason="is not an AMRFinderPlus HMM accession")

    # Check that in combination, at least one of these columns has a value
    invalid_combo_dict = {}
    for index, values in enumerate(zip(nodeID_list, protein_list, nucleotide_list, hmm_list)):
        values = [value.strip() for value in values]
        if all(value in ['NA', '-', 'ENTRY MISSING', ''] for value in values):
            # if all the values are empty, check if variation type is 'Combination' for this row
            # if variation type is 'Combination', then this is a valid value
            if variation_type_list[index].strip() == 'Combination':
                continue
            else:
                invalid_combo_dict[index] = "All ID accessions are empty, NA, or '-'. At least one of these columns must contain a valid accession value."

    check_result = report_check_results(
        check_name="ID accessions",
        invalid_dict={**invalid_node_dict, **invalid_prot_dict, **invalid_nucl_dict, **invalid_hmm_dict, **invalid_combo_dict},
        success_message="All ID accessions are valid",
        failure_message="At least one ID accession must be present, not 'NA', empty or '-'. Node IDs should be in the NCBI Reference Gene Hierarchy node ID list, protein and nucleotide accessions should be in the NCBI Reference Gene Catalog accession lists, and HMM accessions should be in the AMRFinderPlus HMM accession list.\nNOTE: If you have used an accession outside of those reference catalogs (e.g. your gene is not present in the AMRFinderPlus database), then this check will fail. Please double check those accessions exist."
    )

    return check_result, rows

def check_aro(aro_list, aro_terms, rows):

    aro_missing, rows = check_if_col_empty(aro_list, 'ARO accession', rows=rows)

    if aro_missing:
        print("❌ ARO accession column is empty. Please provide values in this column to validate.")
        return False, rows
    
    #check_values_in_list(value_list, allowed_values, col_name, rows, missing_allowed=False, fail_reason=None)
    invalid_dict, rows = check_values_in_list(
        value_list=aro_list,
        allowed_values=aro_terms,
        col_name='ARO accession',
        rows=rows,
        missing_allowed=True,
        fail_reason="is not a valid ARO accession"
    )

    check_result = report_check_results(
        check_name="ARO accession",
        invalid_dict=invalid_dict,
        success_message="All ARO accessions are valid",
        failure_message="ARO accession column must contain a value that is not empty, NA or '-'. ARO accessions should be in the ARO ontology list, as per file resources/aro_terms.tsv."
    )

    return check_result, rows
    

def check_variation(variation_list, rows):

    variation_missing, rows = check_if_col_empty(variation_list, 'variation type', rows=rows)
    if variation_missing:
        print("❌ Variation type column is empty. Please provide values in this column to validate.")
        return False, rows

    variation_allowed_types = [
            "Gene presence detected", "Protein variant detected", 
            "Nucleotide variant detected", "Promoter variant detected", 
            "Inactivating mutation detected", "Gene copy number variant detected", 
            "Nucleotide variant detected in multi-copy gene", 
            "Low frequency variant detected", "Combination"
        ]

    invalid_dict, rows = check_values_in_list(
        value_list=variation_list,
        allowed_values=variation_allowed_types,
        col_name='variation type',
        rows=rows,
        missing_allowed=False,
        fail_reason="is not a valid variation type"
    )

    check_result = report_check_results(
        check_name="variation type",
        invalid_dict=invalid_dict,
        success_message="All variation type values are valid",
        failure_message="Variation type column must contain a value that is not empty, NA or '-'. Variation should be one of the following types: " + ", ".join(variation_allowed_types) + "."
    )

    return check_result, rows


def check_mutation_variation(mutation_list, variation_list, rows):

    # okay, so if variation is entirely 'ENTRY MISSING', that means that this column wasn't present originally
    # therefore no point check if mutation matches, because it just won't
    if all(value.strip() == 'ENTRY MISSING' for value in variation_list):
        print("❌ Variation type column is empty. This column must have values in it to compare with the mutation column.")
        return False, rows
    
    # otherwise we can check. Note that if any variation_type is 'ENTRY MISSING', then we can't check if mutation is in concordance

    invalid_dict = {}

    for index, (mutation, variation) in enumerate(zip(mutation_list, variation_list)):
        reason = None
        mutation = mutation.strip()
        variation = variation.strip()
        if variation == 'ENTRY MISSING' or variation.startswith('CHECK VALUE'):
            reason = "Variation type is either missing or invalid. This column must have valid values in it to compare with the mutation column."
        elif variation == "Gene presence detected" and mutation != '-':
            reason = "Mutation must be '-' if variation type is 'Gene presence detected'. Mutation can never be empty or NA."
        elif variation == "Combination" and mutation != '-':
            reason = "Mutation must be '-' if variation type is 'Combination'. Mutation can never be empty or NA."
        elif variation == "Nucleotide variant detected" and not mutation.startswith("c."):
            reason = "Mutation must start with 'c.' if variation type is 'Nucleotide variant detected'"
        elif variation == "Protein variant detected" and not mutation.startswith("p."):
            reason = "Mutation must start with 'p.' if variation type is 'Protein variant detected'"
        elif variation == "Promoter variant detected" and not re.match(r"^c\.(-|\[-|\(-)", mutation):
            reason = "Mutation must start with 'c.-', 'c.(-', or 'c.[-' if variation type is 'Promoter variant detected'. The - symbol indicates the position before the start of the gene where the mutation occurs."
        elif variation == "Nucleotide variant detected in multi-copy gene" and not mutation.startswith("c."):
            reason = "Mutation must start with 'c.' if variation type is 'Nucleotide variant detected in multi-copy gene'"
        elif variation == "Gene copy number variant detected" and not re.match(r"^c\.\[\d+\]", mutation):
            reason = "Mutation must be in the format 'c.[X]' where X is any number if variation type is 'Gene copy number variant detected'"
        elif variation == "Low frequency variant detected" and not re.match(r"^(c\.|p\.)", mutation):
            reason = "Mutation must start with either 'c.' (for nucleotide variant) or 'p.' (protein variant) if variation type is 'Low frequency variant detected'"
        if reason:
            invalid_dict[index] = reason
            if mutation == '' or mutation == 'NA':
                rows[index]['mutation'] = 'ENTRY MISSING'
            else:
                rows[index]['mutation'] = 'CHECK VALUE: ' + mutation
    
    check_result = report_check_results(
        check_name="mutation variation",
        invalid_dict=invalid_dict,
        success_message="All mutation and variation type value pairs are valid",
        failure_message="Mutation must be concordant with the variation type. For example, if variation type is 'Gene presence detected', mutation must be '-'. If variation type is 'Nucleotide variant detected', mutation must start with 'c.'."
    )

    return check_result, rows

    
def check_context(context_list, variation_list, rows):
    # valid values are core or acquired
    # if variation_type_list isn't None, make sure we check that context is either
    #core or acquired if validation_type isn't 'Combination'
    
    context_missing, rows = check_if_col_empty(context_list, 'gene context', rows)

    if context_missing:
        print("❌ Gene context column is empty. Please provide values in this column to validate.")
        return False, rows
    
    # if context isn't completely empty, but variation type is, we need to check first if it's a valid value of core or acquired
    if all(value.strip() == 'ENTRY MISSING' for value in variation_list):
        invalid_dict, rows = check_values_in_list(context_list, GENE_CONTEXT, 'gene context', rows=rows, missing_allowed=False, fail_reason="must be either 'core' or 'acquired'")
    # otherwise we need to do a more complex check
    else:
        invalid_dict = {}
        for index, (context, variation) in enumerate(zip(context_list, variation_list)):
            context = context.strip()
            variation = variation.strip()
            if context in ['NA', '']:
                rows[index]['gene context'] = 'ENTRY MISSING'
                invalid_dict[index] = "Gene context is empty, 'NA', or '-'."
                continue
            if context not in GENE_CONTEXT and variation != 'Combination':
                reason = "Gene context must be 'core' or 'acquired' if variation type is not 'Combination'."
                invalid_dict[index] = reason
                rows[index]['gene context'] = 'CHECK VALUE: ' + context
            if context != '-' and variation == 'Combination':
                reason = 'If variation type is "Combination", gene context must be "-".'
                invalid_dict[index] = reason
                rows[index]['gene context'] = 'CHECK VALUE: ' + context

    check_result = report_check_results(
        check_name="gene context",
        invalid_dict=invalid_dict,
        success_message="All gene context values are valid",
        failure_message="Gene context must be either 'core' or 'acquired', and cannot be empty or NA. If variation type is 'Combination', gene context must be '-'."
    )

    return check_result, rows


def check_drug_drugclass(drug_list, drug_class_list, rows, rm=None):
    drug_missing, rows = check_if_col_empty(drug_list, 'drug name', rows)
    drug_class_missing, rows = check_if_col_empty(drug_class_list, 'drug class', rows)

    invalid_dict_combo = {}
    missing_values = ['NA', '-', '', 'ENTRY MISSING']

    if not drug_missing and not drug_class_missing:
        # check in combination
        for index, (drug, drug_class) in enumerate(zip(drug_list, drug_class_list)):
            drug = drug.strip()
            drug_class = drug_class.strip()
            if drug in missing_values and drug_class in missing_values:
                invalid_dict_combo[index] = "Both drug and drug class are empty, 'NA', or '-'. At least one of these columns must contain a valid CARD drug or drug class name."
                rows[index]['drug'] = 'ENTRY MISSING'
                rows[index]['drug class'] = 'ENTRY MISSING'
                continue

    # now check each column individually, against the allowed values
    invalid_dict_drug = {}
    invalid_dict_drug_class = {}
    if not drug_missing:
        invalid_dict_drug, rows = check_values_in_list(
            value_list=drug_list,
            allowed_values=rm.drug_names(),
            col_name='drug name',
            rows=rows,
            missing_allowed=True,
            fail_reason="is not a valid CARD drug name"
        )
    if not drug_class_missing:
        invalid_dict_drug_class, rows = check_values_in_list(
            value_list=drug_class_list,
            allowed_values=rm.drug_classes(),
            col_name='drug class',
            rows=rows,
            missing_allowed=True,
            fail_reason="is not a valid CARD drug class name"
        )

    check_result = report_check_results(
        check_name="drug and drug class",
        invalid_dict={**invalid_dict_drug, **invalid_dict_drug_class, **invalid_dict_combo},
        success_message="All drug and drug class values are valid",
        failure_message="Drug and drug class columns must contain a value that is not empty, NA or '-'. Drug names and classes should be in the CARD ontology list, as per file resources/card_drug_names.tsv."
    )

    return check_result, rows


def check_phenotype(phenotype_list, rows):

    phenotype_missing, rows = check_if_col_empty(phenotype_list, 'phenotype', rows=rows)

    if phenotype_missing:
        print("❌ Phenotype column is empty. Please provide values in this column to validate.")
        return False, rows
    
    # check that phenotype is one of the allowable values, wildtype or non wildtype
    invalid_dict, rows = check_values_in_list(
        value_list=phenotype_list,
        allowed_values=PHENOTYPE,
        col_name='phenotype',
        rows=rows
    )

    check_result = report_check_results(
        check_name="phenotype",
        invalid_dict=invalid_dict,
        success_message="All phenotype values are valid",
        failure_message="Phenotype values must be either 'wildtype' or 'non wildtype'."
    )

    return check_result, rows


def check_clinical(clinical_cat_list, rows):
    
    print("\nChecking clinical category column...")

    # Check if clinical category column is empty
    clinical_cat_missing, rows = check_if_col_empty(clinical_cat_list, 'clinical category', rows)
    if clinical_cat_missing:
        print("❌ Clinical category column is empty. Please provide values in this column to validate.")
        return False, rows

    invalid_dict, rows = check_values_in_list(
        value_list=clinical_cat_list,
        allowed_values=CLINICAL_CAT,
        col_name='clinical category',
        rows=rows,
        missing_allowed=False,
        fail_reason="must be one of the following: S, I, or R."
    )
    
    check_result = report_check_results(
        check_name="clinical category",
        invalid_dict=invalid_dict,
        success_message="All clinical category values are valid",
        failure_message="Clinical category must be one of the following: S, I, or R."
    )

    return check_result, rows


def check_breakpoint(breakpoint_list, clinical_cat_list, rows):

    breakpoint_missing, rows = check_if_col_empty(breakpoint_list, 'breakpoint', rows)

    if breakpoint_missing:
        print("❌ Breakpoint column is empty. Please provide values in this column to validate.")
        return False, rows
    
    # otherwise, we can check the column itself, and then in combo with clinical category if that column has at least some values
    # breakpoint can be not applicable, but it can't be empty, NA, or '-'
    invalid_dict_breakpoint = {}
    for index, breakpoint in enumerate(breakpoint_list):
        breakpoint = breakpoint.strip()
        if breakpoint in ['NA', '', '-']:
            rows[index]['breakpoint'] = 'ENTRY MISSING'
            invalid_dict_breakpoint[index] = "Breakpoint is empty, 'NA', or '-'. Breakpoint should contain information about the MIC or disk breakpoint, or 'not applicable' if no breakpoint is available."
            continue
        # if the breakpoint is not empty, we can check that it starts with one of the expected prefixes
        if not any(breakpoint.startswith(prefix) for prefix in ['MIC <', 'MIC <=', 'disk >', 'MIC >', 'MIC >=', 'disk <', 'not applicable']):
            invalid_dict_breakpoint[index] = f"Breakpoint '{breakpoint}' does not start with an expected prefix (MIC <, MIC <=, disk >, MIC >, MIC >=, disk <, not applicable)."
            rows[index]['breakpoint'] = 'CHECK VALUE: ' + breakpoint
            continue
        # now check if breakpoint matches what we'd expected for the clinical category, only if clinical category isn't completely empty
        if not all(value.strip() == 'ENTRY MISSING' for value in clinical_cat_list):
            category = clinical_cat_list[index].strip()
            reason = None
            if category == 'S' and not any(breakpoint.startswith(prefix) for prefix in ['MIC <', 'MIC <=', 'disk >', 'not applicable']):
                reason = "If clinical category is 'S', breakpoint should contain a value of 'MIC <', 'MIC <=', or 'disk >'. 'not applicable' is an allowed value if no breakpoint is available due to expected resistances."
            if category == 'R' and not any(breakpoint.startswith(prefix) for prefix in ['MIC >', 'MIC >=', 'disk <', 'not applicable']):
                reason = "If clinical category is 'R', breakpoint should contain a value of 'MIC >', 'MIC >=', or 'disk <'. 'not applicable' is an allowed value if no breakpoint is available due to expected resistances."
            if category == 'ENTRY MISSING' or category.startswith('CHECK VALUE'):
                reason = "Clinical category is missing or invalid, so breakpoint cannot be validated. Please provide clinical category to validate this row properly."
            if reason:
                invalid_dict_breakpoint[index] = reason
                rows[index]['breakpoint'] = 'CHECK VALUE: ' + breakpoint

    check_result = report_check_results(
        check_name="clinical category and breakpoint",
        invalid_dict=invalid_dict_breakpoint,
        success_message="All clinical category and breakpoint values are concordant",
        failure_message="Clinical category and breakpoint values must be compatible."
    )

    return check_result, rows


def check_bp_standard(breakpoint_standard_list, rows):

    bp_standard_missing, rows = check_if_col_empty(breakpoint_standard_list, 'breakpoint standard', rows)

    if bp_standard_missing:
        print("❌ Breakpoint standard column is empty. Please provide values in this column to validate.")
        return False, rows
    
    # Allowable patterns for breakpoint standards
    suggested_values = [
        r'^ECOFF \(\w+ \d{4}\)$',  # ECOFF (Month Year)
        r'^EUCAST .+ Expert Rules \(\w+ \d{4}\)$',  # EUCAST [organism] Expert Rules (Month year)
        r'^EUCAST Expected Resistant Phenotypes v\d+(\.\d+)? \(\d{4}\)$',  # EUCAST Expected Resistant Phenotypes version (year)
        r'^(EUCAST|CLSI)\s+v\d+(\.\d+)?\s+\(\d{4}\)$'  # EUCAST/CLSI version (year)
    ]

    invalid_dict, rows = validate_pattern(breakpoint_standard_list, suggested_values, rows, 'breakpoint standard', missing_allowed=False)
    unique_values = set(breakpoint_standard_list)
    
    failure_message = ("We check for the following formats: ECOFF (Month Year), "
                      "EUCAST [organism] Expert Rules (Month year), "
                      "EUCAST Expected Resistant Phenotypes vX (year), "
                      "or EUCAST/CLSI vX (year).")
    
    check_result = report_check_results(
        check_name="breakpoint standard",
        invalid_dict=invalid_dict,
        success_message="All breakpoint standard values match expected patterns.",
        failure_message=failure_message,
        unique_values=unique_values
    )

    return check_result, rows


def check_bp_condition(breakpoint_condition_list, rows):

    bp_condition_missing, rows = check_if_col_empty(breakpoint_condition_list, 'breakpoint condition', rows)

    if bp_condition_missing:
        print("❌ Breakpoint condition column is empty. Please provide values in this column to validate.")
        return False, rows

    invalid_dict, rows = check_values_in_list(
        value_list=breakpoint_condition_list,
        allowed_values=BREAKPOINT_CONDITIONS,
        col_name='breakpoint condition',
        rows=rows,
        missing_allowed=True,
        fail_reason="is not a valid breakpoint condition"
    )

    check_result = report_check_results(
        check_name="breakpoint condition",
        invalid_dict=invalid_dict,
        success_message="All breakpoint condition values are valid",
        failure_message="Breakpoint condition must be one of the following: " + ", ".join(BREAKPOINT_CONDITIONS) + "."
    )

    return check_result, rows


def check_PMID(pmid_list, rows):

    pmid_missing, rows = check_if_col_empty(pmid_list, 'PMID', rows)

    if pmid_missing:
        print("❌ PMID column is empty. Please provide values in this column to validate.")
        return False, rows
    
    # Check that PMID values are present and not empty, NA, or '-'
    invalid_dict = {}
    for index, pmid in enumerate(pmid_list):
        pmid = pmid.strip()
        if pmid in ['NA', '']:
            rows[index]['PMID'] = 'ENTRY MISSING'
            invalid_dict[index] = "PMID is empty or 'NA'."
            continue
        if pmid == "-":
            rows[index]['PMID'] = 'CHECK VALUE: ' + pmid
            invalid_dict[index] = "PMID is '-', however most rules should have a PMID associated with them."
            continue

    check_result = report_check_results(
        check_name="PMID",
        invalid_dict=invalid_dict,
        success_message="All PMIDs are valid",
        failure_message="PMID column must contain a value that is not empty, NA or '-'. PMIDs should be positive integers."
    )

    return check_result, rows


def check_evidence_code(evidence_code_list, rows):

    evidence_code_missing, rows = check_if_col_empty(evidence_code_list, 'evidence code', rows)

    if evidence_code_missing:
        print("❌ Evidence code column is empty. Please provide values in this column to validate.")
        return False, rows
    
    # check we have an evidence code that belongs to one of the allowable values
    # if it doesn't, check to see if the code starts with ECO:, and if it does, flag it so it can be checked
    
    
    # can be more than one of those values in this column, so need to split on the , separating them
    invalid_dict = {}
    invalid_codes = []
    
    for index, value in enumerate(evidence_code_list):
        value = value.strip()
        if value == '' or value in ['NA', '-']:
            invalid_dict[index] = "Evidence code must contain a value that is not empty, NA or '-'."
            rows[index]['evidence code'] = 'ENTRY MISSING'
            continue
        if ';' in value:
            invalid_dict[index] = "Evidence codes should be separated by a comma, not a semi-colon."
            rows[index]['evidence code'] = 'CHECK VALUE: ' + value
            continue
        codes = [code.strip() for code in value.split(',')]
        for code in codes:
            if not code.startswith("ECO:"):
                invalid_dict[index] = f"Evidence code '{code}' does not start with 'ECO:', and so may not be valid."
                rows[index]['evidence code'] = 'CHECK VALUE: ' + value
                continue
            if code not in EVIDENCE_CODES:
                invalid_dict[index] = f"Evidence code '{code}' is not in the list of typical evidence codes."
                rows[index]['evidence code'] = 'CHECK VALUE: ' + value
                # only add the code to the list if it's got an eco code prefix
                if code.startswith("ECO:"):
                    invalid_codes.append(code)

    check_result = report_check_results(
        check_name="evidence code",
        invalid_dict=invalid_dict,
        success_message="All evidence codes are valid",
        failure_message="Evidence code column must contain a value that is not empty, NA or '-'. Multiple evidence codes should be separated by a comma, not a semi-colon. Evidence codes should start with 'ECO:'. If the code is not in the list of typical evidence codes, it should be checked manually. Below is a list of unique invalid codes that were found:\n",
        unique_values=set(invalid_codes)
    )

    return check_result, rows


def check_evidence_grade_limitations(evidence_grade_list, evidence_limitations_list, rows):

    grade_missing, rows = check_if_col_empty(evidence_grade_list, 'evidence grade', rows=rows)
    limitations_missing, rows = check_if_col_empty(evidence_limitations_list, 'evidence limitations', rows=rows)

    # if both are missing, then return False
    if grade_missing and limitations_missing:
        print("❌ Evidence grade and limitations columns are empty. Please provide values in these columns to validate.")
        return False, rows
    
    # if grade isn't missing, it's a required col, so we can check that the value is in the list of allowable values
    invalid_grade_dict = {}
    if not grade_missing:
        invalid_grade_dict, rows = check_values_in_list(
            value_list=evidence_grade_list,
            allowed_values=EVIDENCE_GRADES,
            col_name='evidence grade',
            rows=rows,
            missing_allowed=False,
            fail_reason="is not a valid evidence grade"
        )
    # we can check limitations on its own, and make sure it's an allowed value
    invalid_limitations_dict = {}
    if not limitations_missing:
        # we need to split out the values one by one
        for index, value in enumerate(evidence_limitations_list):
            limitations = value.strip()
            if limitations in ['NA', '']:
                rows[index]['evidence limitations'] = 'ENTRY MISSING'
                invalid_limitations_dict[index] = "Evidence limitations is empty or 'NA', should be '-' if no value required"
                continue
            if limitations != '-':
                limitation_values = [lim.strip() for lim in limitations.split(',')]
                # check that all limitations are in the list of allowable limitations
                if not all(lim in EVIDENCE_LIMITATIONS for lim in limitation_values):
                    invalid_limitations_dict[index] = f"Evidence limitations '{limitations}' contains values that are not in the list of allowable limitations. Multiple limitations must be separated by a comma, or this check will fail."
                    rows[index]['evidence limitations'] = 'CHECK VALUE: ' + limitations
        
    # now check that if grade is 'moderate', 'low', or 'very low', then limitations is not empty
    if not grade_missing and not limitations_missing:
        for index, (grade, limitations) in enumerate(zip(evidence_grade_list, evidence_limitations_list)):
            grade = grade.strip()
            limitations = limitations.strip()
            # only check the row if we haven't flagged it as invalid already
            if not limitations.startswith('CHECK VALUE') and index not in invalid_grade_dict.keys():
                if grade in ["moderate", "low", "very low"] and limitations in ['-', 'ENTRY MISSING']:
                        invalid_limitations_dict[index] = f"If evidence grade is '{grade}', evidence limitations must contain a value, not '-' or 'ENTRY MISSING'."
                        rows[index]['evidence limitations'] = 'CHECK VALUE: ' + limitations

    check_result = report_check_results(
        check_name="evidence grade and limitations",
        invalid_dict={**invalid_grade_dict, **invalid_limitations_dict},
        success_message="All evidence grade and limitations values are valid",
        failure_message="Evidence grade must be one of the following: " + ", ".join(EVIDENCE_GRADES) + ". Evidence limitations must be one of the following: " + ", ".join(EVIDENCE_LIMITATIONS) + ". If evidence grade is 'moderate', 'low', or 'very low', then evidence limitations must contain a value, not '-'. If evidence limitations is empty, it should be '-', not 'ENTRY MISSING'."
    )
    
    return check_result, rows
