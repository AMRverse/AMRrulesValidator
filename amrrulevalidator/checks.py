"""Validation check functions extracted from legacy script."""

import csv
import re
from pathlib import Path
from amrrulevalidator.utils.check_helpers import report_check_results, validate_pattern, check_values_in_list, check_if_col_empty



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


def check_id_accessions(nodeID_list, protein_list, nucleotide_list, hmm_list, variation_type_list, refseq_file, node_file, hmm_file, rows):
    
    # parse the refseq file - get the refseq nucl and prot accessions, genbank accessions and hmm accessions to check against
    # Combine the nucl and prot accessions together
    protein_accessions = [] 
    nucleotide_accessions = []
    refseq = csv.DictReader(open(refseq_file, 'r'), delimiter='\t')
    for row in refseq:
        if "refseq_protein_accession" in row:
            protein_accessions.append(row["refseq_protein_accession"])
        if "refseq_nucleotide_accession" in row:
            nucleotide_accessions.append(row["refseq_nucleotide_accession"])
        if "genbank_protein_accession" in row:
            protein_accessions.append(row["genbank_protein_accession"])
        if "genbank_nucleotide_accession" in row:
            nucleotide_accessions.append(row["genbank_nucleotide_accession"])
    # remove any empty strings
    protein_accessions = [value for value in protein_accessions if value != ""]
    nucleotide_accessions = [value for value in nucleotide_accessions if value != ""]

    refseq_node_ids = []
    refseq_hierarchy = csv.DictReader(open(node_file, 'r'), delimiter='\t')
    for row in refseq_hierarchy:
        if "parent_node_id" in row:
            refseq_node_ids.append(row["parent_node_id"])
        if "node_id" in row:
            refseq_node_ids.append(row["node_id"])
    # remove any duplicates and empty strings
    refseq_node_ids = set(refseq_node_ids)
    refseq_node_ids = [value for value in refseq_node_ids if value != ""]

    hmm_accessions = []
    hmm_table = csv.DictReader(open(hmm_file, newline=''), delimiter='\t')
    for row in hmm_table:
        if "#hmm_accession" in row:
            hmm_accessions.append(row["#hmm_accession"])
    # remove any empty strings
    hmm_accessions  = [value for value in hmm_accessions if value != ""]

    # now check individual columns for allowable values
    # this function is actually multiple smaller checks
    # for each list, if any value is NA or empty, we should replace with ENTRY MISSING, as there should be a dash
    # secondly, if any value is not empty, a dash, or ENTRY MISSING, we should check if it's in the relevant accession list
    # finally, any rows that have all values empty, NA, '-' or ENTRY MISSING should be checked to see if they have variation type 'Combination'
    # this would make that row valid. Otherwise, the row is invalid
    
    invalid_node_dict, rows = check_values_in_list(nodeID_list, refseq_node_ids, 'nodeID', rows, missing_allowed=True, fail_reason="is not a valid NCBI Reference Gene Hierarchy node ID")

    invalid_prot_dict, rows = check_values_in_list(protein_list, protein_accessions, 'protein accession', rows, missing_allowed=True, fail_reason="is not an NCBI Reference Gene Catalog protein accession")

    invalid_nucl_dict, rows = check_values_in_list(nucleotide_list, nucleotide_accessions, 'nucleotide accession', rows, missing_allowed=True, fail_reason="is not an NCBI Reference Gene Catalog nucleotide accession")

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

def check_aro(aro_list, aro_terms):
    
    print("\nChecking ARO accession column...")

    # value is invalid if doesn't start with ARO or isn't '-'
    # need to expand this as well to check if the ARO accessions in the list are actually ones that exist in the CARD ontology.
    invalid_indices = []
    for index, value in enumerate(aro_list):
        value = value.strip()
        # a dash is fine
        if value == '-':
            continue
        # if the value is an accession, check it
        if value.startswith("ARO:"):
            if value not in aro_terms:
                invalid_indices.append(index)
        # the cell can't be empty, if it's empty it must be '-' (checked above)
        elif value == '' or value == 'NA':
            invalid_indices.append(index)
        # if it's anything else then it's invalid
        else:
            invalid_indices.append(index)

    if not invalid_indices:
        print("✅ All ARO accession values are valid and exist in the CARD ontology")
        return True
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("ARO accession column must contain a valid ARO accession, starting with 'ARO:', and cannot be empty. The following rows contain invalid or empty accessions:")
        for index in invalid_indices:
            print(f"Row {index + 2}: {aro_list[index]}")
        return False


def check_context(context_list, variation_type_list):
    # valid values are core or acquired
    # if variation_type_list isn't None, make sure we check that context is either
    #core or acquired if validation_type isn't 'Combination'
    print("\nChecking gene context column...")

    invalid_indices = {}
    if variation_type_list is not None:
        for index, (context, variation) in enumerate(zip(context_list, variation_type_list)):
            context = context.strip()
            variation = variation.strip()
            if context not in ['core', 'acquired'] and variation != 'Combination':
                reason = "Gene context must be 'core' or 'acquired' if variation type is not 'Combination'."
                invalid_indices[index] = reason
            if context != '-' and variation == 'Combination':
                reason = 'If variation type is "Combination", gene context must be "-".'
                invalid_indices[index] = reason
    if not variation_type_list:
        for index, context in enumerate(context_list):
            context = context.strip()
            if context not in ['core', 'acquired']:
                reason = "Gene context must be 'core' or 'acquired'."
                invalid_indices[index] = reason
    
    if not invalid_indices:
        print("✅ All gene context values are valid")
        return True
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Gene context column must contain either 'core' or 'acquired' and cannot be empty. If variation type is 'Combination', context must be '-'.")
        for index in invalid_indices:
            print(f"Row {index + 2}: {context_list[index]}; {invalid_indices[index]}")
        return False


def check_mutation(mutation_list):
    
    print("\nChecking mutation column...")
    
    # check that there is either a value or '-' in this column
    invalid_indices = [index for index, value in enumerate(mutation_list) if value.strip() == '']

    if not invalid_indices:
        print("✅ All mutation values are valid")
        return True
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Mutation column must contain either a value or '-' if no mutation required.")
        for index in invalid_indices:
            print(f"Row {index + 2}: {mutation_list[index]}")
        return False


def check_context_mutation(context_list, mutation_list):
    # check that if context is core, and mutation is not '-', then provide a warning
    invalid_indices = []
    for index, (context, mutation) in enumerate(zip(context_list, mutation_list)):
        context = context.strip()
        mutation = mutation.strip()
        if context == 'core' and mutation != '-':
            invalid_indices.append(index)
    if not invalid_indices:
        print("✅ All context and mutation values are concordant")
        return True
    else:
        print(f"⚠️ Warning: {len(invalid_indices)} rows may have incorrect values.")
        print("If the gene context is a 'core' gene, the expected mutation should generally be '-' unless the rule " \
        "refers to a specific variant of the core gene for which there is evidence of a nonwildtype mutation. "
        "Note that a resistance-associated mutation in a core gene (e.g. Ser83Phe in GyrA) should be coded as context 'acquired', " \
        "since the variant the rule applies to is not itself 'core'.")
        for index in invalid_indices:
            print(f"Row {index + 2}: {context_list[index]} and {mutation_list[index]}")
        return False


def check_mutation_variation(mutation_list, variation_list):
    
    # check that the mutation and variation type are compatible
    print("\nChecking mutation and variation type columns are compatible...")

    invalid_indices_dict = {}

    for index, (mutation, variation) in enumerate(zip(mutation_list, variation_list)):
        reason = None
        mutation = mutation.strip()
        variation = variation.strip()
        if variation == "Gene presence detected" and mutation != '-':
            reason = "Mutation must be '-' if variation type is 'Gene presence detected'"
        elif variation == "Combination" and mutation != '-':
            reason = "Mutation must be '-' if variation type is 'Combination'"
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
            invalid_indices_dict[index + 2] = reason

    if not invalid_indices_dict:
        print("✅ All mutation and variation type values are compatible")
        return True
    else:
        print(f"❌ {len(invalid_indices_dict)} rows have failed the check")
        for index, reason in invalid_indices_dict.items():
            print(f"Row {index}: {reason}")
        return False


def extract_card_drug_names(resource_manager=None):
    # read in the file that lists all the drugs and drug classes that are in the current version of the CARD ontology
    drug_names_card = []
    drug_classes_card = []
    
    if resource_manager:
        # Use the ResourceManager to access the file
        drug_names_file_path = resource_manager.dir / "card_drug_names.tsv"
        if not drug_names_file_path.exists():
            print("❌ Cannot find CARD drug names file. Run 'amrrule update-resources' to download it.")
            return [], []
    else:
        # Fallback to direct file access if no ResourceManager provided
        drug_names_file_path = Path('amrrulevalidator/resources/card_drug_names.tsv')
        if not drug_names_file_path.exists():
            drug_names_file_path = Path('card_drug_names.tsv')
            if not drug_names_file_path.exists():
                print("❌ Cannot find CARD drug names file. Run 'amrrule update-resources' to download it.")
                return [], []
    
    with open(drug_names_file_path, newline='') as card_drugs_file:
        reader = csv.DictReader(card_drugs_file, delimiter='\t')
        for row in reader:
            if row['Drug Name'] not in drug_names_card:
                drug_names_card.append(row['Drug Name'])
            if row['Drug Class'] not in drug_classes_card:
                drug_classes_card.append(row['Drug Class'])
    return drug_names_card, drug_classes_card


def check_drug_drugclass(drug_list, drug_class_list, resource_manager=None):
   
    print("\nChecking drug and drug class columns...")

    # read in the file that lists all the drugs and drug classes that are in the current version of the CARD ontology
    card_drugs, card_drug_classes = extract_card_drug_names(resource_manager)
    
    # want to check that there is at least one value in either drug or drug class
    # need to check that if the value isn't '-', its a valid drug or drug class name as per card_drugs and card_drug_classes
    invalid_indices_dict = {}
    for index, (drug, drug_class) in enumerate(zip(drug_list, drug_class_list)):
        drug = drug.strip()
        drug_class = drug_class.strip()
        if (drug == '' or drug == '-') and (drug_class == '' or drug_class == '-'):
            invalid_indices_dict[index] = "Both drug and drug class are empty. At least one of these columns must contain a valid CARD drug or drug class name."
            continue
        if drug != '-' and drug not in card_drugs:
            reason = 'Drug name ' + drug + ' is not a valid CARD drug name.'
            invalid_indices_dict[index] = reason
        if drug_class != '-' and drug_class not in card_drug_classes:
            reason = 'Drug class ' + drug_class + ' is not a valid CARD drug class name.'
            invalid_indices_dict[index] = reason
    
    if not invalid_indices_dict:
        print("✅ All drug and drug class values are valid and listed in the CARD drug name ontology.")
        return True
    else:
        print(f"❌ {len(invalid_indices_dict)} rows have failed the check")
        print("One of drug or drug class must contain a value that is not empty, NA or '-'. Values must be listed in the CARD drug name ontology, as per card_drug_names.tsv. Drugs and their classes should be given in all lower case.")
        for index, reason in invalid_indices_dict.items():
            print(f"Row {index + 2}: {reason}")
        return False


def check_phenotype_context(phenotype_list, context_list):
    
    print("\nChecking phenotype and context columns are concordant...")

    # check that if context is core, phenotype must be wildtype
    invalid_indices = []
    for index, (phenotype, context) in enumerate(zip(phenotype_list, context_list)):
        phenotype = phenotype.strip()
        context = context.strip()
        if context == 'core' and phenotype != 'wildtype':
            invalid_indices.append(index)
    
    if not invalid_indices:
        print("✅ All phenotype and context values are concordant")
        return True
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("If the gene context is a 'core' gene, the expected phenotype should generally be 'wildtype', " \
        "unless the rule refers to a specific variant of the core gene for which there is evidence of a nonwildtype " \
        "phenotype (in which case the variant should be coded as 'acquired' not core)")
        for index in invalid_indices:
            print(f"Row {index + 2}: {phenotype_list[index]} and {context_list[index]}")
        return False


def check_sir_breakpoint(clinical_category_list, breakpoint_list):
    """
    Check that clinical category (S/I/R) and breakpoint values are compatible.
    
    Args:
        clinical_category_list: List of clinical category values (S, I, R)
        breakpoint_list: List of breakpoint values
        
    Returns:
        bool: True if check passed, False otherwise
    """
    
    print("\nChecking clinical category and breakpoint columns...")

    invalid_indices_dict = {}

    for index, (category, breakpoint) in enumerate(zip(clinical_category_list, breakpoint_list)):
        reason = None
        category = category.strip()
        breakpoint = breakpoint.strip()
        
        if category == 'S' and not any(breakpoint.startswith(prefix) for prefix in ['MIC <', 'MIC <=', 'disk >', 'not applicable']):
            reason = "If clinical category is 'S', breakpoint should contain a value of 'MIC <', 'MIC <=', or 'disk >'. 'not applicable' is an allowed value if no breakpoint is available due to expected resistances."
        if category == 'R' and not any(breakpoint.startswith(prefix) for prefix in ['MIC >', 'MIC >=', 'disk <', 'not applicable']):
            reason = "If clinical category is 'R', breakpoint should contain a value of 'MIC >', 'MIC >=', or 'disk <'. 'not applicable' is an allowed value if no breakpoint is available due to expected resistances."
            
        if reason:
            invalid_indices_dict[index + 2] = reason

    return report_check_results(
        check_name="clinical category and breakpoint",
        invalid_dict=invalid_indices_dict,
        success_message="All clinical category and breakpoint values are concordant",
        failure_message="Clinical category and breakpoint values must be compatible."
    )


def check_bp_standard(breakpoint_standard_list):
    """
    Check that breakpoint standard values match expected patterns.
    
    Args:
        breakpoint_standard_list: List of breakpoint standard values to check
        
    Returns:
        bool: True if check passed, False otherwise
    """
    
    print("\nChecking breakpoint standard column...")
    
    # Allowable patterns for breakpoint standards
    suggested_values = [
        r'^ECOFF \(\w+ \d{4}\)$',  # ECOFF (Month Year)
        r'^EUCAST .+ Expert Rules \(\w+ \d{4}\)$',  # EUCAST [organism] Expert Rules (Month year)
        r'^EUCAST Expected Resistant Phenotypes v\d+(\.\d+)? \(\d{4}\)$',  # EUCAST Expected Resistant Phenotypes version (year)
        r'^(EUCAST|CLSI)\s+v\d+(\.\d+)?\s+\(\d{4}\)$'  # EUCAST/CLSI version (year)
    ]
    
    invalid_dict = validate_pattern(breakpoint_standard_list, suggested_values)
    unique_values = set(breakpoint_standard_list)
    
    failure_message = ("We check for the following formats: ECOFF (Month Year), "
                      "EUCAST [organism] Expert Rules (Month year), "
                      "EUCAST Expected Resistant Phenotypes vX (year), "
                      "or EUCAST/CLSI vX (year).")
    
    return report_check_results(
        check_name="breakpoint standard",
        invalid_dict=invalid_dict,
        success_message="All breakpoint standard values match expected patterns.",
        failure_message=failure_message,
        unique_values=unique_values
    )


def check_evidence_code(evidence_code_list):
    
    print("\nChecking evidence code column...")

    allowable_values = ["ECO:0001091 knockout phenotypic evidence", "ECO:0000012 functional complementation evidence", "ECO:0001113 point mutation phenotypic evidence", "ECO:0000024 protein-binding evidence", "ECO:0001034 crystallography evidence", "ECO:0000005 enzymatic activity assay evidence", "ECO:0000042 gain-of-function mutant phenotypic evidence", "ECO:0007000 high throughput mutant phenotypic evidence", "ECO:0001103 natural variation mutant evidence", "ECO:0005027 genetic transformation evidence", "ECO:0000020 protein inhibition evidence", "ECO:0006404 experimentally evolved mutant phenotypic evidence", "ECO:0000054 double mutant phenotype evidence"]

    # can be more than one of those values in this column, so need to split on the , separating them
    invalid_indices = []
    invalid_codes = []
    
    for index, value in enumerate(evidence_code_list):
        value = value.strip()
        if value == '' or value in ['NA', '-']:
            invalid_indices.append(index)
            continue
        codes = [code.strip() for code in value.split(',')]
        for code in codes:
            if code not in allowable_values:
                invalid_indices.append(index)
                # only add the code to the list if it's got an eco code prefix
                if code.startswith("ECO:"):
                    invalid_codes.append(code)

    if not invalid_indices:
        print("✅ All evidence codes are valid")
        return True
    # otherwise we want to check if the ECO code still starts with "ECO:" but is a code that we haven't got in our suggested list
    # users can have other ECO codes but we want to flag these for potential inclusion in the allowable values
    else:
        if invalid_codes:
            unique_codes = set(invalid_codes)
            print("The following evidence codes are new and not currently in the list of suggested values:")
            print(f"{', '.join(unique_codes)}")
            print("If these are valid ECO codes, please ignore this message from the check.")
        print(f"\n❌ {len(invalid_indices)} rows have failed the check. Each rule must have an evidence code and not be empty. If there are multiple evidence codes for a row, they must be separated by a ',', not by a new line. Evidence codes must start with 'ECO:'.")
        for index in invalid_indices:
            print(f"Row {index + 2}: {evidence_code_list[index]}")
        return False


def check_evidence_grade_limitations(evidence_grade_list, evidence_limitations_list):

    print("\nChecking evidence grade and limitations columns...")

    allowable_grades = ["strong", "moderate", "weak"]
    allowable_limitations = ["lacks evidence for this species", "lacks evidence for this genus", "lacks evidence for this allele", "lacks evidence of the degree to which MIC is affected", "low clinical relevance", "unknown clinical relevance", "statistical geno/pheno evidence but no experimental evidence"]
    
    invalid_indices_grades = []
    invalid_indices_limitations = []

    for index, (grade, limitations) in enumerate(zip(evidence_grade_list, evidence_limitations_list)):
        grade = grade.strip()
        limitations = limitations.strip()
        if grade not in allowable_grades:
            invalid_indices_grades.append(index)
            #print(f"Invalid evidence grade at row {index + 1}: {grade}")
        
        if grade in ["moderate", "weak"]:
            if limitations == '' or limitations == '-':
                invalid_indices_limitations.append(index)
            else:
                # Split limitations by comma and check each one
                limitation_values = [lim.strip() for lim in limitations.split(',')]
                if not all(lim in allowable_limitations for lim in limitation_values):
                    invalid_indices_limitations.append(index)

    if not invalid_indices_grades:
        print("✅ All evidence grades are valid")
    if not invalid_indices_limitations:
        print("✅ All evidence limitations are valid")
    
    if invalid_indices_grades or invalid_indices_limitations:
        print(f"❌ {len(invalid_indices_grades) + len(invalid_indices_limitations)} rows have failed the check")
        if invalid_indices_grades:
            print("Evidence grade column must contain one of the following values:\n" + ", ".join(allowable_grades))
            for index in invalid_indices_grades:
                print(f"Row {index + 2}: {evidence_grade_list[index]}")
        if invalid_indices_limitations:
            print("If evidence grade is 'moderate' or 'weak', evidence limitations column must contain one of the following values: " + ", ".join(allowable_limitations))
            for index in invalid_indices_limitations:
                print(f"Row {index + 2}: {evidence_limitations_list[index]}")
    
    if not invalid_indices_grades and not invalid_indices_limitations:
        return True
    else:
        return False
