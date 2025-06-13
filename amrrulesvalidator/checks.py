"""Validation check functions extracted from legacy script."""

import csv
import re
from pathlib import Path


def check_if_allowed_value(value_list, col_name, allowable_values, missing_allowed=False):

    print("\nChecking  " + col_name + " column...")

    if missing_allowed:
        # can be empty, but this must be reflected by a '-', otherwise must be an allowed value (which should include '-')
        invalid_indices = [index for index, value in enumerate(value_list) if value.strip() == '' or value.strip() in ['NA'] or value.strip() not in allowable_values]
    else:
        # disallow '-' as it must be one of the approved values
        invalid_indices = [index for index, value in enumerate(value_list) if value.strip() == '' or value.strip() in ['NA', '-'] or value.strip() not in allowable_values]

    if not invalid_indices:
        print("✅ All " + col_name + " values are valid")
        return True
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print(col_name + " column must contain one of the following values:\n" + ", ".join(allowable_values))
        for index in invalid_indices:
            print(f"Row {index + 2}: {value_list[index]}")
        return False


def check_if_not_missing(value_list, col_name, list_unique=False):

    print("\nChecking  " + col_name + " column...")

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


def check_ruleIDs(id_list):
    
    print("\nChecking ruleID column...")
    
    invalid_indices = []
    
    if len(id_list) != len(set(id_list)):
        print(f"❌ Rule IDs are not unique")
        invalid_indices = [index for index, value in enumerate(id_list) if id_list.count(value) > 1]
    else:
        print("All rule IDs have passed auto validation")

    prefix = id_list[0][:3]
    prefix_options = [prefix]
    for rule_id in id_list:
        if rule_id[:3] != prefix:
            invalid_indices.append(id_list.index(rule_id))
            # add the different prefix to the list
            prefix_options.append(rule_id[:3])
    
    if not invalid_indices and len(prefix_options) == 1:
        print("✅ All values are valid")
        print(f"Rule prefix: {prefix}")
    else:
        print(f"{len(invalid_indices)} rows have failed the check. Rule IDs must be unique and have the same prefix.")
        for index in invalid_indices:
            print(f"Row {index + 2}: {id_list[index]}")
        print(f"Multiple rule ID prefixes: " + ", ".join(prefix_options))
    
    if not invalid_indices:
        return True
    else:
        return False


def check_organism(txid_list, organism_list, resource_manager=None):
    
    print("\nChecking txid and organism columns...")

    # read in the valid NCBI organism names and their corresponding txids
    if resource_manager:
        # Use the ResourceManager to access the file
        taxonomy_file_path = resource_manager.dir / "ncbi_taxonomy.tsv"
        if not taxonomy_file_path.exists():
            print("❌ Cannot find NCBI taxonomy file. Run 'amrrules update-resources' to download it.")
            return False
    
    # Read the taxonomy file
    with open(taxonomy_file_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        ncbi_organism_dict = {row['Accession']: row['Name'] for row in reader}
    
    # Initialize a dictionary to store invalid rows and reasons
    invalid_txid_indices = {}
    invalid_org_indices = {}
    invalid_org_names = []
    invalid_txids = []

    # for each txid, check that its a value (so not empty, NA, or '-'), and that it's in the ncbi_organism_dict keys
    for index, txid in enumerate(txid_list):
        txid = txid.strip()
        if txid in ['NA', '-', '']:
            invalid_txid_indices[index] = "txid is empty, 'NA', or '-'"
        elif txid not in ncbi_organism_dict.keys():
            invalid_txid_indices[index] = "Taxonomic ID " + txid + " is not in the NCBI taxonomy list"
            invalid_txids.append(txid)
    
    for index, organism in enumerate(organism_list):
        organism = organism.strip()
        # check that the organism name is not empty, NA, or '-'
        if organism in ['NA', '-', '']:
            invalid_org_indices[index] = "Organism name is empty, 'NA', or '-'"
        # check that the organism name starts with 's__' and is in the NCBI organism names list
        elif not organism.startswith('s__'):
            invalid_org_indices[index] = "Organism name " + organism + " does not start with 's__'"
            invalid_org_names.append(organism)
        elif organism.replace('s__', '', 1) not in ncbi_organism_dict.values():
            invalid_org_indices[index] = "Organism name " + organism + " is not in the NCBI taxonomy list"
            invalid_org_names.append(organism)
    
    # for all valid taxids, check that the associated organism name matches the one in the list
    for index, txid in enumerate(txid_list):
        txid = txid.strip()
        if txid in ncbi_organism_dict:
            expected_organism = ncbi_organism_dict[txid]
            # we need to split the 's__' prefix from the organism name before checking
            current_organism = organism_list[index].strip().replace('s__', '', 1)
            if index not in invalid_org_indices and current_organism != expected_organism:
                invalid_org_indices[index] = f"Organism name {current_organism} does not match expected name {expected_organism} for taxid {txid}"
                invalid_org_names.append(current_organism)
    
    if not invalid_txid_indices and not invalid_org_indices:
        print("✅ All txid and organism names passed validation")
    else:
        print(f"❌ {len(invalid_txid_indices) + len(invalid_org_indices)} rows have failed the check")
        print("txids and organism names must be present, not 'NA' or '-'. Organism names should start with 's__'. Both txids and organism names should be in the NCBI taxonomic list, as per file card_ontology/ncbi_taxonomy.tsv.")
        for index in invalid_txid_indices:
            print(f"Row {index + 2}: {txid_list[index]}")
        for index in invalid_org_indices:
            print(f"Row {index + 2}: {organism_list[index]}")
    
    unique_organisms = set(organism_list)
    # if we have some invalid organism names that aren't because the value is empty, and instead
    # becase the value isn't in the GTDB list, go through the GTDB list and extract anything
    # where the genus is the same, and provide those as options for the user to consider
    if len(invalid_org_names) > 0:
        print("\nThe following organism names are not in the NCBI list:\n")
        unique_invalid_org_names = set(invalid_org_names)
        for org_name in unique_invalid_org_names:
            print(org_name)

    unique_organisms_str = ', '.join(map(str, unique_organisms))
    print(f"\nUnique organism names: {unique_organisms_str}")

    if not invalid_txid_indices and not invalid_org_indices:
        return True
    else:
        return False


def check_gene(gene_list, rule_list):
    
    print("\nChecking gene column...")
    
    # want to check if any values are missing, NA, or '-', and return which index number in the list where that's the case
    invalid_indices = [index for index, value in enumerate(gene_list) if value == '' or value in ['NA', '-']]

    if not invalid_indices:
        print("✅ All gene values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Gene names must be present, not 'NA' or '-'")
        for index in invalid_indices:
            print(f"Row {index + 2}: {gene_list[index]}")

    # now we want to check for gene names that are actually combo rules - if there are any, we want to check that any rule IDs mentioned here are present in the file already
    # if there is a value in gene list that follows the format of three capital letters followed by a string of four numbers, this is one to compare against rule ids
    
    print("\nNow checking for combinatorial rules in gene column...")

    invalid_indices_dict = {}

    pattern = re.compile(r'[A-Z]{3}\d{4}')
    for index, gene in enumerate(gene_list):
        if index not in invalid_indices:
            matches = pattern.findall(gene)
            for match in matches:
                if match not in rule_list:
                    invalid_indices_dict[index + 1] = match
                    break
    
    if not invalid_indices_dict:
        print("✅ All gene combinatorial rule IDs are valid")
    
    else:
        print(f"❌ {len(invalid_indices_dict)} rows have failed the check")
        for index, rule_id in invalid_indices_dict.items():
            print(f"Row {index}: ruleID {rule_id} is not present in the list of rules")
    
    if not invalid_indices and not invalid_indices_dict:
        return True
    else:
        return False


def check_id_accessions(nodeID_list, protein_list, nucleotide_list, hmm_list, variation_type_list, refseq_file, node_file, hmm_file):
    
    print("\nChecking nodeID, refseq accession, GenBank accession and HMM accession columns...")

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

    # Check that in combination, at least one of these columns has a value
    invalid_indices = []
    for index, values in enumerate(zip(nodeID_list, protein_list, nucleotide_list, hmm_list)):
        if all(value == '' or value in ['NA', '-'] for value in values):
            # if all the values are empty, check if variation type is 'Combination' for this row
            # if variation type is 'Combination', then this is a valid value
            # only applies if we have the variation type column
            if variation_type_list is not None:
                if variation_type_list[index].strip() == 'Combination':
                    continue
                else:
                    invalid_indices.append(index)
            else:
                invalid_indices.append(index)
    
    # now we need to assess each list individually against the relevant accession lists
    # fine in the value is '-' as this is allowed in individual columns, just not all in combo
    # store the output so we can easily print what rows and columns are the issues
    invalid_node = []
    for index, value in enumerate(nodeID_list):
        value = value.strip()
        if index not in invalid_indices and value not in refseq_node_ids and value != '-':
            invalid_node.append(index)
    invalid_prot = []
    for index, value in enumerate(protein_list):
        value = value.strip()
        if index not in invalid_indices and value not in protein_accessions and value != '-':
            invalid_prot.append(index)
    invalid_nucl = []
    for index, value in enumerate(nucleotide_list):
        value = value.strip()
        if index not in invalid_indices and value not in nucleotide_accessions and value != '-':
            invalid_nucl.append(index)
    invalid_hmm = []
    for index, value in enumerate(hmm_list):
        value = value.strip()
        if index not in invalid_indices and value not in hmm_accessions and value != '-':
            invalid_hmm.append(index)
    
    if not invalid_indices:
        print("✅ All rows contain at least one value in one of these columns.")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check because at least one of either nodeID, refseq accession, GenBank accession and HMM accession must contain a value. The only exception to this is if the variation type is 'Combination', in which case all of these columns can be '-'.")
        for index in invalid_indices:
            print(f"Row {index + 2}")
    if invalid_node or invalid_prot or invalid_nucl or invalid_hmm:
        print(f"❌ One or more accessions aren't present in either the NCBI Reference Gene Catalog (for nodeID, refseq accession and genbank accession) or the NCBI Reference HMM Catalog (for HMM accession). Empty cells must be specified by '-'.")
        print("\nInvalid nodeID accessions values:")
        for index in invalid_node:
            print(f"Row {index + 2}: {nodeID_list[index]}")
        print("\nInvalid protein accession values:")
        for index in invalid_prot:
            print(f"Row {index + 2}: {protein_list[index]}")
        print("\nInvalid nucleotide accession values:")
        for index in invalid_nucl:
            print(f"Row {index + 2}: {nucleotide_list[index]}")
        print("\nInvalid HMM accession values:")
        for index in invalid_hmm:
            print(f"Row {index + 2}: {hmm_list[index]}")
    else:
        print("✅ All accessions are present in the relevant catalogues.")

    if not invalid_indices and not invalid_node and not invalid_prot and not invalid_nucl and not invalid_hmm:
        return True
    else:
        return False


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
            print("❌ Cannot find CARD drug names file. Run 'amrrules update-resources' to download it.")
            return [], []
    else:
        # Fallback to direct file access if no ResourceManager provided
        drug_names_file_path = Path('amrrulesvalidator/resources/card_drug_names.tsv')
        if not drug_names_file_path.exists():
            drug_names_file_path = Path('card_drug_names.tsv')
            if not drug_names_file_path.exists():
                print("❌ Cannot find CARD drug names file. Run 'amrrules update-resources' to download it.")
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

    if not invalid_indices_dict:
        print("✅ All clinical category and breakpoint values are concordant")
        return True
    else:
        print(f"❌ {len(invalid_indices_dict)} rows have failed the check")
        for index, reason in invalid_indices_dict.items():
            print(f"Row {index}: {reason}")
        return False


def check_bp_standard(breakpoint_standard_list):
    
    print("\nChecking breakpoint standard column...")

    # allowable values are ECOFF (Month Year), Name version (year), EUCAST Expected resistant phenotypes [version] ([date]), EUCAST [organism] Expert Rules [version (year])]
    # eg EUCAST v14.0 (2024), ECOFF (May 2025), EUCAST Expected Resistant Phenotypes v1.2 (2023)
    # we need regex to check these

    suggested_values = [
        r'^ECOFF \(\w+ \d{4}\)$',  # ECOFF (Month Year)
        r'^EUCAST .+ Expert Rules \(\w+ \d{4}\)$',  # EUCAST [organism] Expert Rules (Month year)
        r'^EUCAST Expected Resistant Phenotypes v\d+(\.\d+)? \(\d{4}\)$',  # EUCAST Expected Resistant Phenotypes version (year)
        r'^(EUCAST|CLSI)\s+v\d+(\.\d+)?\s+\(\d{4}\)$'  # EUCAST/CLSI version (year)
    ]
    invalid_indices_dict = {}
    unique_values = set(breakpoint_standard_list)
    for index, value in enumerate(breakpoint_standard_list):
        value = value.strip()
        if value == '' or value in ['NA', '-']:
            continue
        if not any(re.match(pattern, value) for pattern in suggested_values):
            invalid_indices_dict[index + 2] = value
    if not invalid_indices_dict:
        print("✅ All breakpoint standard values match expected patterns.")
        print("Here is a list of the unique values found in the column:")
        unique_values_str = '\n'.join(unique_values)
        print(unique_values_str)
        return True
    else:
        print(f"❌ {len(invalid_indices_dict)} rows have failed the check")
        print("We check for the following formats: ECOFF (Month Year), EUCAST [organism] Expert Rules (Month year), EUCAST Expected Resistant Phenotypes vX (year), or EUCAST/CLSI vX (year).")
        print("Breakpoint standard values that didn't match this format are:")
        for index, value in invalid_indices_dict.items():
            print(f"Row {index}: {value}")
        print("Please double check these entries to ensure they are valid.")

        print("Here is a list of the unique values found in the column:")
        unique_values_str = '\n'.join(unique_values)
        print(unique_values_str)

        return False


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
