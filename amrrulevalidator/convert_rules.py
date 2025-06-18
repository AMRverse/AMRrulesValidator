"""Aids with conversion of rules to the latest version of the AMRrules spec"""

from pathlib import Path
import csv

from amrrulevalidator.constants import CANONICAL_COLUMNS, SPEC_VERSION
from amrrulevalidator.utils.io import read_tsv, write_tsv
from amrrulevalidator.validate import get_column
from amrrulevalidator.utils.resources import ResourceManager

def run_convert_to_latest_spec(input_p: Path, output_p: Path, rm: ResourceManager) -> bool:
    """
    Convert an AMRrules file to the latest specification.
    
    Args:
        input_p: Path to the input TSV file
        output_p: Path to the output TSV file
        
    Returns:
        bool: True if conversion was successful, False otherwise
    """
    
    print(f"\nConverting rules file to latest spec: {input_p}")
    
    # Read rows from input file
    rows = read_tsv(input_p)
    
    # Currently converts spec v0.5 to v0.6
    # first add the txid column
    print(f"\nAdding txid column")
    # get the organism list
    organism_list = get_column("organism", rows)
    # get the ncbi taxonomy dict
    taxonomy_file_path = rm.dir / "ncbi_taxonomy.tsv"
    if not taxonomy_file_path.exists():
        print("❌ Cannot find NCBI taxonomy file. Run 'amrrule update-resources' to download it.")
        ncbi_organism_dict = None

    # Read the taxonomy file and store all the organisms and their txids in a dictionary
    with open(taxonomy_file_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        ncbi_organism_dict = {row['Name']: row['Accession'] for row in reader}
    
   # now add the txid column to each row, make UNKNOWN if we can't find it
    unknown_orgs = []
    for index, organism in enumerate(organism_list):
        organism = organism.strip().replace('s__', '')
        if organism in ncbi_organism_dict:
            rows[index]['txid'] = ncbi_organism_dict[organism]
        else:
            unknown_orgs.append(organism)
            rows[index]['txid'] = "UNKNOWN"

    if unknown_orgs:
        print(f"❌ Some organisms not found in NCBI taxonomy. Setting txid to UNKNOWN for these organisms: {', '.join(set(unknown_orgs))}")
    
    # now we need to fix the accession colummns. First, refseq and genbank accession need to be split into nucelotide and protein accessions.
    print("Adding protein and nucleotide accession columns")
    # by default, if we can't find out if the accession is nucleotide or protein, we're going to put it in the protein accession column but with the note to CHECK ACCESSION TYPE
    refseq_accessions = get_column("refseq accession", rows)
    genbank_accessions = get_column("GenBank accession", rows)

    # get the refgene prot and nucl accessions to check against
    refseq_prot_accessions, refseq_nucl_accessions = rm.refseq_accessions()
    for index, (ref_acc, gen_acc) in enumerate(zip(refseq_accessions, genbank_accessions)):
        prot_found = False
        nucl_found = False
        ref_assigned = False
        gen_assigned = False
        if ref_acc == "-" and gen_acc == "-":
            # if both accessions are missing, set both to empty string
            rows[index]['protein accession'] = "-"
            rows[index]['nucleotide accession'] = "-"
            continue
        # if we have a non-empty value for both accessions, we
        # check refseq accessions
        if ref_acc in refseq_prot_accessions:
            rows[index]['protein accession'] = ref_acc
            prot_found = True
            ref_assigned = True
        elif ref_acc in refseq_nucl_accessions:
            rows[index]['nucleotide accession'] = ref_acc
            nucl_found = True
            ref_assigned = True

        # check genbank accessions
        if gen_acc in refseq_prot_accessions:
            rows[index]['protein accession'] = gen_acc
            prot_found = True
            gen_assiged = True
        elif gen_acc in refseq_nucl_accessions:
            rows[index]['nucleotide accession'] = gen_acc
            nucl_found = True
            gen_assigned = True

        # if we didn't find a protein accession, we still need to add a value here
        if not prot_found:
            # let's build our possible string
            prot_string = []
            if not ref_assigned and ref_acc != "-":
                prot_string.append(ref_acc)
            if not gen_assigned and gen_acc != "-":
                prot_string.append(gen_acc)
            if prot_string:
                rows[index]['protein accession'] = "CHECK ACCESSION TYPE: " + ", ".join(prot_string)
        # do the same for nucleotide accession
        if not nucl_found:
            nucl_string = []
            if not ref_assigned and ref_acc != "-":
                nucl_string.append(ref_acc)
            if not gen_assigned and gen_acc != "-":
                nucl_string.append(gen_acc)
            if nucl_string:
                rows[index]['nucleotide accession'] = "CHECK ACCESSION TYPE: " + ", ".join(nucl_string)
    
    # now we can remove the refseq and genbank accession columns
    print("Removing refseq and genbank accession columns")
    for row in rows:
        del row['refseq accession']
        del row['GenBank accession']
    
    # change the name of the context column to gene context
    print("Renaming 'context' column to 'gene context'")
    gene_context_col = get_column("context", rows)
    for index, context in enumerate(gene_context_col):
        rows[index]['gene context'] = context
        del rows[index]['context']

    print("Adding breakpoint condition column")
    for row in rows:
        row['breakpoint condition'] = "ADD VALUE"

    # update evidece grade to be 'high' instead of strong, and for 'weak' made a note that this could be low or very low
    print("Updating evidence grade values")
    evidence_grade_col = get_column("evidence grade", rows)
    for index, grade in enumerate(evidence_grade_col):
        if grade.lower() == "strong":
            rows[index]['evidence grade'] = "high"
        elif grade.lower() == "weak":
            rows[index]['evidence grade'] = "UPDATE TO low or very low"
    
    # Write converted rows to output file
    write_tsv(rows, output_p, columns=CANONICAL_COLUMNS)
    
    print(f"Conversion complete. Output written to {output_p}")
    return True