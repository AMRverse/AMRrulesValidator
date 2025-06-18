# AMRrulevalidator

The AMRrulevalidator package provides tools for validating [AMRrules files](https://github.com/interpretAMR/AMRrules) according to the current specification ([v0.6](https://github.com/AMRverse/AMRrulesCuration)).

As part of the install, AMRrulevalidator will download the current CARD ontology and AMRFinderPlus resources to validate against.

The `validate` subcommand will print to stdout a summary of checks that have been completed, and whether they've passed or failed, and will write out a version of the rules file where cells are flagged with values that need to be checked. 

The `clean` subcommand will write out a cleaned version of a rules file after all values have been checked, which will be ready for integration into the AMRrules interpretation engine.

## Installation

AMRrulevalidator is compatible with Python >= 3.8. The only dependency is [obonet v1.1.1](https://pypi.org/project/obonet/).

The easiest installation method is via pip from GitHub, which will install all required dependencies for you.

### Installation with pip

```bash
# Optional: create a conda environment to install package into
conda create amrrulevalidator

conda activate amrrulevalidator

# Install with pip from GitHub
pip install git+https://github.com/AMRverse/AMRrulevalidator.git

# Download CARD and AMRFinderPlus ontology files
amrrule update-resources
```

### Installation for development purposes

For development or local installation. This will:
1. Install the package in editable mode
2. Download and set up all required resource files

```bash
# Clone the repository
git clone https://github.com/AMRverse/AMRrulevalidator.git
cd AMRrulevalidator

# Install in development mode
make dev
```

## Updating Resources

The validator relies on several external files including the CARD ontology and files from the AMRFinderPlus database. To update these files:

```bash
amrrule update-resources
```

This will download the latest versions of:
- CARD ontology (ARO) (v4.0.1)
- CARD drug names and drug classes (v4.0.1)
- NCBI taxonomy data (from CARD (v4.0.1))
- AMRFinderPlus Reference Gene Hierarchy, Reference Gene Accessions and HMM accessions (using `latest` version)

## Usage

### Validating a rules file

To validate a draft AMRrules file:

```bash
amrrule validate --input path/to/draft_rules.tsv --output path/to/validated_rules.tsv
```

This will:
1. Check the input file against the current AMRrules specification
2. Generate a validated output file with annotations for any errors
3. Print a summary of validation results to the console

### Validation checks

During validation, the script annotates problematic values in the output file to help identify and fix issues:

- `ENTRY MISSING`: Indicates that a required value is missing in a field where a value is expected.
- `CHECK VALUE: [value]`: Indicates that the existing value doesn't match the expected format or isn't in the list of allowed values.

The rules files must contain the following columns:

- ruleID
- txid
- organism
- gene
- nodeID
- protein accession
- HMM accession
- nucleotide accession
- ARO accession
- mutation
- variation type
- gene context
- drug
- drug class
- phenotype
- clinical category
- breakpoint
- breakpoint standard
- breakpoint condition
- PMID
- evidence code
- evidence grade
- evidence limitations
- rule curation note

Any columns which do not exist in the file will be added, with all values sent to `ENTRY MISSING`.

#### Validation check logic

The validator performs a series of checks, with each focusing on specific columns:

1. **ruleID**: Checks that rule IDs are unique and have a consistent prefix.

2. **txid**: Validates that taxonomic IDs exist in the NCBI taxonomy database.

3. **organism**: Confirms that organism names are valid NCBI taxonomy names and follow the format `s__[organism name]`.

4. **txid-organism pairs**: Ensures that each txid is correctly paired with its corresponding organism name.

5. **gene**: Checks that the gene column is not empty. For combination rules, it verifies that referenced rule IDs exist.

6. **Accession checks**:
   - nodeID: Verifies node IDs against the AMRFinderPlus Reference Gene Hierarchy.
   - protein accession: Verifies accessions against the AMRFinderPlus Reference Gene Catalog.
   - nucleotide accession: Verifies accessions against the AMRFinderPlus Reference Gene Catalog.
   - HMM accession: Verifies accessions against the AMRFinderPlus HMM accessions list.
   - At least one of these accessions must be present unless the variation type is "Combination".

   _Note:_ Accessions are checked only against the listed reference files - if you have an accession that has come from elsewhere, this value may be flagged as `CHECK VALUE:`.

7. **ARO accession**: Validates ARO accessions against the CARD ontology.

8. **variation type**: Confirms values match one of the allowed variation types, as per the spec. Value must be supplied, cannot be empty.

9. **mutation and variation type compatibility**: Ensures the mutation format is compatible with the specified variation type.

10. **gene context**: Validates against allowed values (`core` or `acquired`). Value must be supplied, cannot be empty.

11. **drug and drug class**: Confirms values exist in the CARD drug and drug class lists. A value in one of these columns must be supplied, cannot be empty.

12. **phenotype**: Validates against allowed values (`wildtype` or `nonwildtype`). Value must be supplied, cannot be empty.

13. **clinical category**: Confirms values match the allowed clinical categories (`S`, `I`, or `R`). Value must be supplied, cannot be empty.

14. **breakpoint**: Checks to see if the breakpoint value is consistent with the clinical category. If not, flag for checking. A value must be supplied, if no breakpoint is required then `not applicable` is a valid entry.

15. **breakpoint standard**: Checks to see if the given breakpoint standard source includes information about version, or month/year when standard was set. Flags as a value to check if it doesn't match.

16. **breakpoint condition**: If provided, confirms values match the allowed breakpoint conditions.

17. **PMID**: Checks only that there is an entry in this column, as most rules should have a paper associated with them.

18. **evidence code**: Checks that the codes provided start with the `ECO:` prefix. Will check if they are listed as one of the suggested Evidence Code Ontology codes in the spec - if not, flags as a code to check manually. Checks to make sure multiple codes are separated by a comma and not some other delimiter.

19. **evidence grade and limitations**: Validates evidence grades against allowed values. For evidence limitations, checks that mulitple limitations are separated by a comma, and not a different delimiter. Checks limitations are one of the allowed values. Checks that if evidence grade is not `high`, an evidence limitation is provided.

Allowable values for some columns are specified within `constants.py`, and should match the current version of the AMRrules spec:

### Cleaning a rules file

After validation, you can clean a rules file to prepare it for import into the interpretation engine:

```bash
amrrule clean --input path/to/validated_rules.tsv --output path/to/clean_rules.tsv
```

## License

This project is licensed under the GNU General Public License v3.0.

## Contributors

Code was developed by Jane Hawkey, with input from Kat Holt and Natacha Couto.

## Contact

For issues or questions, please use the [GitHub issue tracker](https://github.com/AMRverse/AMRrulevalidator/issues).
