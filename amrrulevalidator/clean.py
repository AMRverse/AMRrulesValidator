"""AMR rules cleaning module."""

import re
from pathlib import Path

from amrrulevalidator.utils.io import read_tsv, write_tsv

def run_clean(input_p: Path, output_p: Path) -> bool:
    """
    Clean an AMRrules file by checking for validation markers and removing
    trailling whitespaces or blank lines.
    
    Args:
        input_p: Path to the input TSV file (should be validated first)
        output_p: Path to the output TSV file (clean version)
        
    Returns:
        bool: True if cleaning was successful, False otherwise
    """
    print(f"\nCleaning rules file: {input_p}")

    # Read rows from input file
    rows = read_tsv(input_p)
    
    if not rows:
        print("❌ Input file is empty or contains no valid data.")
        return False
    
    # Check for validation markers
    validation_issues = []
    for i, row in enumerate(rows):
        for col, value in row.items():
            # Look for CHECK or ENTRY MISSING markers in any cell
            if re.search(r'(CHECK|ENTRY MISSING)', str(value)):
                validation_issues.append({
                    'row': i + 1,  # 1-based row index for user-friendly reporting
                    'column': col,
                    'value': value
                })
    
    # If validation issues exist, report them and exit without writing output
    if validation_issues:
        print("\n❌ Validation markers still prsent in the file. Please resolve and remove these before cleaning:")
        for issue in validation_issues:
            print(f"  Row {issue['row']}, Column '{issue['column']}': {issue['value']}")
        print("\nNo output file was written.")
        return False
    
    # Clean the data:
    # 1. Strip trailing whitespace from all cells
    # 2. Keep all columns in their original order
    # 3. Remove rows where all values are empty or only whitespace
    columns = list(rows[0].keys()) if rows else []
    
    for row in rows:
        for col in row:
            if isinstance(row[col], str):
                row[col] = row[col].rstrip()
    
    # Filter out rows where all values are empty or whitespace
    rows = [row for row in rows if any(str(value).strip() for value in row.values())]
    
    # Write cleaned data to output file
    write_tsv(rows, output_p, columns)
    print(f"✅ Cleaned file written to {output_p}")
    return True
