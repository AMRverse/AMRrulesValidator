"""AMR rules validator module."""

from pathlib import Path

from amrrulesvalidator.constants import CANONICAL_COLUMNS
from amrrulesvalidator.utils.io import read_tsv, write_tsv
from amrrulesvalidator.utils.resources import ResourceManager
from amrrulesvalidator.checks import *


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
                row[col] = "EMPTY: MISSING VALUE"
    
    # Write the processed rows to the output file
    write_tsv(rows, output_p, CANONICAL_COLUMNS)
    
    return True
