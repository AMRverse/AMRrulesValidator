"""I/O utilities for reading and writing TSV files."""

import csv
import subprocess
from pathlib import Path
from typing import Any

from ..constants import CANONICAL_COLUMNS


def detect_encoding(file_path):
    """
    Detect the encoding of a file using the `file -I` command.
    Args:
        file_path (str): Path to the file.
    Returns:
        str: Detected encoding or None if the encoding cannot be determined.
    """
    try:
        # Run the `file -I` command
        result = subprocess.run(['file', '-I', file_path], capture_output=True, text=True, check=True)
        output = result.stdout.strip()

        # Extract the charset from the output
        for part in output.split():
            if "charset=" in part:
                return part.split("=")[1]

        print("❌ Unable to determine encoding from `file -I` output.")
        return None
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running `file -I` command: {e}")
        return None
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return None


def read_tsv(path: Path) -> list[dict[str, Any]]:
    """
    Read a TSV file and return list of dictionaries.
    
    Mimics the behavior from the legacy script:
    - Uses detect_encoding to determine file encoding
    - Strips whitespace from column names
    - Filters out completely empty rows
    - Removes the second row if it contains 'required' in the first column
    
    Args:
        path: Path to the TSV file
        
    Returns:
        List of dictionaries where each dict represents a row
    """
    encoding_type = detect_encoding(str(path))
    
    with open(path, newline='', encoding=encoding_type) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        columns = [col.rstrip() for col in reader.fieldnames]
        
        # Filter out rows where all values are empty
        rows = [row for row in reader if any(value.strip() for value in row.values())]
        
        # Remove the second row if the first column contains "required"
        # This is common in AMR rules files to indicate required vs optional columns
        if len(rows) > 1 and rows[0][columns[0]].strip().lower() == "required":
            del rows[0]
    
    return rows


def write_tsv(rows: list[dict[str, Any]], path: Path, columns: list[str] = None) -> None:
    """
    Write rows to a TSV file in canonical column order.
    
    Args:
        rows: List of dictionaries to write
        path: Output path for the TSV file
        columns: List of column names in desired order. If None, uses CANONICAL_COLUMNS
    """
    if not rows:
        # Create empty file if no rows
        path.touch()
        return
    
    if columns is None:
        columns = CANONICAL_COLUMNS
    
    # Use UTF-8 encoding for output
    with open(path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
