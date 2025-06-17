"""Helper functions for validation checks."""
import re


def report_check_results(check_name, invalid_indices=None, invalid_dict=None, 
                        success_message="All values are valid", 
                        failure_message=None, 
                        unique_values=None):
    """
    Unified reporting function for check results.
    
    Args:
        check_name: Name of the check being performed
        invalid_indices: List of invalid indices (optional)
        invalid_dict: Dictionary of invalid indices with reasons (optional)
        success_message: Message to display on success
        failure_message: Message to display on failure
        unique_values: Set of unique values to display (optional)
    
    Returns:
        bool: True if check passed, False otherwise
    """
    # Convert None to empty collections
    invalid_indices = invalid_indices or []
    invalid_dict = invalid_dict or {}
    
    total_invalid = len(invalid_indices) + len(invalid_dict)
    
    if total_invalid == 0:
        print(f"✅ {success_message}")
        if unique_values:
            unique_values_str = ', '.join(map(str, unique_values))
            print(f"\nUnique {check_name} values: {unique_values_str}")
        return True
    else:
        print(f"❌ {total_invalid} rows have failed the check")
        if failure_message:
            print(failure_message)
        
        for index in invalid_indices:
            print(f"Row {index + 2}")
            
        for index, reason in invalid_dict.items():
            print(f"Row {index + 2}: {reason}")
                
        if unique_values:
            unique_values_str = ', '.join(map(str, unique_values))
            print(f"\nUnique {check_name} values: {unique_values_str}")
            
        return False


def validate_pattern(value_list, patterns, rows, col_name, missing_allowed=False):

    invalid_indices = {}
    for index, value in enumerate(value_list):
        value = value.strip()
        if not missing_allowed and (value == '' or value == 'NA'):
            rows[index][col_name] = 'ENTRY MISSING'
            invalid_indices[index] = f"Value in {col_name}"
        if not any(re.match(pattern, value) for pattern in patterns) or value == '-':
            invalid_indices[index] = f"Value '{value}' in {col_name} does not match the expected patterns"
            rows[index][col_name] = f'CHECK VALUE: {value}'
    
    return invalid_indices, rows


def check_values_in_list(value_list, allowed_values, col_name, rows, missing_allowed=False, fail_reason=None):
    """
    Check if all values in a list are in the set of allowed values.
    
    Args:
        value_list: List of values to check
        allowed_values: List of allowed values
        col_name: Name of the column being checked
        missing_allowed: Whether missing values are allowed
        
    Returns:
        dict: Dictionary of invalid indices with reasons
        rows: List of rows with invalid values flagged
    """
    invalid_dict = {}

    if missing_allowed:
        allowed_values = set(allowed_values) | {'-'}

    for index, value in enumerate(value_list):
        value = value.strip()
        if not missing_allowed and value in ['NA', '-', '']:
            invalid_dict[index] = f"Missing value in {col_name}"
            rows[index][col_name] = 'ENTRY MISSING'
            continue
        if value not in allowed_values:
            invalid_dict[index] = f"{value} {fail_reason}"
            rows[index][col_name] = 'CHECK VALUE: ' + value
    
    return invalid_dict, rows


def check_if_col_empty(col_values, col_name, rows):
    """
    Check if a column is empty or contains only missing values.
    
    Args:
        col_values: List of values in the column
        col_name: Name of the column
        
    Returns:
        tuple: (bool indicating if the column is empty, list of rows with missing values flagged in that column)
    """
    if all(value.strip() in ['NA', '-', ''] for value in col_values):
        print(f"Column '{col_name}' is empty or contains only missing values.")
        for row in rows:
            row[col_name] = "ENTRY MISSING"
        return True, rows
    else:
        return False, rows 