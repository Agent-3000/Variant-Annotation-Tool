import pandas as pd

def filter_variants(variants, filters):
    """
    Filter the annotated variants based on the specified filters.
    
    Args:
        variants (pandas.DataFrame): DataFrame containing the annotated variants.
        filters (dict): A dictionary containing the filter conditions.
    
    Returns:
        pandas.DataFrame: DataFrame containing the filtered variants.
    """
    filtered_variants = variants.copy()
    
    # Apply filters based on the provided conditions
    for column, condition in filters.items():
        if condition.startswith(('>=', '>', '<=', '<', '==')):
            operator = condition[:2]
            value = float(condition[2:])
            filtered_variants = filtered_variants[eval(f'filtered_variants["{column}"] {operator} {value}')]
        else:
            filtered_variants = filtered_variants[filtered_variants[column].isin(condition.split(','))]
    
    return filtered_variants

def sort_variants(variants, sort_by, ascending=True):
    """
    Sort the annotated variants based on the specified column(s).
    
    Args:
        variants (pandas.DataFrame): DataFrame containing the annotated variants.
        sort_by (str or list): Column(s) to sort the DataFrame by.
        ascending (bool, optional): Whether to sort in ascending or descending order. Default is True (ascending).
    
    Returns:
        pandas.DataFrame: DataFrame containing the sorted variants.
    """
    sorted_variants = variants.sort_values(by=sort_by, ascending=ascending)
    return sorted_variants

def export_variants(variants, file_format, output_file):
    """
    Export the annotated variants to the specified file format.
    
    Args:
        variants (pandas.DataFrame): DataFrame containing the annotated variants.
        file_format (str): File format to export the data (e.g., 'csv', 'tsv', 'xlsx').
        output_file (str): Path to the output file.
    """
    if file_format == 'csv':
        variants.to_csv(output_file, index=False)
    elif file_format == 'tsv':
        variants.to_csv(output_file, sep='\t', index=False)
    elif file_format == 'xlsx':
        variants.to_excel(output_file, index=False)
    else:
        raise ValueError(f'Unsupported file format: {file_format}')