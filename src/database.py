import sqlite3
from typing import List, Dict, Any

# Database connection
conn = sqlite3.connect('annotations.db')
cursor = conn.cursor()

def create_tables():
    """
    Create the necessary tables in the database if they don't exist.
    """
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS variants (
            id INTEGER PRIMARY KEY,
            chrom TEXT,
            pos INTEGER,
            ref TEXT,
            alt TEXT,
            gene TEXT,
            transcript TEXT,
            annotation TEXT,
            af_gnomad REAL,
            af_1000g REAL,
            pathogenicity_score REAL,
            clinvar_significance TEXT
        )
    """)
    conn.commit()

def insert_variant(variant: Dict[str, Any]):
    """
    Insert a new variant into the database.
    
    Args:
        variant (Dict[str, Any]): A dictionary containing the variant information.
    """
    columns = ', '.join(variant.keys())
    placeholders = ', '.join('?' * len(variant))
    query = f'INSERT INTO variants ({columns}) VALUES ({placeholders})'
    cursor.execute(query, tuple(variant.values()))
    conn.commit()

def get_annotations(chrom: str, pos: int) -> List[Dict[str, Any]]:
    """
    Retrieve annotations for a specific variant from the database.
    
    Args:
        chrom (str): Chromosome for the variant.
        pos (int): Position for the variant.
    
    Returns:
        List[Dict[str, Any]]: A list of dictionaries containing the variant annotations.
    """
    query = 'SELECT * FROM variants WHERE chrom = ? AND pos = ?'
    cursor.execute(query, (chrom, pos))
    rows = cursor.fetchall()
    annotations = [dict(zip([col[0] for col in cursor.description], row)) for row in rows]
    return annotations

# Create the tables if they don't exist
create_tables()