import pandas as pd
import allel
import requests
import xml.etree.ElementTree as ET
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio import SeqIO


# Allele Frequency Retrieval
def get_allele_frequencies(variants, database):
    """
    Retrieve allele frequencies for the given variants from the specified database.
    """
    if database == 'gnomad':
        # Query the gnomAD API for allele frequencies
        for variant in variants:
            chrom, pos, ref, alt = variant.CHROM, variant.POS, variant.REF, variant.ALT
            url = f"https://gnomad.broadinstitute.org/variant/{chrom}-{pos}-{ref}-{alt}?dataset=gnomad_r4_1"
            response = requests.get(url)
            if response.ok:
                data = response.json()
                af = data.get('gnomad_exome_af', 0.0)
                variant.set_annotations({'gnomad_af': af})

    return variants

# Pathogenicity Score Calculation
def calculate_pathogenicity_scores(variants):
    """
    Calculate pathogenicity scores (e.g., SIFT, PolyPhen) for the given variants.
    """
    for variant in variants:
        chrom, pos, ref, alt = variant.CHROM, variant.POS, variant.REF, variant.ALT
        seq = get_sequence(chrom, pos - 100, pos + 100)  # Retrieve sequence around the variant
        sift_score = calculate_sift_score(seq, pos, ref, alt)
        polyphen_score = calculate_polyphen_score(seq, pos, ref, alt)
        variant.set_annotations({'sift_score': sift_score, 'polyphen_score': polyphen_score})

    return variants

# Clinical Annotation Retrieval
def get_clinical_annotations(variants, database):
    """
    Retrieve clinical annotations and disease associations for the given variants from the specified database.
    """
    if database == 'clinvar':
        # Query the ClinVar API for clinical annotations
        for variant in variants:
            chrom, pos, ref, alt = variant.CHROM, variant.POS, variant.REF, variant.ALT
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&id={chrom}:{pos}:{ref}:{alt}"
            response = requests.get(url)
            if response.ok:
                root = ET.fromstring(response.content)
                clinical_significance = root.find('ClinicalSignificance/Description').text
                variant.set_annotations({'clinvar_significance': clinical_significance})

    return variants

# Gene and Protein Information Lookup
def get_gene_protein_info(variants, database):
    """
    Retrieve gene and protein information (e.g., gene names, transcript IDs, protein domains) for the given variants from the specified database.
    """
    if database == 'ensembl':
        # Query the Ensembl API for gene and protein information
        for variant in variants:
            chrom, pos = variant.CHROM, variant.POS
            url = f"https://rest.ensembl.org/overlap/region/{chrom}:{pos}-{pos}?feature=gene"
            response = requests.get(url, headers={"Content-Type": "application/json"})
            if response.ok:
                data = response.json()
                if data:
                    gene_name = data[0]['external_name']
                    transcript_id = data[0]['transcript_id']
                    variant.set_annotations({'gene_name': gene_name, 'transcript_id': transcript_id})

    return variants

def calculate_sift_score(seq, pos, ref, alt):
    # Convert DNA sequences to protein sequences
    protein_seq = Seq(str(seq)).translate()
    protein_ref = str(protein_seq[pos - 1])
    protein_alt = str(Seq(str(alt)).translate())

    # Calculate SIFT score using BLOSUM62 matrix
    blosum62 = substitution_matrices.load("BLOSUM62")
    sift_score = blosum62[(protein_ref, protein_alt)]

    return sift_score

def calculate_polyphen_score(seq, pos, ref, alt):
    # Prepare data for PolyPhen server
    protein_seq = Seq(str(seq)).translate()
    protein_ref = str(protein_seq[pos - 1])
    protein_alt = str(Seq(str(alt)).translate())
    data = {
        'sequence': str(protein_seq),
        'position': pos,
        'refaa': protein_ref,
        'aaalt': protein_alt
    }

    # Send request to PolyPhen server
    response = requests.get('http://genetics.bwh.harvard.edu/pph2/bgi.shtml', params=data)

    # Extract PolyPhen score from response
    polyphen_score = re.search(r'PolyPhen Score: ([\d.]+)', response.text)
    if polyphen_score:
        return float(polyphen_score.group(1))
    else:
        return None

# Helper function to retrieve sequence
def get_sequence(chrom, start, end):
    # Load the reference genome sequence from a FASTA file
    with open('reference_genome.fa', 'r') as fasta_file:
        sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

    # Extract the sequence for the given chromosome and coordinates
    ref_seq = sequences[chrom].seq
    sequence = str(ref_seq[start - 1:end])

    return sequence