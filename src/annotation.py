import pandas as pd
import requests
from json import JSONDecodeError
import time
import io
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio import SeqIO

def get_allele_frequencies(variants, database):
    if database == 'ensembl':
        server = "https://rest.ensembl.org"
        ext = "/vep/human/hgvs"

        for variant in variants:
            chrom = variant.CHROM
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            pos = variant.POS
            ref = variant.REF
            alt = variant.ALT
            
            hgvs = f"{chrom}:g.{pos}{ref}>{alt}"
            headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
            
            r = requests.post(server+ext, headers=headers, data='{ "hgvs_notations" : ["'+hgvs+'" ]}')
            if not r.ok:
                r.raise_for_status()
                sys.exit()

            decode = r.json()
            allele_freqs = decode[0]["colocated_variants"][0]["frequencies"]
            variant.set_annotations({"allele_freqs": allele_freqs})

    return variants


def calculate_polyphen_score(seq, pos, ref, alt):
    protein_seq = str(Seq(str(seq)).translate())
    protein_ref = protein_seq[pos - 1]
    protein_alt = str(Seq(str(alt)).translate())

    url = f"http://genetics.bwh.harvard.edu/cgi-bin/pph2?query={protein_ref},{protein_alt}&query_aa=1&transcript_id={seq.id}&output_format=json"
    response = requests.get(url)
    if response.ok:
        data = response.json()
        polyphen_score = data['polyphen2_hdiv_score']
        return polyphen_score
    else:
        return None

def annotate_variants(vcf_data, annotation_types):
    vcf_df = pd.read_table(io.StringIO(vcf_data), comment='#', header=None)

    num_cols = len(vcf_df.columns)
    column_map = {
        vcf_df.columns[0]: "CHROM",
        vcf_df.columns[1]: "POS",
        vcf_df.columns[2]: "ID",
        vcf_df.columns[3]: "REF",
        vcf_df.columns[4]: "ALT",
    }

    # Add the rest of the columns to the dictionary if necessary
    if num_cols > 5:
        column_map[vcf_df.columns[5]] = "QUAL"
    if num_cols > 6:
        column_map[vcf_df.columns[6]] = "FILTER"
    if num_cols > 7:
        column_map[vcf_df.columns[7]] = "INFO"
    if num_cols > 8:
        column_map[vcf_df.columns[8]] = "FORMAT"
    if num_cols > 9:
        column_map[vcf_df.columns[9]] = "SAMPLE"

    vcf_df = vcf_df.rename(columns=column_map)

    variants = []
    for _, row in vcf_df.iterrows():
        variant = Variant()
        variant.CHROM = row["CHROM"]
        variant.POS = row["POS"]
        variant.REF = row["REF"]
        variant.ALT = row["ALT"]

        if "allele_frequencies" in annotation_types:
            get_allele_frequencies([variant], 'gnomad')

        variants.append(variant)

    return variants


class Variant:
    def __init__(self):
        self.annotations = {}

    def set_annotations(self, annotations):
        self.annotations.update(annotations)