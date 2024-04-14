import streamlit as st
import pandas as pd
from src import annotation, visualization

def main():
    st.title("Variant Annotation Tool")

    # Upload VCF file
    st.sidebar.header("Upload VCF File")
    vcf_file = st.sidebar.file_uploader("Upload VCF file", type=["vcf"])

    if vcf_file is not None:
        # Load VCF data
        vcf_df = pd.read_csv(vcf_file, sep="\t", header=None)
        vcf_df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]

        # Perform variant annotation
        annotated_df = annotation.annotate_variants(vcf_df)

        # Display annotated variants
        st.subheader("Annotated Variants")
        st.write(annotated_df)

        # Visualize annotated variants
        visualization.plot_variant_histogram(annotated_df)

if __name__ == "__main__":
    main()