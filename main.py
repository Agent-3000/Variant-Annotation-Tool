import streamlit as st
import pandas as pd
from src import annotation, visualization, data_processing, database



def main():
    # Load custom CSS
    with open("src/style.css", "r") as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

    st.title("Variant Annotation Tool")

    # Upload VCF file
    st.sidebar.header("Upload VCF File")
    vcf_file = st.sidebar.file_uploader("Upload VCF file", type=["vcf"])

    if vcf_file is not None:
        # Load VCF data
        vcf_df = pd.read_csv(vcf_file, sep="\t", header=None, skiprows=1)
        vcf_df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]

        # Perform variant annotation
        annotated_df = annotation.annotate_variants(vcf_df, ['allele_frequencies', 'pathogenicity_scores', 'clinical_annotations', 'gene_protein_info'])

        # Display annotated variants
        st.subheader("Annotated Variants")
        st.write(annotated_df)

        # Filter and sort variants
        st.sidebar.header("Filter Variants")
        filters = st.sidebar.text_input("Enter filters (e.g., af_gnomad>=0.01,pathogenicity_score<0.5)")
        filtered_df = data_processing.filter_variants(annotated_df, filters)

        st.sidebar.header("Sort Variants")
        sort_by = st.sidebar.text_input("Enter column(s) to sort by (e.g., af_gnomad, pathogenicity_score)")
        ascending = st.sidebar.checkbox("Sort in ascending order", value=True)
        sorted_df = data_processing.sort_variants(filtered_df, sort_by.split(","), ascending)

        st.subheader("Filtered and Sorted Variants")
        st.write(sorted_df)

        # Visualize annotated variants
        st.subheader("Visualizations")
        visualization_type = st.selectbox("Select visualization type", ["Histogram", "Scatter Plot"])
        if visualization_type == "Histogram":
            column = st.selectbox("Select column for histogram", sorted_df.columns)
            visualization.plot_variant_histogram(sorted_df, column)
        elif visualization_type == "Scatter Plot":
            x_column = st.selectbox("Select column for x-axis", sorted_df.columns)
            y_column = st.selectbox("Select column for y-axis", sorted_df.columns)
            hue_column = st.selectbox("Select column for color (optional)", sorted_df.columns + [None], index=len(sorted_df.columns))
            visualization.plot_variant_scatter(sorted_df, x_column, y_column, hue_column)

        # Export annotated variants
        st.sidebar.header("Export Data")
        file_format = st.sidebar.selectbox("Select file format", ["csv", "tsv", "xlsx"])
        if st.sidebar.button("Export"):
            output_file = st.sidebar.text_input("Enter output file path")
            data_processing.export_variants(sorted_df, file_format, output_file)
            st.success(f"Data exported to {output_file}")

if __name__ == "__main__":
    main()