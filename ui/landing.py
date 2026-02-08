import streamlit as st
from Bio import SeqIO
from io import StringIO
import time

from bwt_builder.fmi import FMIndex
from bwa_mapper.mapping import map_reads 
from bwa_mapper.snp import call_snps
from stats.coverage import compute_coverage
from stats.summary import calculate_snp_stats, calculate_mapping_stats

def landing_page():
    st.title("BWA/FM-index SNP Caller")

    reference_file = st.file_uploader(
        "Upload Reference FASTA", type=["fasta", "fa", "fna"], key="ref_uploader"
    )
    query_file = st.file_uploader(
        "Upload Query FASTQ", type=["fastq", "fq"], key="query_uploader"
    )
    col1, col2 = st.columns(2)
    with col1:
        min_depth = st.slider(
            "Minimum Depth", min_value=1, max_value=20, value=3, step=1
        )
    with col2:
        min_alt_frac = st.slider(
            "Minimum Alternate Allele Fraction", min_value=0.0, max_value=1.0, value=0.2, step=0.01
        )

    st.session_state["min_depth"] = min_depth
    st.session_state["min_alt_frac"] = min_alt_frac

    if st.button("Run Analysis"):
        if reference_file is None or query_file is None:
            st.warning("Please upload both reference and query files before running.")
            return

        reference_seq = ""
        for record in SeqIO.parse(StringIO(reference_file.getvalue().decode()), "fasta"):
            reference_seq += str(record.seq)
        st.session_state["reference_seq"] = reference_seq

        query_reads = {}
        for record in SeqIO.parse(StringIO(query_file.getvalue().decode()), "fastq"):
            query_reads[record.id] = str(record.seq)
        st.session_state["query_reads"] = query_reads

        st.info("Building BWT/FM-index...")
        start_bwt = time.time()
        fm = FMIndex(reference_seq)
        end_bwt = time.time()
        st.session_state["reference_bwt"] = fm.bwt
        st.session_state["runtime_bwt"] = f"{(end_bwt - start_bwt):.2f}"

        st.info("Mapping reads to reference...")
        start_mapping = time.time()
        alignments = map_reads(reference_seq, query_reads, fm)
        end_mapping = time.time()
        st.session_state["mapping"] = alignments
        st.session_state["runtime_mapping"] = f"{(end_mapping - start_mapping):.2f}"

        st.info("Calling SNPs...")
        start_snp = time.time()
        snp_list = call_snps(alignments, reference_seq, min_depth=min_depth, min_alt_frac=min_alt_frac)
        
        coverage = compute_coverage(alignments, len(reference_seq))
        st.session_state["coverage"] = coverage

        snp_stats = calculate_snp_stats(snp_list)
        st.session_state["snp_stats"] = snp_stats

        mapping_stats = calculate_mapping_stats(alignments, len(query_reads))
        st.session_state["mapping_stats"] = mapping_stats

        end_snp = time.time()
        st.session_state["results"] = snp_list
        st.session_state["runtime_snp"] = f"{(end_snp - start_snp):.2f}"

        st.session_state["view"] = "workspace"
        st.rerun()