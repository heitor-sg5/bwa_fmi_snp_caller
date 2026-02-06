import streamlit as st
import pandas as pd
import plotly.express as px
import numpy as np

def workspace():
    st.title("SNP Analysis")

    col1, col2 = st.columns(2)
    with col1:
        if "reference_bwt" in st.session_state:
            st.download_button(
                label="Download Reference BWT",
                data=st.session_state["reference_bwt"],
                file_name="reference_bwt.txt",
            )

    with col2:
        if "results" in st.session_state:
            snp_df = pd.DataFrame(
                st.session_state["results"], columns=["Position", "Ref", "Alt", "Count", "Depth"]
            )
            st.download_button(
                label="Download SNP Table",
                data=snp_df.to_csv(index=False),
                file_name="snps.csv",
            )

    stats_col, snp_col = st.columns(2)

    with stats_col:
        st.subheader("Statistics")
        stats = st.session_state.get("stats", {})
        stats["BWT Runtime"] = st.session_state.get("runtime_bwt", None)
        stats["Mapping Runtime"] = st.session_state.get("runtime_mapping", None)
        stats["SNP Calling Runtime"] = st.session_state.get("runtime_snp", None)
        stats_df = pd.DataFrame(stats.items(), columns=["Metric", "Value"])
        st.table(stats_df)

    with snp_col:
        st.subheader("SNP Table")
        if "results" in st.session_state:
            snp_df = pd.DataFrame(
                st.session_state["results"], columns=["Position", "Ref", "Alt", "Count", "Depth"]
            )
            def mutation_type(row):
                transitions = {("A","G"), ("G","A"), ("C","T"), ("T","C")}
                if (row["Ref"], row["Alt"]) in transitions:
                    return "Transition"
                return "Transversion"
            snp_df["Type"] = snp_df.apply(mutation_type, axis=1)
            sort_col = st.selectbox("Sort by:", ["Position", "Type"], key="sort_col")
            sort_order = st.selectbox("Order:", ["Ascending", "Descending"], key="sort_order")
            ascending = sort_order == "Ascending"
            snp_df = snp_df.sort_values(by=sort_col, ascending=ascending)
            st.dataframe(snp_df)

    tabs = st.tabs(["Coverage", "Depth", "Density", "Mutation"])
    reference_length = len(st.session_state.get("reference_seq", ""))
    coverage = st.session_state.get("coverage", np.zeros(reference_length, dtype=int))

    with tabs[0]:
        st.subheader("Query Read Coverage")
        fig = px.line(y=coverage, labels={"x": "Reference Position", "y": "Coverage"})
        fig.update_layout(height=300, margin=dict(l=40, r=20, t=30, b=40))
        st.plotly_chart(fig, use_container_width=True)

    with tabs[1]:
        st.subheader("Read Depth Distribution")
        fig = px.histogram(x=coverage, nbins=50, labels={"x": "Depth", "y": "Frequency"})
        fig.update_layout(height=300, margin=dict(l=40, r=20, t=30, b=40))
        st.plotly_chart(fig, use_container_width=True)

    with tabs[2]:
        st.subheader("SNP Density per 10% of Reference")
        if "results" in st.session_state:
            snp_positions = [r[0] for r in st.session_state["results"]]
            bins = np.linspace(0, reference_length, 11)
            counts, _ = np.histogram(snp_positions, bins=bins)
            density_df = pd.DataFrame({
                "Reference Segment": [f"{i*10}-{(i+1)*10}%" for i in range(10)],
                "SNP Count": counts
            })
            fig = px.bar(density_df, x="Reference Segment", y="SNP Count")
            fig.update_layout(height=300, margin=dict(l=40, r=20, t=30, b=40))
            st.plotly_chart(fig, use_container_width=True)

    with tabs[3]:
        st.subheader("Mutation Heatmap")
        if "results" in st.session_state:
            bases = ["A", "C", "G", "T"]
            mutation_matrix = pd.DataFrame(0, index=bases, columns=bases)
            for row in st.session_state["results"]:
                ref, alt = row[1], row[2]
                if ref in bases and alt in bases:
                    mutation_matrix.loc[ref, alt] += 1
            fig = px.imshow(
                mutation_matrix,
                text_auto=True,
                labels=dict(x="Alt Base", y="Ref Base", color="Count"),
                aspect="equal",
            )
            fig.update_layout(height=400, margin=dict(l=60, r=40, t=40, b=60))
            st.plotly_chart(fig, use_container_width=True)

    if st.button("Return to start page"):
        st.session_state["view"] = "landing"
        st.rerun()