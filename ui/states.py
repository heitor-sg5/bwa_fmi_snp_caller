import streamlit as st

def ensure_session_defaults():
    st.session_state.setdefault("view", "landing")
    st.session_state.setdefault("results", [])
    st.session_state.setdefault("runtime_bwt", 0)
    st.session_state.setdefault("runtime_mapping", 0)
    st.session_state.setdefault("runtime_snp", 0)
    st.session_state.setdefault("reference_seq", "")
    st.session_state.setdefault("reference_bwt", "")
    st.session_state.setdefault("query_reads", {})
    st.session_state.setdefault("stats", {})
    st.session_state.setdefault("coverage", [])