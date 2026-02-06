import streamlit as st

from ui.states import ensure_session_defaults 
from ui.landing import landing_page   
from ui.workspace import workspace 

def main():
    """Main application entry point."""
    st.set_page_config(
        page_title="BWA/FM-index SNP Caller",
        layout="wide",
        initial_sidebar_state="collapsed",
    )

    ensure_session_defaults()

    if st.session_state["view"] == "landing":
        landing_page()
    else:
        workspace()

if __name__ == "__main__":
    main()