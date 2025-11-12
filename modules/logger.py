# modules/logger.py
from datetime import datetime
import streamlit as st

def add_to_log(message):
    """Add message to log content and update display"""
    if 'log_content' not in st.session_state:
        st.session_state.log_content = ""
    
    timestamp = datetime.now().strftime("%H:%M:%S")
    st.session_state.log_content += f"[{timestamp}] {message}\n"
    
    # Update log display by recreating the text area
    with st.session_state.log_container:
        st.markdown("### 📝 Analysis Log")
        st.text_area(
            "Log Output", 
            value=st.session_state.log_content, 
            height=300,
            key=f"log_display_{datetime.now().timestamp()}"
        )