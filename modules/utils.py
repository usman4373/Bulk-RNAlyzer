# modules/utils.py
import subprocess
from pathlib import Path
import streamlit as st


def check_tool_available(tool_name):
    """Check if a required tool is available in the system PATH"""
    try:
        cmd = [tool_name, "--version"]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="ignore",
            timeout=5
        )
        return result.returncode == 0
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return False

def run_command(cmd, description):
    """Run a shell command and handle errors"""
    st.info(f"🔄 Running: {description}")
    st.code(f"Command: {cmd}")
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        st.error(f"❌ {description} failed!")
        st.code(f"Error: {result.stderr}")
        return False
    
    st.success(f"✅ {description} completed successfully")
    return True

def setup_logging(output_dir):
    """Set up logging directory"""
    log_dir = output_dir / "logs"
    log_dir.mkdir(exist_ok=True)
    return log_dir

def find_files(directory, pattern):
    """Find files matching pattern in directory"""
    path = Path(directory)
    return list(path.glob(pattern))

def validate_directory(directory):
    """Validate that directory exists and is accessible"""
    path = Path(directory)
    return path.exists() and path.is_dir()

def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string into dictionary"""
    attributes = {}
    for field in attr_string.split(';'):
        field = field.strip()
        if field and ' ' in field:
            key, value = field.split(' ', 1)
            value = value.strip('"')
            attributes[key] = value
    return attributes