# modules/qc_module.py
import streamlit as st
from pathlib import Path
from .utils import find_files, validate_directory
from .logger import add_to_log
import subprocess

# modules/quality_control.py
def get_qc_params(qc_type="before_trim"):
    """Collect QC parameters - ONLY ask for input directory if starting step"""
    if qc_type == "before_trim":
        st.markdown("#### 📁 Input Parameters")
        
        # Only ask for input directory if this is the starting step (Step 1)
        if st.session_state.start_step == 1:
            st.session_state.params['input_dir'] = st.text_input(
                "**Input Directory with Raw FASTQ Files**",
                help="Directory containing raw FASTQ files (.fastq.gz)",
                placeholder="/path/to/raw_fastq_files",
                key="qc_before_input_dir"
            )
            
            if st.session_state.params.get('input_dir'):
                input_path = Path(st.session_state.params['input_dir'])
                if validate_directory(input_path):
                    fastq_files = find_files(input_path, "*.fastq.gz")
                    if fastq_files:
                        st.success(f"✅ Found {len(fastq_files)} FASTQ files")
                        sample_files = [f.name for f in fastq_files[:3]]
                        st.write(f"Sample files: {', '.join(sample_files)}")
                        if len(fastq_files) > 3:
                            st.write(f"... and {len(fastq_files) - 3} more files")
                    else:
                        st.warning("⚠️ No FASTQ files found in directory")
                else:
                    st.error("❌ Invalid directory path")
        else:
            # For downstream steps, auto-detect from previous step output
            st.info("🔄 Input will be auto-detected from previous step output")
    
    else:  # after_trim
        st.markdown("#### 📁 Input Parameters")
        
        # Only ask for input directory if this is the starting step (Step 3)
        if st.session_state.start_step == 3:
            st.session_state.params['qc_after_trim_input_dir'] = st.text_input(
                "**Input Directory with Trimmed FASTQ Files**",
                help="Directory containing trimmed FASTQ files",
                placeholder="/path/to/trimmed_fastq_files",
                key="qc_after_input_dir"
            )
            
            if st.session_state.params.get('qc_after_trim_input_dir'):
                input_path = Path(st.session_state.params['qc_after_trim_input_dir'])
                if validate_directory(input_path):
                    fastq_files = find_files(input_path, "*.fastq.gz")
                    if fastq_files:
                        st.success(f"✅ Found {len(fastq_files)} trimmed FASTQ files")
                        sample_files = [f.name for f in fastq_files[:3]]
                        st.write(f"Sample files: {', '.join(sample_files)}")
                        if len(fastq_files) > 3:
                            st.write(f"... and {len(fastq_files) - 3} more files")
                    else:
                        st.warning("⚠️ No trimmed FASTQ files found in directory")
                else:
                    st.error("❌ Invalid directory path")
        else:
            # For downstream steps, auto-detect from trimming output
            st.info("🔄 Input will be auto-detected from trimming output")
    
    st.markdown("#### ⚙️ QC Parameters")
    
    col1, = st.columns(1)
    
    with col1:
        st.session_state.params['qc_threads'] = st.number_input(
            "**FastQC Threads**",
            min_value=2,
            max_value=200,
            value=16,
            step=2,
            help="Number of threads for FastQC processing",
            key=f"qc_threads_{qc_type}"
        )

def run_qc_before_trim_step(params, output_dir):
    """Run FastQC quality control on raw reads"""
    return run_qc_step(params, output_dir, "before_trim")

def run_qc_after_trim_step(params, output_dir):
    """Run FastQC quality control on trimmed reads"""
    return run_qc_step(params, output_dir, "after_trim")

def run_qc_step(params, output_dir, qc_type="before_trim"):
    """Run FastQC quality control - handles both before and after trim with proper auto-detection"""
    try:
        qc_dir = output_dir / "quality_control"
        
        if qc_type == "before_trim":
            qc_output_dir = qc_dir / "qc_before_trim"
            # Use the input directory from params (user provided for starting step)
            input_dir = params.get('input_dir')
            log_prefix = "Before trim"
        else:  # after_trim
            qc_output_dir = qc_dir / "qc_after_trim"
            # Handle auto-detection for step 3 when starting from step 1
            if st.session_state.start_step == 3:
                # Starting from step 3, use manual input
                input_dir = params.get('qc_after_trim_input_dir')
            else:
                # Auto-detect from trimming output when starting from earlier steps
                input_dir = str(output_dir / "trimmed_reads")
            log_prefix = "After trim"
        
        qc_output_dir.mkdir(parents=True, exist_ok=True)
        
        if not input_dir:
            add_to_log(f"❌ {log_prefix}: Input directory not specified")
            return False
        
        input_path = Path(input_dir)
        if not validate_directory(input_path):
            add_to_log(f"❌ {log_prefix}: Invalid input directory: {input_dir}")
            return False
        
        # Find FASTQ files
        fastq_files = find_files(input_path, "*.fastq.gz")
        if not fastq_files:
            add_to_log(f"❌ {log_prefix}: No FASTQ files found in input directory: {input_dir}")
            return False
        
        add_to_log(f"📊 {log_prefix}: Processing {len(fastq_files)} FASTQ files from {input_dir}...")
        
        # Build FastQC command
        quiet_flag = "--quiet" if params.get('fastqc_quiet', True) else ""
        threads = params.get('qc_threads', 4)
        
        cmd = f"fastqc {quiet_flag} --threads {threads} --outdir {qc_output_dir} {input_dir}/*.fastq.gz"
        
        add_to_log(f"Running: {cmd}")
        
        # Run command and capture output
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            add_to_log(f"✅ {log_prefix}: FastQC completed successfully")
            # Check if output was generated
            qc_outputs = list(qc_output_dir.glob("*_fastqc.*"))
            if qc_outputs:
                add_to_log(f"✅ {log_prefix}: QC reports generated: {len(qc_outputs)} files")
                add_to_log(f"📁 Reports saved to: {qc_output_dir}")
                return True
            else:
                add_to_log(f"❌ {log_prefix}: No QC output files generated")
                return False
        else:
            add_to_log(f"❌ {log_prefix}: FastQC failed with error: {result.stderr}")
            return False
        
    except Exception as e:
        add_to_log(f"❌ {log_prefix}: QC step failed: {str(e)}")
        return False