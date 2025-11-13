# modules/trimming_module.py
import streamlit as st
from pathlib import Path
from .utils import validate_directory, find_files
from .config import Config
from .logger import add_to_log
import subprocess
import re
import os

# modules/trimming.py
def get_trimming_params():
    """Collect trimming parameters - ONLY ask for input directory if starting step"""
    st.markdown("#### 📁 Input Parameters")
    
    # Only ask for input directory if this is the starting step (Step 2)
    if st.session_state.start_step == 2:
        st.session_state.params['input_dir'] = st.text_input(
            "**Input Directory with FASTQ Files**",
            help="Directory containing raw FASTQ files for trimming",
            placeholder="/path/to/raw_fastq_files"
        )
        
        if st.session_state.params.get('input_dir'):
            input_path = Path(st.session_state.params['input_dir'])
            if validate_directory(input_path):
                fastq_files = find_files(input_path, "*.fastq.gz")
                if fastq_files:
                    st.success(f"✅ Found {len(fastq_files)} FASTQ files")
                else:
                    st.warning("⚠️ No FASTQ files found")
    else:
        # For downstream steps, auto-detect from QC output
        st.info("🔄 Input will be auto-detected from QC output")
    
    # ALWAYS ask for adapter file (required for trimming)
    st.markdown("#### 🔧 Adapter Selection")

    # Get the absolute path of the directory where trimming.py is located
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))

    # Path to the adapters folder (relative to this file)
    ADAPTER_DIR = os.path.join(BASE_DIR, '..', 'Trimmomatic-0.39', 'adapters')

    # Define available adapter files
    adapter_options = {
        "TruSeq3-PE (Paired-End)": os.path.join(ADAPTER_DIR, 'TruSeq3-PE.fa'),
        "TruSeq3-PE-2 (Paired-End)": os.path.join(ADAPTER_DIR, 'TruSeq3-PE-2.fa'),
        "TruSeq3-SE (Single-End)": os.path.join(ADAPTER_DIR, 'TruSeq3-SE.fa'),
        "TruSeq2-PE (Paired-End)": os.path.join(ADAPTER_DIR, 'TruSeq2-PE.fa'),
        "TruSeq2-SE (Single-End)": os.path.join(ADAPTER_DIR, 'TruSeq2-SE.fa'),
        "NexteraPE-PE": os.path.join(ADAPTER_DIR, 'NexteraPE-PE.fa')
    }
    
    # Auto-select based on mode, but allow override
    mode = st.session_state.params.get('mode', 'PE')
    default_adapter = "TruSeq3-PE (Paired-End)" if mode == "PE" else "TruSeq3-SE (Single-End)"
    
    selected_adapter = st.selectbox(
        "**Select Adapter File**",
        options=list(adapter_options.keys()),
        index=list(adapter_options.keys()).index(default_adapter) if default_adapter in adapter_options else 0,
        help="Select the appropriate adapter file for your sequencing protocol"
    )
    
    # Show custom path option
    use_custom = st.checkbox("Use custom adapter file path")
    
    if use_custom:
        # When using custom, clear the current adapter file and ask for new path
        custom_adapter_path = st.text_input(
            "**Custom Adapter File Path**",
            value="",
            help="Path to custom adapter FASTA file",
            placeholder="add/path/to/custom/adapter/file"
        )
        
        if custom_adapter_path:
            st.session_state.params['adapter_file'] = custom_adapter_path
            st.info(f"**Using custom adapter file:** {custom_adapter_path}")
        else:
            st.session_state.params['adapter_file'] = ""  # Clear if no custom path provided
            st.warning("⚠️ Please provide a custom adapter file path")
    else:
        # Use the selected adapter from dropdown
        st.session_state.params['adapter_file'] = adapter_options[selected_adapter]
        st.info(f"**Selected adapter file:** {st.session_state.params['adapter_file']}")
    
    # Validate adapter file exists (only if path is provided)
    if st.session_state.params.get('adapter_file'):
        adapter_path = Path(st.session_state.params['adapter_file'])
        if adapter_path.exists():
            st.success("✅ Adapter file found")
        else:
            st.warning("⚠️ Adapter file not found at specified path")
    else:
        st.warning("⚠️ No adapter file selected")
    
    # Rest of the parameters (always ask for these)
    st.markdown("#### ⚙️ Trimming Parameters")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.session_state.params['phred'] = st.selectbox(
            "**PHRED Encoding**", 
            ["33", "64"], 
            index=0,
            help="Quality score encoding (33 for Sanger, 64 for Illumina)"
        )
        
        st.session_state.params['trim_threads'] = st.number_input(
            "**Number of Threads**", 
            min_value=2, 
            max_value=200,
            value=16, 
            step=2,
            help="Number of threads for Trimmomatic"
        )
    
    with col2:
        st.session_state.params['leading_quality'] = st.number_input(
            "**LEADING Quality**",
            min_value=3, 
            max_value=200,
            value=3, 
            step=3,
            help="Remove leading low quality bases (below this quality)"
        )
        
        st.session_state.params['trailing_quality'] = st.number_input(
            "**TRAILING Quality**",
            min_value=3, 
            max_value=200,
            value=3, 
            step=3,
            help="Remove trailing low quality bases (below this quality)"
        )
    
    st.markdown("#### 🔧 Advanced Parameters")
    
    use_advanced = st.checkbox(
        "Use Advanced Trimming Parameters", 
        value=False,
        help="Enable sliding window trimming and minimum length filtering"
    )
    
    if use_advanced:
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.session_state.params['minlen'] = st.number_input(
                "**MINLEN Length**",
                min_value=3, 
                max_value=200,
                value=36, 
                step=1,
                help="Minimum length of reads to keep"
            )
        
        with col2:
            st.session_state.params['window_size'] = st.number_input(
                "**SLIDINGWINDOW Size**",
                min_value=2, 
                max_value=200,
                value=4, 
                step=1,
                help="Window size for sliding window trimming"
            )
        
        with col3:
            st.session_state.params['window_quality'] = st.number_input(
                "**SLIDINGWINDOW Quality**",
                min_value=5, 
                max_value=200,
                value=20, 
                step=5,
                help="Average quality required in sliding window"
            )
    else:
        # Set advanced parameters to None when not used
        st.session_state.params['minlen'] = None
        st.session_state.params['window_size'] = None
        st.session_state.params['window_quality'] = None

def find_pairs(fastq_files):
    """Find PE pairs in FASTQ files using more flexible pattern matching"""
    pairs = {}
    
    # First, identify all possible R1 files
    r1_files = []
    for fq in fastq_files:
        name = fq.name
        # More flexible pattern to match R1 files
        if re.search(r'(_R1_|_R1\.|_1_|_1\.)', name) or re.search(r'[._][Rr]1[._]', name):
            r1_files.append(fq)
    
    # Now try to find matching R2 files
    for r1_file in r1_files:
        r1_name = r1_file.name
        
        # Try different patterns to find the base name and R2 counterpart
        base_name = None
        r2_name = None
        
        # Pattern 1: _R1_ / _R2_
        if '_R1_' in r1_name:
            base_name = r1_name.replace('_R1_', '_')
            r2_name = r1_name.replace('_R1_', '_R2_')
        # Pattern 2: _1_ / _2_
        elif '_1_' in r1_name:
            base_name = r1_name.replace('_1_', '_')
            r2_name = r1_name.replace('_1_', '_2_')
        # Pattern 3: _R1. / _R2.
        elif '_R1.' in r1_name:
            base_name = r1_name.replace('_R1.', '.')
            r2_name = r1_name.replace('_R1.', '_R2.')
        # Pattern 4: _1. / _2.
        elif '_1.' in r1_name:
            base_name = r1_name.replace('_1.', '.')
            r2_name = r1_name.replace('_1.', '_2.')
        # Pattern 5: .R1. / .R2.
        elif '.R1.' in r1_name:
            base_name = r1_name.replace('.R1.', '.')
            r2_name = r1_name.replace('.R1.', '.R2.')
        # Pattern 6: .1. / .2.
        elif '.1.' in r1_name:
            base_name = r1_name.replace('.1.', '.')
            r2_name = r1_name.replace('.1.', '.2.')
        # Pattern 7: _R1.fastq / _R2.fastq
        elif re.search(r'_R1\.fastq', r1_name):
            base_name = re.sub(r'_R1\.fastq', '.fastq', r1_name)
            r2_name = r1_name.replace('_R1.fastq', '_R2.fastq')
        # Pattern 8: _1.fastq / _2.fastq
        elif re.search(r'_1\.fastq', r1_name):
            base_name = re.sub(r'_1\.fastq', '.fastq', r1_name)
            r2_name = r1_name.replace('_1.fastq', '_2.fastq')
        
        if base_name and r2_name:
            r2_path = r1_file.parent / r2_name
            if r2_path.exists():
                # Use the original R1 filename as base for output naming
                pairs[r1_name] = (r1_file, r2_path, r2_name)
            else:
                add_to_log(f"⚠️ Warning: no matching R2 found for {r1_name}")
    
    return pairs

def extract_sample_base(filename):
    """Extract clean sample base name by removing read identifiers and extensions"""
    # Remove common file extensions
    base = filename.replace('.fastq.gz', '').replace('.fq.gz', '').replace('.fastq', '').replace('.fq', '')
    
    # Remove common read identifier patterns (in order of specificity)
    patterns_to_remove = [
        r'_R1$', r'_R2$',        # _R1, _R2 at end
        r'_1$', r'_2$',          # _1, _2 at end  
        r'_R1_', r'_R2_',        # _R1_, _R2_ in middle
        r'_1_', r'_2_',          # _1_, _2_ in middle
        r'\.R1\.', r'\.R2\.',    # .R1., .R2. in middle
        r'\.1\.', r'\.2\.',      # .1., .2. in middle
        r'_R1\.', r'_R2\.',      # _R1., _R2. before extension
        r'_1\.', r'_2\.',        # _1., _2. before extension
    ]
    
    for pattern in patterns_to_remove:
        base = re.sub(pattern, '_', base)
    
    # Remove any trailing underscores or dots
    base = base.rstrip('_.')
    
    return base

def run_trimming_step(params, output_dir):
    """Run read trimming - PROCESS FILES INDIVIDUALLY"""
    try:
        trimmed_dir = output_dir / "trimmed_reads"
        trimmed_dir.mkdir(parents=True, exist_ok=True)
        
        # Validate inputs
        input_dir = params.get('input_dir')
        adapter_file = params.get('adapter_file')
        
        if not input_dir or not validate_directory(Path(input_dir)):
            add_to_log("❌ Invalid input directory")
            return False
            
        if not adapter_file or not Path(adapter_file).exists():
            add_to_log("❌ Adapter file not found")
            return False
        
        # Find input files
        fastq_files = find_files(Path(input_dir), "*.fastq.gz")
        if not fastq_files:
            add_to_log("❌ No FASTQ files found in input directory")
            return False
        
        add_to_log(f"✂️ Found {len(fastq_files)} FASTQ files for processing...")
        
        # Build Trimmomatic command based on mode
        mode = params.get('mode', 'PE')
        success_count = 0
        
        if mode == "PE":
            # Find paired-end files
            pairs = find_pairs(fastq_files)
            if not pairs:
                add_to_log("❌ No paired-end file pairs found")
                return False
            
            if params.get('minlen') or params.get('window_size'):
                add_to_log("🔧 Using advanced trimming parameters")
            else:
                add_to_log("🔧 Using basic trimming parameters (no advanced options)")
            
            add_to_log(f"🔍 Found {len(pairs)} paired-end file pairs")
            
            # Process each pair individually
            for r1_name, (r1_file, r2_file, r2_name) in pairs.items():
                # FIXED: Use extract_sample_base function instead of simple string replacement
                sample_base = extract_sample_base(r1_name)
                add_to_log(f"🔧 Processing pair: {r1_file.name} + {r2_file.name}")

                cmd = build_pe_trimming_command_for_pair(r1_file, r2_file, trimmed_dir, params, sample_base)

                add_to_log(f"Running: {cmd}")

                # Run command
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

                if result.returncode == 0:
                    success_count += 1
                    add_to_log(f"✅ Successfully trimmed {sample_base}")
                else:
                    add_to_log(f"❌ Failed to trim {sample_base}: {result.stderr}")
        else:
            # Single-end mode - process each file individually
            add_to_log(f"🔧 Processing {len(fastq_files)} single-end files")
            
            for fastq_file in fastq_files:
                # FIXED: Also use extract_sample_base for single-end files for consistency
                sample_base = extract_sample_base(fastq_file.name)
                add_to_log(f"🔧 Processing: {fastq_file.name}")
                
                cmd = build_se_trimming_command_for_file(fastq_file, trimmed_dir, params, sample_base)
                
                add_to_log(f"Running: {cmd}")
                
                # Run command
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                
                if result.returncode == 0:
                    success_count += 1
                    add_to_log(f"✅ Successfully trimmed {sample_base}")
                else:
                    add_to_log(f"❌ Failed to trim {sample_base}: {result.stderr}")
        
        # Check results
        if success_count > 0:
            add_to_log(f"✅ Trimming completed for {success_count} files")
            trimmed_files = list(trimmed_dir.glob("*.fastq.gz"))
            if trimmed_files:
                add_to_log(f"✅ Trimmed files generated: {len(trimmed_files)} files")
                add_to_log(f"📁 Output directory: {trimmed_dir}")
                return True
            else:
                add_to_log("❌ No trimmed output files generated")
                return False
        else:
            add_to_log("❌ Trimming failed for all files")
            return False
        
    except Exception as e:
        add_to_log(f"❌ Trimming step failed: {str(e)}")
        return False

def build_pe_trimming_command_for_pair(r1_file, r2_file, output_dir, params, clean_base):
    """Build paired-end trimming command for a specific file pair"""
    threads = params.get('trim_threads', 8)
    phred = params.get('phred', '33')
    leading = params.get('leading_quality', 3)
    trailing = params.get('trailing_quality', 3)
    minlen = params.get('minlen')
    window_size = params.get('window_size')
    window_quality = params.get('window_quality')
    adapter_file = params.get('adapter_file')
    
    # Base command without advanced parameters
    base_cmd = (
        f"java -jar {Config.TRIMMOMATIC_JAR} PE "
        f"-threads {threads} "
        f"-phred{phred} "
        f"{r1_file} {r2_file} "
        f"{output_dir}/{clean_base}_1_paired.fastq.gz {output_dir}/{clean_base}_1_unpaired.fastq.gz "
        f"{output_dir}/{clean_base}_2_paired.fastq.gz {output_dir}/{clean_base}_2_unpaired.fastq.gz "
        f"ILLUMINACLIP:{adapter_file}:2:30:10 "
        f"LEADING:{leading} "
        f"TRAILING:{trailing}"
    )
    
    # Conditionally add advanced parameters
    if window_size and window_quality:
        base_cmd += f" SLIDINGWINDOW:{window_size}:{window_quality}"
    
    if minlen:
        base_cmd += f" MINLEN:{minlen}"
    
    return base_cmd

def build_se_trimming_command_for_file(fastq_file, output_dir, params, sample_base):
    """Build single-end trimming command for a specific file"""
    threads = params.get('trim_threads', 8)
    phred = params.get('phred', '33')
    leading = params.get('leading_quality', 3)
    trailing = params.get('trailing_quality', 3)
    minlen = params.get('minlen')
    window_size = params.get('window_size')
    window_quality = params.get('window_quality')
    adapter_file = params.get('adapter_file')
    
    # Base command without advanced parameters
    base_cmd = (
        f"java -jar {Config.TRIMMOMATIC_JAR} SE "
        f"-threads {threads} "
        f"-phred{phred} "
        f"{fastq_file} "
        f"{output_dir}/{sample_base}_trimmed.fastq.gz "
        f"ILLUMINACLIP:{adapter_file}:2:30:10 "
        f"LEADING:{leading} "
        f"TRAILING:{trailing}"
    )
    
    # Conditionally add advanced parameters
    if window_size and window_quality:
        base_cmd += f" SLIDINGWINDOW:{window_size}:{window_quality}"
    
    if minlen:
        base_cmd += f" MINLEN:{minlen}"
    
    return base_cmd
