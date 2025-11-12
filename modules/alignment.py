# modules/alignment_module.py
import streamlit as st
from pathlib import Path
import re
import subprocess
from datetime import datetime
from .utils import check_tool_available, validate_directory, find_files
from .logger import add_to_log

# modules/alignment.py
def get_alignment_params():
    """Collect alignment parameters - ONLY ask for input directory if starting step"""
    st.markdown("#### 📁 Input Parameters")
    
    # Only ask for input directory if this is the starting step (Step 4)
    if st.session_state.start_step == 4:
        st.session_state.params['input_dir'] = st.text_input(
            "**Input Directory with Trimmed FASTQ Files**",
            help="Directory containing trimmed FASTQ files",
            placeholder="/path/to/trimmed_fastq_files",
            key="alignment_input_dir"
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
        # For downstream steps, auto-detect from trimming output
        st.info("🔄 Input will be auto-detected from trimming output")
    
    # ALWAYS ask for index base (required for alignment)
    st.session_state.params['index_base'] = st.text_input(
        "**HISAT2 Index Base Path**",
        help="Base path to HISAT2 index files (without .ht2 extension)",
        placeholder="/path/to/genome_index",
        key="index_base"
    )
    
    # Rest of the parameters (always ask for these)
    st.markdown("#### ⚙️ Alignment Parameters")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.session_state.params['align_threads'] = st.number_input(
            "**Number of Threads**",
            min_value=2,
            max_value=200,
            value=16,
            step=2,
            help="Number of threads for HISAT2 alignment",
            key="align_threads"
        )
        
        st.session_state.params['phred_alignment'] = st.selectbox(
            "**PHRED Encoding for Alignment**",
            ["33", "64"],
            index=0,
            help="Quality score encoding for alignment",
            key="phred_alignment"
        )

        st.session_state.params['min_intron_length'] = st.number_input(
            "**Minimum Intron Length**",
            min_value=10,
            max_value=1000,
            value=20,
            step=10,
            help="Minimum intron length",
            key="min_intron_length"
        )
    
    with col2:
        st.session_state.params['max_alignments'] = st.number_input(
            "**Maximum Alignments per Read**",
            min_value=1,
            max_value=1000,
            value=5,
            step=1,
            help="Maximum number of alignments to report per read",
            key="max_alignments"
        )

        st.session_state.params['max_intron_length'] = st.number_input(
            "**Maximum Intron Length**",
            min_value=1000,
            max_value=10000000,
            value=500000,
            step=5000,
            help="Maximum intron length",
            key="max_intron_length"
        )
        
        st.session_state.params['keep_intermediate'] = st.checkbox(
            "**Keep Intermediate Files**", 
            value=False,
            help="Keep intermediate SAM and sorted BAM files",
            key="keep_intermediate"
        )

def find_paired_files(input_dir):
    """Find paired-end files with _paired.fastq.gz pattern - FIXED VERSION"""
    paired_files = []
    
    # Get all paired files
    all_paired_files = sorted(Path(input_dir).glob("*_paired.fastq.gz"))
    
    # Group by sample name
    sample_dict = {}
    
    for file_path in all_paired_files:
        filename = file_path.name
        
        # Extract sample name and read type using more flexible patterns
        base_patterns = [
            r'(.+)_[12]_paired\.fastq\.gz',
            r'(.+)_R[12]_paired\.fastq\.gz',
            r'(.+)_[12]\.fastq\.gz_paired',
            r'(.+)_R[12]\.fastq\.gz_paired'
        ]
        
        sample_name = None
        read_type = None
        
        for pattern in base_patterns:
            match = re.match(pattern, filename)
            if match:
                sample_name = match.group(1)
                # Determine read type based on the actual matched pattern
                if '_R1' in filename or '_1' in filename:
                    # More precise detection for R1 files
                    if re.search(r'_R1_|_R1\.|_1_|_1\.', filename):
                        read_type = "R1"
                    elif re.search(r'_R1|_1', filename):
                        # Check if it's actually R1 by looking for the pattern before _paired
                        if '_1_paired' in filename or '_R1_paired' in filename:
                            read_type = "R1"
                elif '_R2' in filename or '_2' in filename:
                    # More precise detection for R2 files  
                    if re.search(r'_R2_|_R2\.|_2_|_2\.', filename):
                        read_type = "R2"
                    elif re.search(r'_R2|_2', filename):
                        if '_2_paired' in filename or '_R2_paired' in filename:
                            read_type = "R2"
                break
        
        # If pattern matching failed, try simpler approach
        if not sample_name:
            # Extract base name by removing known suffixes
            base_name = filename.replace('_paired.fastq.gz', '')
            if '_R1' in base_name:
                sample_name = base_name.replace('_R1', '')
                read_type = "R1"
            elif '_R2' in base_name:
                sample_name = base_name.replace('_R2', '')
                read_type = "R2" 
            elif '_1' in base_name:
                sample_name = base_name.replace('_1', '')
                read_type = "R1"
            elif '_2' in base_name:
                sample_name = base_name.replace('_2', '')
                read_type = "R2"
            else:
                # If we still can't determine, skip this file
                add_to_log(f"⚠️ Warning: Cannot determine sample name and read type for {filename}. Skipping.")
                continue
        
        if sample_name and read_type:
            if sample_name not in sample_dict:
                sample_dict[sample_name] = {"R1": None, "R2": None}
            
            sample_dict[sample_name][read_type] = file_path
            add_to_log(f"🔍 Found {read_type} file for sample {sample_name}: {filename}")
    
    # Only include samples that have both R1 and R2
    for sample_name, reads in sample_dict.items():
        if reads["R1"] and reads["R2"]:
            paired_files.append((sample_name, reads["R1"], reads["R2"]))
            add_to_log(f"✅ Complete pair found for sample: {sample_name}")
        else:
            add_to_log(f"⚠️ Incomplete pair for sample {sample_name}: R1={reads['R1'] is not None}, R2={reads['R2'] is not None}")
    
    return paired_files

def find_single_end_files(input_dir):
    """Find single-end files with _trimmed.fastq.gz pattern"""
    single_files = []
    
    input_path = Path(input_dir)
    
    # First, look for files with _trimmed.fastq.gz pattern (from trimming step)
    trimmed_files = sorted(input_path.glob("*_trimmed.fastq.gz"))
    
    for file_path in trimmed_files:
        filename = file_path.name
        
        # Skip files that look like paired-end
        if any(pe_indicator in filename.lower() for pe_indicator in ['_r1', '_r2', '_1', '_2']):
            continue
            
        # Extract sample name by removing _trimmed.fastq.gz
        sample_name = filename.replace('_trimmed.fastq.gz', '')
        single_files.append((sample_name, file_path))
        add_to_log(f"🔍 Found single-end trimmed file: {filename} → sample: {sample_name}")
    
    # If no trimmed files found, look for any fastq files that don't look paired
    if not single_files:
        add_to_log("⚠️ No specifically named single-end files found. Looking for any FASTQ files...")
        all_files = sorted(input_path.glob("*.fastq.gz"))
        for file_path in all_files:
            filename = file_path.name
            
            # Skip obvious paired-end files
            if any(pe_indicator in filename.lower() for pe_indicator in ['_r1', '_r2', '_1', '_2']):
                continue
                
            sample_name = filename.replace('.fastq.gz', '')
            single_files.append((sample_name, file_path))
            add_to_log(f"🔍 Found potential single-end file: {filename} → sample: {sample_name}")
    
    add_to_log(f"📁 Total single-end files found: {len(single_files)}")
    return single_files

def debug_file_patterns(input_dir):
    """Debug function to see what files are actually present"""
    input_path = Path(input_dir)
    add_to_log(f"🔍 Debug: Looking for files in {input_dir}")
    
    # List all fastq files
    all_files = sorted(input_path.glob("*.fastq.gz"))
    add_to_log(f"📁 All FASTQ files found: {len(all_files)}")
    for f in all_files:
        add_to_log(f"   - {f.name}")
    
    # List paired files specifically
    paired_files = sorted(input_path.glob("*_paired.fastq.gz"))
    add_to_log(f"📁 Paired FASTQ files found: {len(paired_files)}")
    for f in paired_files:
        add_to_log(f"   - {f.name}")

def check_hisat2_index(index_base):
    """Check if HISAT2 index files exist"""
    index_extensions = ['.1.ht2', '.2.ht2', '.3.ht2', '.4.ht2', '.5.ht2', '.6.ht2', 
                       '.7.ht2', '.8.ht2', '.rev.1.ht2', '.rev.2.ht2']
    
    index_files_exist = []
    for ext in index_extensions:
        index_file = Path(f"{index_base}{ext}")
        if index_file.exists():
            index_files_exist.append(True)
        else:
            index_files_exist.append(False)
    
    # Check if at least the basic index files exist
    if sum(index_files_exist) >= 6:
        return True
    else:
        add_to_log("❌ HISAT2 index files not found or incomplete.")
        add_to_log(f"Expected index base: {index_base}")
        add_to_log("Looking for files with extensions: .1.ht2, .2.ht2, ..., .8.ht2, .rev.1.ht2, .rev.2.ht2")
        return False

def run_samtools_sort(sam_file, sorted_bam_file, threads):
    """Convert SAM to sorted BAM using samtools"""
    # samtools view -bS sample_aligned.sam | samtools sort -o sample_aligned.sorted.bam -@ threads
    view_cmd = ["samtools", "view", "-bS", str(sam_file)]
    sort_cmd = ["samtools", "sort", "-o", str(sorted_bam_file), "-@", str(threads)]
    
    add_to_log(f"   Converting SAM to sorted BAM...")
    view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)
    sort_process = subprocess.Popen(sort_cmd, stdin=view_process.stdout)
    view_process.stdout.close()
    sort_process.communicate()
    
    if sort_process.returncode == 0:
        add_to_log(f"   ✅ Successfully created sorted BAM: {sorted_bam_file.name}")
        # Remove the intermediate SAM file
        sam_file.unlink()
        add_to_log(f"   ✅ Removed intermediate SAM file")
        return True
    else:
        add_to_log(f"   ❌ Error converting SAM to sorted BAM")
        return False

def run_sambamba_dedup(sorted_bam_file, dedup_bam_file, threads):
    """Perform deduplication using sambamba"""
    # sambamba markdup sample_aligned.sorted.bam sample_dedup.bam -t threads
    cmd = ["sambamba", "markdup", str(sorted_bam_file), str(dedup_bam_file), "-t", str(threads)]
    
    add_to_log(f"   Performing deduplication...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        add_to_log(f"   ✅ Successfully created deduplicated BAM: {dedup_bam_file.name}")
        return True
    else:
        add_to_log(f"   ❌ Error during deduplication")
        add_to_log(f"   STDERR: {result.stderr}")
        return False

def create_bam_index(bam_file):
    """Create BAM index file using samtools"""
    cmd = ["samtools", "index", str(bam_file)]
    
    add_to_log(f"   Creating BAM index...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        add_to_log(f"   ✅ Successfully created BAM index: {bam_file.name}.bai")
        return True
    else:
        add_to_log(f"   ❌ Error creating BAM index")
        add_to_log(f"   STDERR: {result.stderr}")
        return False

def extract_alignment_stats(hisat2_output):
    """Extract alignment statistics from HISAT2 output"""
    stats_lines = []
    for line in hisat2_output.split('\n'):
        if line.strip() and ('reads' in line.lower() or 'aligned' in line.lower() or '%' in line):
            stats_lines.append(line.strip())
    return '\n'.join(stats_lines)

def write_alignment_summary_file(summary_data, output_dir):
    """Write comprehensive summary file"""
    summary_file = output_dir / "alignment_summary.txt"
    
    with open(summary_file, 'w') as f:
        # Header section
        f.write("=" * 80 + "\n")
        f.write("HISAT2 ALIGNMENT SUMMARY REPORT\n")
        f.write("=" * 80 + "\n")
        f.write(f"Generated on: {summary_data['timestamp']}\n")
        f.write(f"Alignment mode: {summary_data['alignment_mode']}\n")
        f.write(f"Threads used: {summary_data['threads']}\n")
        f.write(f"Reference index: {summary_data['reference_index']}\n")
        f.write(f"Input directory: {summary_data['input_dir']}\n")
        f.write(f"Output directory: {summary_data['output_dir']}\n")
        f.write("=" * 80 + "\n\n")
        
        # Summary statistics
        f.write(f"TOTAL SAMPLES PROCESSED: {summary_data['total_samples']}\n")
        f.write(f"SUCCESSFULLY PROCESSED: {summary_data['successful_samples']}\n")
        f.write(f"FAILED: {summary_data['failed_samples']}\n")
        f.write("\n" + "=" * 80 + "\n\n")
        
        # Detailed sample information
        for sample in summary_data['samples']:
            f.write(f"SAMPLE: {sample['name']}\n")
            f.write("-" * 50 + "\n")
            
            # Input files
            f.write("Input files:\n")
            if summary_data['alignment_mode'] == 'Paired-End':
                f.write(f"  R1: {sample['r1_file']}\n")
                f.write(f"  R2: {sample['r2_file']}\n")
            else:
                f.write(f"  Single-end: {sample['input_file']}\n")
            
            # Output files
            f.write("Output files:\n")
            f.write(f"  Deduplicated BAM: {sample['dedup_bam']}\n")
            f.write(f"  BAM index: {sample['bam_index']}\n")
            f.write(f"Status: {sample['status']}\n\n")
            
            # Alignment statistics
            if sample['status'] == 'SUCCESS' and sample['alignment_stats']:
                f.write("ALIGNMENT STATISTICS:\n")
                f.write("-" * 30 + "\n")
                f.write(sample['alignment_stats'] + "\n")
            elif sample['status'] == 'FAILED':
                f.write("ERROR INFORMATION:\n")
                f.write("-" * 30 + "\n")
                f.write(f"Error during: {sample['error_step']}\n")
                if sample['error_message']:
                    f.write(f"Error details: {sample['error_message'][:500]}...\n")  # Truncate long error messages
            
            f.write("=" * 80 + "\n\n")
    
    add_to_log(f"✅ Summary report saved to: {summary_file}")

def run_alignment_step(params, output_dir):
    """Run HISAT2 alignment with detailed sample processing - WITH PROPER MODE DETECTION"""
    try:
        # Check required tools
        tools = ["hisat2", "samtools", "sambamba"]
        missing_tools = []
        for tool in tools:
            if not check_tool_available(tool):
                missing_tools.append(tool)
        
        if missing_tools:
            add_to_log(f"❌ Required tools not found: {', '.join(missing_tools)}")
            return False
        
        aligned_dir = output_dir / "aligned_reads"
        aligned_dir.mkdir(parents=True, exist_ok=True)
        
        # Handle input directory
        if st.session_state.start_step == 4:
            # Starting from alignment step, use provided input directory
            input_dir = params.get('input_dir')
        else:
            # Auto-detect from trimming output
            input_dir = str(output_dir / "trimmed_reads")
            add_to_log(f"🔄 Auto-detected trimmed reads directory: {input_dir}")

        # Ensure mode is properly set for all cases
        if 'mode' not in st.session_state.params:
            st.session_state.params['mode'] = mode

        index_base = params.get('index_base')
        
        if not input_dir or not validate_directory(Path(input_dir)):
            add_to_log(f"❌ Invalid input directory: {input_dir}")
            return False
            
        if not index_base:
            add_to_log("❌ HISAT2 index base not specified")
            return False
        
        # Debug: show what files are actually present
        debug_file_patterns(input_dir)

        # Check if index files exist
        if not check_hisat2_index(index_base):
            return False
        
        # Get alignment parameters - USE THE MODE FROM USER SELECTION
        threads = params.get('align_threads', 8)
        
        # CRITICAL FIX: Get mode from params, with proper fallback
        mode = st.session_state.params.get('mode', params.get('mode', 'PE'))
        add_to_log(f"🔧 Using sequencing mode: {mode}")
        
        keep_intermediate = params.get('keep_intermediate', False)
        max_alignments = params.get('max_alignments', 10)
        min_intron = params.get('min_intron_length', 20)
        max_intron = params.get('max_intron_length', 500000)
        
        add_to_log(f"🧬 Starting alignment with HISAT2 in {mode} mode...")
        
        # Prepare summary data structure
        summary_data = {
            'timestamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'alignment_mode': 'Paired-End' if mode == 'PE' else 'Single-End',
            'threads': threads,
            'reference_index': index_base,
            'input_dir': str(input_dir),
            'output_dir': str(aligned_dir),
            'total_samples': 0,
            'successful_samples': 0,
            'failed_samples': 0,
            'samples': []
        }

        # Find appropriate files based on the actual mode
        if mode == "PE":
            # Find paired-end files
            paired_samples = find_paired_files(input_dir)
            
            if not paired_samples:
                add_to_log("❌ No paired-end file pairs found!")
                add_to_log("Looking for files with pattern: *_[1|R1]_paired.fastq.gz and *_[2|R2]_paired.fastq.gz")
                return False
            
            summary_data['total_samples'] = len(paired_samples)
            
            add_to_log(f"Found {len(paired_samples)} paired-end samples:")
            for sample_name, r1, r2 in paired_samples:
                add_to_log(f"  {sample_name}:")
                add_to_log(f"    R1: {r1.name}")
                add_to_log(f"    R2: {r2.name}")
            
            add_to_log(f"Starting alignment of {len(paired_samples)} samples...")
            
            for sample_name, r1_file, r2_file in paired_samples:
                # Initialize sample data for summary
                sample_data = {
                    'name': sample_name,
                    'r1_file': r1_file.name,
                    'r2_file': r2_file.name,
                    'status': 'FAILED',
                    'alignment_stats': '',
                    'dedup_bam': '',
                    'bam_index': '',
                    'error_step': '',
                    'error_message': ''
                }
                
                # Define output files
                temp_sam = aligned_dir / f"{sample_name}_aligned.sam"
                sorted_bam = aligned_dir / f"{sample_name}_aligned.sorted.bam"
                dedup_bam = aligned_dir / f"{sample_name}_dedup.bam"
                bam_index = aligned_dir / f"{sample_name}_dedup.bam.bai"
                
                sample_data['dedup_bam'] = dedup_bam.name
                sample_data['bam_index'] = bam_index.name
                
                # Build HISAT2 command with all parameters
                cmd = [
                    "hisat2",
                    "-p", str(threads),
                    "-x", index_base,
                    "-1", str(r1_file),
                    "-2", str(r2_file),
                    "-S", str(temp_sam),
                    "-k", str(max_alignments),
                    "--min-intronlen", str(min_intron),
                    "--max-intronlen", str(max_intron)
                ]
                
                add_to_log(f"▶ Processing sample: {sample_name}")
                add_to_log(f"   Read 1: {r1_file.name}")
                add_to_log(f"   Read 2: {r2_file.name}")
                add_to_log(f"   Intermediate files: {temp_sam.name} → {sorted_bam.name} → {dedup_bam.name}")
                
                # Run HISAT2
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                # Capture alignment statistics
                alignment_stats = extract_alignment_stats(result.stdout)
                sample_data['alignment_stats'] = alignment_stats
                
                # Print HISAT2 output
                if result.stdout:
                    add_to_log("HISAT2 ALIGNMENT STATISTICS:")
                    for line in result.stdout.split('\n'):
                        if "aligned" in line and ("%" in line or "overall" in line.lower()):
                            add_to_log(f"   {line.strip()}")
                
                if result.returncode == 0:
                    add_to_log(f"   ✅ Successfully aligned {sample_name}")
                    
                    # Convert SAM to sorted BAM
                    if run_samtools_sort(temp_sam, sorted_bam, threads):
                        # Perform deduplication
                        if run_sambamba_dedup(sorted_bam, dedup_bam, threads):
                            # Create BAM index
                            if create_bam_index(dedup_bam):
                                add_to_log(f"   ✅ Successfully processed {sample_name}")
                                sample_data['status'] = 'SUCCESS'
                                summary_data['successful_samples'] += 1
                                
                                # Remove intermediate sorted BAM if requested
                                if not keep_intermediate:
                                    sorted_bam.unlink()
                                    add_to_log(f"   ✅ Removed intermediate sorted BAM file")
                                else:
                                    add_to_log(f"   💾 Kept intermediate sorted BAM file: {sorted_bam.name}")
                            else:
                                add_to_log(f"   ❌ Failed BAM indexing for {sample_name}")
                                sample_data['error_step'] = 'BAM indexing'
                                sample_data['error_message'] = 'Failed to create BAM index'
                                summary_data['failed_samples'] += 1
                        else:
                            add_to_log(f"   ❌ Failed deduplication for {sample_name}")
                            sample_data['error_step'] = 'Deduplication'
                            sample_data['error_message'] = 'Failed during sambamba markdup'
                            summary_data['failed_samples'] += 1
                    else:
                        add_to_log(f"   ❌ Failed BAM conversion for {sample_name}")
                        sample_data['error_step'] = 'BAM conversion'
                        sample_data['error_message'] = 'Failed during samtools sort'
                        summary_data['failed_samples'] += 1
                else:
                    add_to_log(f"❌ Error aligning {sample_name}")
                    add_to_log(f"STDERR: {result.stderr}")
                    sample_data['error_step'] = 'HISAT2 alignment'
                    sample_data['error_message'] = result.stderr
                    summary_data['failed_samples'] += 1
                
                # Add sample data to summary
                summary_data['samples'].append(sample_data)

        else:  # SE mode
            # Find single-end files
            single_samples = find_single_end_files(input_dir)
            
            if not single_samples:
                add_to_log("❌ No single-end files found!")
                add_to_log("Looking for files with pattern: *_trimmed.fastq.gz")
                return False
            
            summary_data['total_samples'] = len(single_samples)
            
            add_to_log(f"Found {len(single_samples)} single-end samples:")
            for sample_name, file_path in single_samples:
                add_to_log(f"  {sample_name}: {file_path.name}")
            
            add_to_log(f"Starting alignment of {len(single_samples)} samples...")
            
            for sample_name, input_file in single_samples:
                # Initialize sample data for summary
                sample_data = {
                    'name': sample_name,
                    'input_file': input_file.name,
                    'status': 'FAILED',
                    'alignment_stats': '',
                    'dedup_bam': '',
                    'bam_index': '',
                    'error_step': '',
                    'error_message': ''
                }
                
                # Define output files
                temp_sam = aligned_dir / f"{sample_name}_aligned.sam"
                sorted_bam = aligned_dir / f"{sample_name}_aligned.sorted.bam"
                dedup_bam = aligned_dir / f"{sample_name}_dedup.bam"
                bam_index = aligned_dir / f"{sample_name}_dedup.bam.bai"
                
                sample_data['dedup_bam'] = dedup_bam.name
                sample_data['bam_index'] = bam_index.name
                
                # Build HISAT2 command with all parameters
                cmd = [
                    "hisat2",
                    "-p", str(threads),
                    "-x", index_base,
                    "-U", str(input_file),
                    "-S", str(temp_sam),
                    "-k", str(max_alignments),
                    "--min-intronlen", str(min_intron),
                    "--max-intronlen", str(max_intron)
                ]
                
                add_to_log(f"▶ Processing sample: {sample_name}")
                add_to_log(f"   Input: {input_file.name}")
                add_to_log(f"   Intermediate files: {temp_sam.name} → {sorted_bam.name} → {dedup_bam.name}")
                
                # Run HISAT2
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                # Capture alignment statistics
                alignment_stats = extract_alignment_stats(result.stdout)
                sample_data['alignment_stats'] = alignment_stats
                
                # Print HISAT2 output
                if result.stdout:
                    add_to_log("HISAT2 ALIGNMENT STATISTICS:")
                    for line in result.stdout.split('\n'):
                        if "aligned" in line and ("%" in line or "overall" in line.lower()):
                            add_to_log(f"   {line.strip()}")
                
                if result.returncode == 0:
                    add_to_log(f"   ✅ Successfully aligned {sample_name}")
                    
                    # Convert SAM to sorted BAM
                    if run_samtools_sort(temp_sam, sorted_bam, threads):
                        # Perform deduplication
                        if run_sambamba_dedup(sorted_bam, dedup_bam, threads):
                            # Create BAM index
                            if create_bam_index(dedup_bam):
                                add_to_log(f"   ✅ Successfully processed {sample_name}")
                                sample_data['status'] = 'SUCCESS'
                                summary_data['successful_samples'] += 1
                                
                                # Remove intermediate sorted BAM if requested
                                if not keep_intermediate:
                                    sorted_bam.unlink()
                                    add_to_log(f"   ✅ Removed intermediate sorted BAM file")
                                else:
                                    add_to_log(f"   💾 Kept intermediate sorted BAM file: {sorted_bam.name}")
                            else:
                                add_to_log(f"   ❌ Failed BAM indexing for {sample_name}")
                                sample_data['error_step'] = 'BAM indexing'
                                sample_data['error_message'] = 'Failed to create BAM index'
                                summary_data['failed_samples'] += 1
                        else:
                            add_to_log(f"   ❌ Failed deduplication for {sample_name}")
                            sample_data['error_step'] = 'Deduplication'
                            sample_data['error_message'] = 'Failed during sambamba markdup'
                            summary_data['failed_samples'] += 1
                    else:
                        add_to_log(f"   ❌ Failed BAM conversion for {sample_name}")
                        sample_data['error_step'] = 'BAM conversion'
                        sample_data['error_message'] = 'Failed during samtools sort'
                        summary_data['failed_samples'] += 1
                else:
                    add_to_log(f"❌ Error aligning {sample_name}")
                    add_to_log(f"STDERR: {result.stderr}")
                    sample_data['error_step'] = 'HISAT2 alignment'
                    sample_data['error_message'] = result.stderr
                    summary_data['failed_samples'] += 1
                
                # Add sample data to summary
                summary_data['samples'].append(sample_data)

        # Write summary file
        write_alignment_summary_file(summary_data, aligned_dir)
        
        add_to_log(f"✅ Alignment completed!")
        add_to_log(f"Final deduplicated BAM files saved to: {aligned_dir}")
        add_to_log(f"Summary report saved to: {aligned_dir}/alignment_summary.txt")
        
        if summary_data['successful_samples'] > 0:
            return True
        else:
            add_to_log("❌ No samples were successfully processed")
            return False
        
    except Exception as e:
        add_to_log(f"❌ Alignment step failed: {str(e)}")
        return False
