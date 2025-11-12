# modules/quantification_module.py
import streamlit as st
from pathlib import Path
from .utils import run_command, check_tool_available, validate_directory, find_files
from .config import Config
from .logger import add_to_log
import subprocess

# modules/quantification.py
def get_quantification_params():
    """Collect quantification parameters - ONLY ask for input directory if starting step"""
    st.markdown("#### 📁 Input Parameters")
    
    # Only ask for BAM directory if this is the starting step (Step 5)
    if st.session_state.start_step == 5:
        st.session_state.params['bam_dir'] = st.text_input(
            "**BAM Files Directory**",
            help="Directory containing BAM files for quantification",
            placeholder="/path/to/bam_files",
            key="bam_dir"
        )
        
        if st.session_state.params.get('bam_dir'):
            bam_path = Path(st.session_state.params['bam_dir'])
            if validate_directory(bam_path):
                bam_files = find_files(bam_path, "*.bam")
                if bam_files:
                    st.success(f"✅ Found {len(bam_files)} BAM files")
                else:
                    st.warning("⚠️ No BAM files found")
    else:
        # For downstream steps, auto-detect from alignment output
        st.info("🔄 BAM directory will be auto-detected from alignment output")
    
    # ALWAYS ask for annotation file
    st.session_state.params['annotation_file'] = st.text_input(
        "**Annotation File (GTF/GFF3)**",
        help="Path to gene annotation file",
        placeholder="/path/to/annotation.gtf",
        key="annotation_file"
    )
    
    if st.session_state.params.get('annotation_file'):
        anno_path = Path(st.session_state.params['annotation_file'])
        if not anno_path.exists():
            st.warning("⚠️ Annotation file not found at specified path")
    
    # Rest of the parameters (always ask for these)
    st.markdown("#### ⚙️ Quantification Parameters")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.session_state.params['quant_threads'] = st.number_input(
            "**Number of Threads**",
            min_value=2,
            max_value=200,
            value=16,
            step=2,
            help="Number of threads for StringTie",
            key="quant_threads"
        )
    
        st.session_state.params['min_transcript_cov'] = st.number_input(
            "**Minimum Transcript Coverage**",
            min_value=1.0,
            max_value=100.0,
            value=2.5,
            step=0.5,
            help="Minimum coverage for predicted transcripts",
            key="min_transcript_cov"
        )
        
        st.session_state.params['min_isoform_abundance'] = st.number_input(
            "**Minimum Isoform Abundance**",
            min_value=0.01,
            max_value=10.0,
            value=0.01,
            step=0.01,
            help="Minimum isoform abundance",
            key="min_isoform_abundance"
        )
    
    with col2:
        st.session_state.params['junction_cov'] = st.number_input(
            "**Junction Coverage**",
            min_value=1.0,
            max_value=10.0,
            value=1.0,
            step=0.5,
            help="Minimum junction coverage",
            key="junction_cov"
        )

        st.session_state.params['min_anchor_length'] = st.number_input(
            "**Minimum Anchor Length**",
            min_value=1,
            max_value=100,
            value=10,
            help="Minimum anchor length for junctions",
            key="min_anchor_length"
        )

def run_quantification_step(params, output_dir):
    """Run StringTie quantification - WITH PROPER AUTO-DETECTION"""
    try:
        # Check required tools
        tools = ["stringtie", "gffread"]
        missing_tools = []
        for tool in tools:
            if not check_tool_available(tool):
                missing_tools.append(tool)
        
        if missing_tools:
            add_to_log(f"❌ Required tools not found: {', '.join(missing_tools)}")
            return False
        
        # Create quantification directory
        quantification_dir = output_dir / "quantification"
        quantification_dir.mkdir(parents=True, exist_ok=True)
        
        # Validate inputs
        annotation_file = params.get('annotation_file')
        
        # Handle BAM directory - AUTO-DETECT from previous step when not starting step
        if st.session_state.start_step == 5:
            # Starting from quantification step, use provided BAM directory
            bam_dir = params.get('bam_dir')
        else:
            # Auto-detect from alignment output
            bam_dir = str(output_dir / "aligned_reads")
            add_to_log(f"🔄 Auto-detected BAM directory: {bam_dir}")
        
        add_to_log(f"🔍 Checking BAM directory: {bam_dir}")
        
        bam_path = Path(bam_dir)
        if not bam_path.exists():
            add_to_log(f"❌ BAM directory does not exist: {bam_dir}")
            return False
            
        if not annotation_file or not Path(annotation_file).exists():
            add_to_log("❌ Annotation file not found")
            return False
        
        # Check for BAM files
        bam_files = find_files(bam_path, "*.bam")
        if not bam_files:
            add_to_log("❌ No BAM files found in directory")
            return False
        
        add_to_log(f"📊 Found {len(bam_files)} BAM files for quantification")
        
        # Convert GFF3 to GTF if needed
        if annotation_file.lower().endswith(('.gff3', '.gff')):
            add_to_log("🔄 Converting GFF3 to GTF format...")
            gtf_file = convert_gff_to_gtf(annotation_file, quantification_dir)
            if not gtf_file:
                return False
        else:
            gtf_file = annotation_file
        
        # Run StringTie quantification for each BAM file
        success_count = 0
        threads = params.get('quant_threads', 16)
        mode = '-eB'  # Default mode for reference-based analysis
        min_cov = params.get('min_transcript_cov', 2.5)
        min_abund = params.get('min_isoform_abundance', 0.01)
        
        for bam_file in bam_files:
            sample_name = bam_file.stem.replace('_dedup', '').replace('_aligned', '')
            output_gtf = quantification_dir / f"{sample_name}.gtf"
            gene_abundance = quantification_dir / f"{sample_name}_gene_abundance.tab"
            
            cmd = (
                f"stringtie {bam_file} {mode} -G {gtf_file} "
                f"-A {gene_abundance} -o {output_gtf} -p {threads} "
                f"-c {min_cov} -f {min_abund}"
            )
            
            add_to_log(f"Processing {bam_file.name}...")
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                success_count += 1
                add_to_log(f"✅ Successfully processed {sample_name}")
            else:
                add_to_log(f"❌ Failed to process {sample_name}: {result.stderr}")
        
        # Generate count matrices using prepDE.py
        if success_count > 0:
            add_to_log("📈 Generating count matrices...")
            count_matrix_dir = quantification_dir / "count_matrix"
            count_matrix_dir.mkdir(exist_ok=True)
            
            # Create sample list for prepDE.py
            sample_list_file = quantification_dir / "sample_list.txt"
            with open(sample_list_file, 'w') as f:
                for bam_file in bam_files:
                    sample_name = bam_file.stem.replace('_dedup', '').replace('_aligned', '')
                    gtf_file = quantification_dir / f"{sample_name}.gtf"
                    if gtf_file.exists():
                        f.write(f"{sample_name}\t{gtf_file}\n")
            
            prepde_cmd = (
                f"python {Config.PREPDE_SCRIPT} -i {sample_list_file} "
                f"-g {count_matrix_dir}/gene_count_matrix.csv "
                f"-t {count_matrix_dir}/transcript_count_matrix.csv"
            )
            
            add_to_log(f"Running: {prepde_cmd}")
            result = subprocess.run(prepde_cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                add_to_log(f"✅ Quantification completed for {success_count}/{len(bam_files)} samples")
                add_to_log(f"📁 Output directory: {quantification_dir}")
                
                # Verify count matrices were created
                gene_matrix = count_matrix_dir / "gene_count_matrix.csv"
                transcript_matrix = count_matrix_dir / "transcript_count_matrix.csv"
                
                if gene_matrix.exists():
                    add_to_log(f"✅ Gene count matrix created: {gene_matrix}")
                if transcript_matrix.exists():
                    add_to_log(f"✅ Transcript count matrix created: {transcript_matrix}")
                    
                return True
            else:
                add_to_log(f"❌ Count matrix generation failed: {result.stderr}")
        
        add_to_log(f"❌ Quantification failed for {len(bam_files) - success_count} samples")
        return False
        
    except Exception as e:
        add_to_log(f"❌ Quantification step failed: {str(e)}")
        import traceback
        add_to_log(f"🔍 Detailed error: {traceback.format_exc()}")
        return False

def convert_gff_to_gtf(gff_file, output_dir):
    """Convert GFF3 file to GTF format"""
    gtf_file = output_dir / "annotation.gtf"
    
    cmd = f"gffread {gff_file} -T -o {gtf_file}"
    
    if run_command(cmd, "GFF to GTF Conversion"):
        return gtf_file
    return None