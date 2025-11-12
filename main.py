# streamlit_app.py
import streamlit as st
import sys
from pathlib import Path
import subprocess
from pathlib import Path
import base64

# Add modules to path
sys.path.append('modules')

# Import modular components
from modules.index_genome import run_index_genome_step, get_index_genome_params, get_index_advanced_params, validate_index_genome_params
from modules.utils import check_tool_available
from modules.quality_control import run_qc_before_trim_step, run_qc_after_trim_step, get_qc_params
from modules.trimming import run_trimming_step, get_trimming_params
from modules.alignment import run_alignment_step, get_alignment_params
from modules.quantification import run_quantification_step, get_quantification_params
from modules.deseq2 import run_deseq2_step, get_deseq2_params
from modules.logger import add_to_log
from modules.config import Config

def init_session_state():
    """Initialize session state variables"""
    if 'current_page' not in st.session_state:
        st.session_state.current_page = "home"
    if 'selected_analysis' not in st.session_state:
        st.session_state.selected_analysis = None
    if 'start_step' not in st.session_state:
        st.session_state.start_step = None
    if 'params' not in st.session_state:
        st.session_state.params = {}
    if 'pipeline_status' not in st.session_state:
        st.session_state.pipeline_status = {}
    if 'output_dir' not in st.session_state:
        st.session_state.output_dir = None
    if 'analysis_running' not in st.session_state:
        st.session_state.analysis_running = False
    if 'log_container' not in st.session_state:
        st.session_state.log_container = None

def show_header():
    """Display main header"""
    
with open("images/icon.png", "rb") as f:
    data = f.read()
    data_base64 = base64.b64encode(data).decode()

st.markdown(
    f"""
    <div style='display: flex; align-items: center;'>
        <h1 style='margin-bottom: 0;'>Bulk RNAlyzer</h1>
        <img src='data:image/png;base64,{data_base64}' width='30' style='margin-bottom: 16px;'>
    </div>
    """,
    unsafe_allow_html=True
)

st.markdown("###### A complete workflow for analyzing and visualizing bulk RNA-seq data")
st.markdown("---")

def home_page():
    """Display home page with main options"""

    # Main options as radio buttons
    analysis_type = st.radio(
        "### Select an analysis type to get started:",
        options=["Index Genome", "Bulk RNA-Seq Analysis", "Data Visualizations"],
        index=1,
    )
    
    st.session_state.selected_analysis = analysis_type
    
    # Show description based on selection
    if analysis_type == "Index Genome":
        st.info("**Index Genome**: Create reference genome indexes for alignment tools")
    elif analysis_type == "Bulk RNA-Seq Analysis":
        st.info("**Bulk RNA-Seq Analysis**: Complete pipeline from raw reads to differential expression")
    else:
        st.info("**Data Visualizations**: Create interactive plots and visualizations from analysis results")
    
    # Proceed button
    if st.button("🚀 Proceed", type="primary", use_container_width=True):
        if st.session_state.selected_analysis == "Bulk RNA-Seq Analysis":
            st.session_state.current_page = "select_step"
        elif st.session_state.selected_analysis == "Index Genome":
            st.session_state.current_page = "index_genome"
        else:
            st.session_state.current_page = "visualization"
        st.rerun()

def select_step_page():
    """Page for selecting starting step in Bulk RNA-Seq analysis"""
    st.markdown("## 🧬 Bulk RNA-Seq Analysis")
    st.markdown("Select which step to start the analysis from:")
    
    # Step selection radio buttons
    start_step = st.radio(
        "### 📋 Select Starting Step",
        options=[1, 2, 3, 4, 5, 6],
        format_func=lambda x: {
            1: "1. 🔍 Quality Control (Before Trimming)",
            2: "2. ✂️ Trimming of Raw Reads", 
            3: "3. 🔍 Quality Control (After Trimming)",
            4: "4. 🧬 Alignment of Trimmed Reads to Reference Genome",
            5: "5. 📊 Quantification of Aligned Reads",
            6: "6. 📈 Differential Gene Expression Analysis"
        }[x],
        help="Choose the step to begin your analysis. Subsequent steps will run automatically."
    )
    
    st.session_state.start_step = start_step
    
    # Show step descriptions
    st.markdown("### 📝 Step Descriptions")
    step_descriptions = {
        1: "**Quality Control (Before Trimming)**: Check read quality with FastQC before processing",
        2: "**Read Trimming**: Remove adapters and low-quality bases with Trimmomatic",
        3: "**Quality Control (After Trimming)**: Check read quality after trimming with FastQC",
        4: "**Alignment**: Map reads to reference genome using HISAT2",
        5: "**Quantification**: Generate count matrices with StringTie",
        6: "**Differential Expression**: Identify differentially expressed genes with DESeq2"
    }
    
    st.info(step_descriptions[start_step])
    
    col1, col2 = st.columns([1, 4])
    
    with col1:
        if st.button("❮ Back", use_container_width=True):
            st.session_state.current_page = "home"
            st.rerun()
    
    with col2:
        if st.button("🔧 Configure Parameters", type="primary", use_container_width=True):
            st.session_state.current_page = "configure_params"
            st.rerun()

def index_genome_page():
    """Page for genome indexing configuration"""
    st.markdown("## 🧬 Genome Indexing")
    st.markdown("Build reference genome indexes for alignment tools")
    
    # Get indexing parameters
    get_index_genome_params()
    
    # Navigation buttons
    col1, col2 = st.columns([1, 4])
    
    with col1:
        if st.button("❮ Back", use_container_width=True):
            st.session_state.current_page = "home"
            st.rerun()
    
    with col2:
        if st.button("⚙️ Configure Parameters", type="primary", use_container_width=True):
            if validate_index_genome_params():
                st.session_state.current_page = "index_genome_advanced"
                st.rerun()
            else:
                st.error("❌ Please fix the validation errors above")

def index_genome_advanced_page():
    """Page for advanced genome indexing parameters"""
    st.markdown("## ⚙️ Genome Indexing - Parameters")
    
    # Get advanced parameters
    get_index_advanced_params()
    
    # Navigation buttons
    col1, col2, col3 = st.columns([1, 2, 1])
    
    with col1:
        if st.button("❮ Back", use_container_width=True):
            st.session_state.current_page = "index_genome"
            st.rerun()
    
    with col3:
        if st.button("🚀 Index Genome", type="primary", use_container_width=True):
            st.session_state.current_page = "running_index_genome"
            st.session_state.analysis_running = True
            st.rerun()

def running_index_genome_page():
    """Page showing real-time genome indexing progress"""
    st.markdown("## 🚀 Running Genome Indexing")
    st.markdown("---")
    
    # Run the indexing process
    if st.session_state.analysis_running:
        success = run_index_genome_step(st.session_state.params, Path(st.session_state.params['index_output_dir']))
        
        if success:
            st.success("✅ Genome indexing completed successfully!")
            st.balloons()
            st.session_state.analysis_running = False
        else:
            st.error("❌ Genome indexing failed!")
            st.session_state.analysis_running = False
    
    # Show completion options
    if not st.session_state.analysis_running:
        if st.button("Run a New Analysis", use_container_width=True):
            reset_session_state()
            st.rerun()

def configure_parameters_page():
    """Page for configuring ALL parameters for the entire pipeline
       — renders all inputs inside a single placeholder so we can
       remove them instantly when Start is clicked.
    """
    st.markdown(f"## ⚙️ Pipeline Configuration - Starting from Step {st.session_state.start_step}")

    # If analysis already started (safety), switch to running page immediately
    if st.session_state.get("analysis_running", False):
        st.session_state.current_page = "running_analysis"
        st.rerun()
        return

    # We'll render the entire configuration UI inside this placeholder
    config_placeholder = st.empty()

    # Use the placeholder's container to draw everything. We keep a reference
    # to start_clicked so we can remove the whole placeholder immediately.
    with config_placeholder.container():
        # Get base output directory
        base_output_dir = st.text_input(
            "📁 **Base Output Directory**",
            help="Directory where all analysis results will be saved",
            placeholder="/path/to/output_directory"
        )

        if base_output_dir:
            st.session_state.params['base_output_dir'] = base_output_dir
            output_dir = Path(base_output_dir)
            try:
                output_dir.mkdir(parents=True, exist_ok=True)
                st.session_state.output_dir = output_dir
                st.success(f"✅ Output directory: {output_dir}")
            except Exception as e:
                st.error(f"❌ Cannot create output directory: {e}")

        start_step = st.session_state.start_step

        # Sequencing type
        if start_step <= 4:
            mode = st.radio(
                "🔬 **Sequencing Type**",
                options=["PE", "SE"],
                format_func=lambda x: "Paired-End (PE)" if x == "PE" else "Single-End (SE)",
                horizontal=True
            )
            st.session_state.params['mode'] = mode

        # Step-by-step configuration (kept exactly as before)
        st.markdown("### 🔧 Step-by-Step Configuration")

        if start_step <= 1:
            st.markdown("---")
            st.markdown("#### 🔍 Step 1: Quality Control Parameters (Before Trimming)")
            get_qc_params("before_trim")

        if start_step <= 2:
            st.markdown("---")
            st.markdown("#### ✂️ Step 2: Trimming Parameters")
            get_trimming_params()

        if start_step <= 3:
            st.markdown("---")
            st.markdown("#### 🔍 Step 3: Quality Control Parameters (After Trimming)")
            get_qc_params("after_trim")

        if start_step <= 4:
            st.markdown("---")
            st.markdown("#### 🧬 Step 4: Alignment Parameters")
            get_alignment_params()

        if start_step <= 5:
            st.markdown("---")
            st.markdown("#### 📊 Step 5: Quantification Parameters")
            get_quantification_params()

        if start_step <= 6:
            st.markdown("---")
            st.markdown("#### 📈 Step 6: DESeq2 Analysis Parameters")
            get_deseq2_params()

        # Navigation buttons placed inside the same placeholder so they are removed together
        col1, col2, col3 = st.columns([1, 2, 1])

        with col1:
            if st.button("❮ Back", use_container_width=True):
                st.session_state.current_page = "select_step"
                st.rerun()

        # NOTE: we capture click result in a variable so we can act on it AFTER removing the placeholder.
        with col3:
            start_clicked = st.button("🚀 Start Analysis", type="primary", use_container_width=True)

    # --- AFTER the placeholder container scope ---
    # If the Start button was clicked, immediately remove the placeholder UI
    # so the configuration inputs do not remain dimmed. Then set the running flags and rerun.
    if start_clicked:
        # Validate parameters first (we reuse your validate function)
        if validate_all_parameters():
            # Immediately clear the configuration UI from the page (same render)
            config_placeholder.empty()

            # Prepare pipeline
            st.session_state.pipeline_status = {}
            st.session_state.analysis_running = True
            st.session_state.current_page = "running_analysis"

            # Force a rerun so running_analysis_page() renders right away
            st.rerun()
        else:
            # If validation failed, keep the UI (we re-render since we did not empty placeholder)
            st.error("❌ Please fix the validation errors above")
            # no rerun here; user sees validation errors in the same render


def visualization_page():
    """Page for data visualizations"""
    st.markdown("## 📊 Data Visualizations")
    st.markdown("Generate heatmaps and boxplots from your count data")
    
    # Visualization type selection
    viz_type = st.radio(
        "### 🔍 Choose Visualization Type",
        options=["Heatmap", "Box-plot"],
        help="Select the type of visualization to generate"
    )
    
    # Common inputs
    st.markdown("### 📁 Input Files")
    
    # Count matrix type
    count_type = st.radio(
        "Count Matrix Type",
        options=["Gene", "Transcript"],
        horizontal=True
    )
    
    count_file = st.file_uploader(
        f"Upload {count_type} Count Matrix (CSV/TSV)",
        type=["csv", "tsv", "txt"],
        help="Count matrix file with genes/transcripts as rows and samples as columns"
    )
    
    metadata_file = st.file_uploader(
        "Upload Metadata File (CSV/TSV)",
        type=["csv", "tsv", "txt"],
        help="Metadata file with sample information and conditions"
    )
    
    # Gene input - only comma-separated names now
    genes_input = st.text_input(
        "Enter Gene Names (comma-separated)",
        placeholder="Gene1, Gene2, Gene3, ...",
        help="Comma-separated list of gene names"
    )
    
    # Output settings
    st.markdown("### ⚙️ Output Settings")
    output_dir = st.text_input(
        "Output Directory",
        placeholder="/path/to/output/directory",
        help="Directory where plots will be saved"
    )
    
    col1, col2 = st.columns(2)
    with col1:
        width = st.number_input("Width (pixels)", min_value=100, max_value=10000, value=3000)
    with col2:
        height = st.number_input("Height (pixels)", min_value=100, max_value=10000, value=2600)
    
    # Heatmap-specific inputs - only show for heatmap
    color_palette = "Green_Teal_Purple"
    gene_annotation_file = None
    
    if viz_type == "Heatmap":
        st.markdown("### 🔥 Heatmap Settings")
        
        # Color palette selection
        color_palette = st.selectbox(
            "Color Palette",
            options=["Green_Teal_Purple", "blue_teal", "orange_blue_white", "red_blue_diverging", "warmcool"],
            help="Choose color scheme for the heatmap"
        )
        
        gene_annotation_file = st.file_uploader(
            "Upload Gene Annotation File (Optional)",
            type=["csv", "tsv", "txt"],
            help="Optional file with additional gene annotations"
        )
    
    # Generate button
    if st.button("🚀 Generate Visualizations", type="primary", use_container_width=True):
        # Validate inputs
        if not count_file:
            st.error("❌ Please upload a count matrix file")
            return
        if not metadata_file:
            st.error("❌ Please upload a metadata file")
            return
        if not genes_input:
            st.error("❌ Please enter gene names")
            return
        if not output_dir:
            st.error("❌ Please specify an output directory")
            return
        
        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save uploaded files
        count_path = output_path / count_file.name
        with open(count_path, "wb") as f:
            f.write(count_file.getbuffer())
        
        metadata_path = output_path / metadata_file.name
        with open(metadata_path, "wb") as f:
            f.write(metadata_file.getbuffer())
        
        # Handle gene input - only comma-separated now
        genes_path = output_path / "genes_list.txt"
        # Write genes with proper newlines
        genes_list = [gene.strip() for gene in genes_input.split(',')]
        with open(genes_path, "w") as f:
            f.write("\n".join(genes_list) + "\n")
        
        # Handle gene annotation file - only for heatmap
        gene_annotation_path = "None"
        if viz_type == "Heatmap" and gene_annotation_file is not None:
            gene_annotation_path = output_path / gene_annotation_file.name
            with open(gene_annotation_path, "wb") as f:
                f.write(gene_annotation_file.getbuffer())
            gene_annotation_path = str(gene_annotation_path)
        
        # Run visualizations
        plot_type = "gene" if count_type == "Gene" else "transcript"
        
        # Get modules directory
        modules_dir = Path(__file__).parent / "modules"
        
        try:
            if viz_type == "Heatmap":
                st.info("🔄 Generating heatmap...")
                
                # Build command arguments - DPI is now hardcoded to 300
                heatmap_args = [
                    str(count_path), str(metadata_path), str(output_path), plot_type,
                    str(genes_path), str(width), str(height), "300", color_palette
                ]
                if gene_annotation_path != "None":
                    heatmap_args.append(gene_annotation_path)
                
                # Show the command being run for debugging
                st.code(f"Running: Rscript {Config.HEATMAP_SCRIPT_PATH} {' '.join(heatmap_args)}")
                
                result = subprocess.run(
                    ["Rscript", Config.HEATMAP_SCRIPT_PATH] + heatmap_args,
                    capture_output=True, text=True, cwd=str(modules_dir)
                )
                
                # Display output for debugging
                if result.stdout:
                    st.text("Heatmap output:")
                    st.code(result.stdout)
                if result.stderr:
                    st.text("Heatmap errors:")
                    st.code(result.stderr)
                
                if result.returncode == 0:
                    st.success("✅ Heatmap generated successfully!")
                    # Display heatmap
                    heatmap_file = output_path / f"heatmap_{plot_type}.png"
                    if heatmap_file.exists():
                        st.image(str(heatmap_file), caption="Generated Heatmap", width=500)
                    else:
                        # Try alternative filename
                        heatmap_files = list(output_path.glob("heatmap*.png"))
                        if heatmap_files:
                            st.image(str(heatmap_files[0]), caption="Generated Heatmap", width=500)
                        else:
                            st.warning("⚠️ Heatmap file was not found in output directory")
                else:
                    st.error("❌ Heatmap generation failed!")
            
            elif viz_type == "Box-plot":
                st.info("🔄 Generating boxplots...")
                
                # Build command arguments - DPI is now hardcoded to 300
                boxplot_args = [
                    str(count_path), str(metadata_path), str(output_path), plot_type,
                    str(genes_path), str(width), str(height), "300"
                ]
                
                # Show the command being run for debugging
                st.code(f"Running: Rscript {Config.BOXPLOT_SCRIPT_PATH} {' '.join(boxplot_args)}")
                
                result = subprocess.run(
                    ["Rscript", Config.BOXPLOT_SCRIPT_PATH] + boxplot_args,
                    capture_output=True, text=True, cwd=str(modules_dir)
                )
                
                # Display output for debugging
                if result.stdout:
                    st.text("Boxplot output:")
                    st.code(result.stdout)
                if result.stderr:
                    st.text("Boxplot errors:")
                    st.code(result.stderr)
                
                if result.returncode == 0:
                    st.success("✅ Boxplots generated successfully!")
                    # Display boxplots
                    boxplot_files = list(output_path.glob("*_boxplot.png"))
                    if boxplot_files:
                        for boxplot_file in boxplot_files:
                            st.image(str(boxplot_file), caption=boxplot_file.name, width=400)
                    else:
                        st.warning("⚠️ No boxplot files were found in output directory")
                else:
                    st.error("❌ Boxplot generation failed!")
        
        except Exception as e:
            st.error(f"❌ Visualization generation failed: {str(e)}")
            st.info("💡 Check that all R packages are installed: ComplexHeatmap, circlize, ggplot2, dplyr, tidyr, paletteer")

def validate_all_parameters():
    """Validate all parameters for the entire pipeline"""
    if not st.session_state.params.get('base_output_dir'):
        st.error("❌ Please specify a base output directory")
        return False
    
    start_step = st.session_state.start_step
    
    # Validate each step's parameters
    validation_checks = {
        1: validate_qc_params,
        2: validate_trimming_params,
        3: validate_qc_after_trim_params,
        4: validate_alignment_params,
        5: validate_quantification_params,
        6: validate_deseq2_params
    }
    
    # Only validate steps from start_step onward
    for step in range(start_step, 7):
        if step in validation_checks:
            if not validation_checks[step]():
                st.error(f"❌ Validation failed for Step {step}")
                return False
    
    return True

def validate_qc_params():
    """Validate QC parameters"""
    if not st.session_state.params.get('input_dir'):
        st.error("❌ Please specify input directory with FASTQ files for QC")
        return False
    return True

def validate_trimming_params():
    """Validate trimming parameters"""
    if not st.session_state.params.get('input_dir'):
        st.error("❌ Please specify input directory with FASTQ files for trimming")
        return False
    if not st.session_state.params.get('adapter_file'):
        st.error("❌ Please specify adapter file for trimming")
        return False
    return True

def validate_qc_after_trim_params():
    """Validate QC after trim parameters"""
    # If starting from step 3, require manual input
    if st.session_state.start_step == 3:
        if not st.session_state.params.get('qc_after_trim_input_dir'):
            st.error("❌ Please specify input directory with trimmed FASTQ files for QC")
            return False
    return True

def validate_alignment_params():
    """Validate alignment parameters"""
    if not st.session_state.params.get('input_dir'):
        st.error("❌ Please specify input directory for alignment")
        return False
    if not st.session_state.params.get('index_base'):
        st.error("❌ Please specify HISAT2 index base path")
        return False
    return True

def validate_quantification_params():
    """Validate quantification parameters"""
    start_step = st.session_state.start_step
    
    # If starting from step 5, require BAM directory
    if start_step == 5:
        if not st.session_state.params.get('bam_dir'):
            st.error("❌ Please specify BAM files directory for quantification")
            return False
        # Validate that the BAM directory exists
        bam_dir = Path(st.session_state.params['bam_dir'])
        if not bam_dir.exists():
            st.error(f"❌ BAM directory does not exist: {bam_dir}")
            return False
    else:
        # For steps 1-4, BAM directory will be auto-detected from alignment output
        st.session_state.params['bam_dir'] = "auto_detect"
    
    # Always require annotation file
    if not st.session_state.params.get('annotation_file'):
        st.error("❌ Please specify annotation file for quantification")
        return False
    
    # Validate annotation file exists
    annotation_file = Path(st.session_state.params['annotation_file'])
    if not annotation_file.exists():
        st.error(f"❌ Annotation file not found: {annotation_file}")
        return False
    
    return True

def validate_deseq2_params():
    """Validate DESeq2 parameters"""
    if not st.session_state.params.get('metadata_file'):
        st.error("❌ Please specify metadata file for DESeq2")
        return False
    
    test_type = st.session_state.params.get('test', 'Wald')
    
    if test_type == 'Wald':
        if not st.session_state.params.get('contrasts'):
            st.error("❌ Please add at least one condition comparison for Wald test")
            return False
    else:  # LRT test
        if not st.session_state.params.get('lrt_comparisons'):
            st.error("❌ Please select at least one model combination for LRT test")
            return False
    
    return True

def running_analysis_page():
    """Page showing real-time analysis progress"""
    st.markdown("## 🚀 Running Analysis")
    st.markdown("---")
    
    # Initialize pipeline status if not exists
    if not st.session_state.pipeline_status:
        initialize_pipeline_status()

    # Show the progress section FIRST
    display_progress()
    
    # Initialize log content if not exists
    if 'log_content' not in st.session_state:
        st.session_state.log_content = "Starting analysis...\n"
    
    # Then show the log section below
    st.markdown("#### 📝 Analysis Log")
    if st.session_state.log_container is None:
        st.session_state.log_container = st.empty()
    log_container = st.session_state.log_container

    with log_container.container():
        st.text_area(
            "Log Output",
            value=st.session_state.log_content,
            height=400,
            key="analysis_log_display"
        )
    
    # Create a container for the log
    if st.session_state.log_container is None:
        st.session_state.log_container = st.empty()
    log_container = st.session_state.log_container
    
    # Run the pipeline if analysis is running
    if st.session_state.analysis_running:
        run_pipeline_with_progress()
        # Force a rerun to update the UI
        st.rerun()
    
    # Show completion message when done but KEEP the log visible
    if not st.session_state.analysis_running and st.session_state.pipeline_status:
        # Check if all steps are completed
        all_completed = all(
            step_info["status"] == "completed"
            for step_info in st.session_state.pipeline_status.values()
        )

        # Keep showing progress + log even after completion
        st.markdown("#### 📝 Analysis Log")
        with st.session_state.log_container.container():
            st.text_area(
                "Log Output",
                value=st.session_state.log_content,
                height=400,
                key="analysis_log_display_final"
            )

        if all_completed:
            st.balloons()
            st.success("🎉 Analysis completed successfully!")
        else:
            st.error("❌ Analysis completed with errors!")



        # Option to restart
        if st.button("Run a New Analysis", use_container_width=True):
            reset_session_state()
            st.rerun()

def initialize_pipeline_status():
    """Initialize the pipeline status tracking"""
    start_step = st.session_state.start_step
    
    # Define all steps with their dependencies
    all_steps = {
        1: {"name": "Quality Control (Before Trimming)", "status": "pending"},
        2: {"name": "Trimming of Raw Reads", "status": "pending"},
        3: {"name": "Quality Control (After Trimming)", "status": "pending"},
        4: {"name": "Alignment to Reference Genome", "status": "pending"},
        5: {"name": "Quantification", "status": "pending"},
        6: {"name": "Differential Gene Expression", "status": "pending"}
    }
    
    # Only show steps from start_step onward
    st.session_state.pipeline_status = {
        step: info for step, info in all_steps.items() 
        if step >= st.session_state.start_step
    }

def run_pipeline_with_progress():
    """Run the pipeline with real-time progress updates"""
    if not st.session_state.analysis_running:
        return

    start_step = st.session_state.start_step

    step_functions = {
        1: run_qc_before_trim_step,
        2: run_trimming_step,
        3: run_qc_after_trim_step,
        4: run_alignment_step,
        5: run_quantification_step,
        6: run_deseq2_step,
    }

    # Find the next step to run
    current_step = None
    for step_num in range(start_step, 7):
        step_info = st.session_state.pipeline_status[step_num]
        if step_info["status"] not in ["completed", "error"]:
            current_step = step_num
            break

    if current_step is None:
        st.session_state.analysis_running = False
        return

    step_info = st.session_state.pipeline_status[current_step]

    # --- 🟡 Mark step as running and refresh UI ---
    if step_info["status"] != "running":
        step_info["status"] = "running"
        st.session_state.pipeline_status[current_step] = step_info
        add_to_log(f"🔄 Starting Step {current_step}: {step_info['name']}")
        st.session_state._force_refresh = True  # dummy flag to trigger rerun
        st.rerun()  # immediately re-render UI showing "Running 🔄"
        return  # stop here so rerun happens before the step starts

    # --- 🧩 Actually run the step now ---
    success = run_step_with_logging(current_step, step_functions[current_step], step_info["name"])

    if success:
        step_info["status"] = "completed"
        add_to_log(f"✅ Step {current_step} completed successfully!")
    else:
        step_info["status"] = "error"
        add_to_log(f"❌ Pipeline failed at step {current_step}")
        st.session_state.analysis_running = False
        return

    # Continue automatically to next step if analysis still running
    # (this will show the next "Running 🔄" once rerendered)
    run_pipeline_with_progress()


def run_step_with_logging(step_num, step_function, step_name):
    """Run a single step with proper logging - FIXED"""
    # Add step header to log
    add_to_log(f"=== Step {step_num}: {step_name} ===")
    
    try:
        # Run the step function
        success = step_function(st.session_state.params, st.session_state.output_dir)
        
        if success:
            add_to_log(f"✅ {step_name} completed successfully!")
        else:
            add_to_log(f"❌ {step_name} failed!")
            
        return success
            
    except Exception as e:
        add_to_log(f"💥 Error in {step_name}: {str(e)}")
        return False

def display_progress():
    """Display a compact progress table as a clean bullet list"""
    st.markdown("#### 📊 Analysis Progress")

    progress_lines = []
    for step_num, step_info in st.session_state.pipeline_status.items():
        status_icon = {
            "pending": "⏳ In Queue",
            "running": "🔄 Running",
            "completed": "✅ Completed",
            "error": "❌ Error(s) Occurred",
        }.get(step_info["status"], "⏳ In Queue")
        
        # Use columns for compact display
        col1, col2 = st.columns([2, 4])
        with col1:
            st.write(f"● **{step_info['name']}**")
        with col2:
            st.write(status_icon)
    
    st.markdown("---")

def reset_session_state():
    """Reset session state for a new analysis"""
    st.session_state.current_page = "home"
    st.session_state.selected_analysis = None
    st.session_state.start_step = None
    st.session_state.params = {}
    st.session_state.pipeline_status = {}
    st.session_state.output_dir = None
    st.session_state.analysis_running = False
    
    if 'log_content' in st.session_state:
        del st.session_state.log_content
    if 'log_container' in st.session_state:
        del st.session_state.log_container

def main():
    """Main Streamlit application"""
    # Page configuration
    st.set_page_config(
        page_title="Bulk RNAlyzer",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    # Initialize session state
    init_session_state()
    
    # Show sidebar with tool status
    show_sidebar()
    
    # Show header
    show_header()
    
    # Route to appropriate page
    if st.session_state.current_page == "home":
        home_page()
    elif st.session_state.current_page == "select_step":
        select_step_page()
    elif st.session_state.current_page == "configure_params":
        configure_parameters_page()
    elif st.session_state.current_page == "running_analysis":
        running_analysis_page()
    elif st.session_state.current_page == "index_genome":
        index_genome_page()
    elif st.session_state.current_page == "index_genome_advanced":
        index_genome_advanced_page()
    elif st.session_state.current_page == "running_index_genome":
        running_index_genome_page()
    elif st.session_state.current_page == "visualization":
        visualization_page()

def show_sidebar():
    """Display sidebar with tool information"""
    st.sidebar.title("Bulk RNAlyzer")
    st.sidebar.markdown("---")
    
    st.sidebar.markdown("### 🛠️ Tool Status")
    
    tools = {
        "hisat2-build": "hisat2-build",
        "FastQC": "fastqc",
        "Trimmomatic": "java", 
        "HISAT2": "hisat2",
        "samtools": "samtools",
        "sambamba": "sambamba",
        "StringTie": "stringtie",
        "gffread": "gffread",
        "R": "Rscript"
    }
    
    for tool, cmd in tools.items():
        status = "✅" if check_tool_available(cmd) else "❌"
        st.sidebar.markdown(f"{status} {tool}")
    
    st.sidebar.markdown("---")
    st.sidebar.markdown("### 📋 Current Session")
    if st.session_state.current_page != "home":
        st.sidebar.info(f"**Analysis**: {st.session_state.selected_analysis}")
        if st.session_state.start_step:
            st.sidebar.info(f"**Starting Step**: {st.session_state.start_step}")

if __name__ == "__main__":
    main()