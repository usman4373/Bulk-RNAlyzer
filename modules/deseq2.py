# modules/deseq2_module.py
import streamlit as st
import json
from pathlib import Path
from .config import Config
from .logger import add_to_log
import subprocess
import pandas as pd

def generate_model_combinations(columns):
    """Generate common model combinations for LRT testing based on available columns"""
    combinations = []
    
    # Remove common non-factor columns
    non_factor_cols = ['sample', 'sample_id', 'id', 'name', 'filename']
    factor_columns = [col for col in columns if col.lower() not in non_factor_cols]
    
    if not factor_columns:
        return combinations
    
    # Basic combinations
    for col in factor_columns:
        # Simple: ~col vs ~1
        combinations.append({
            'name': f"test_{col}_effect",
            'full': f"~ {col}",
            'reduced': "~ 1",
            'description': f"Test effect of {col}"
        })
    
    # Two-factor combinations
    if len(factor_columns) >= 2:
        for i, col1 in enumerate(factor_columns):
            for col2 in factor_columns[i+1:]:
                # Control for batch effect: ~col1 + col2 vs ~col1
                combinations.append({
                    'name': f"test_{col2}_controlling_{col1}",
                    'full': f"~ {col1} + {col2}",
                    'reduced': f"~ {col1}",
                    'description': f"Test {col2} effect controlling for {col1}"
                })
                
                # Control for batch effect: ~col1 + col2 vs ~col2
                combinations.append({
                    'name': f"test_{col1}_controlling_{col2}",
                    'full': f"~ {col1} + {col2}",
                    'reduced': f"~ {col2}",
                    'description': f"Test {col1} effect controlling for {col2}"
                })
                
                # Interaction effects (if we have at least 2 samples per combination)
                combinations.append({
                    'name': f"test_{col1}_{col2}_interaction",
                    'full': f"~ {col1} + {col2} + {col1}:{col2}",
                    'reduced': f"~ {col1} + {col2}",
                    'description': f"Test interaction between {col1} and {col2}"
                })
    
    # Three-factor combinations (more complex)
    if len(factor_columns) >= 3:
        for i, col1 in enumerate(factor_columns):
            for j, col2 in enumerate(factor_columns[i+1:]):
                for col3 in factor_columns[i+j+2:]:
                    combinations.append({
                        'name': f"test_{col3}_controlling_{col1}_{col2}",
                        'full': f"~ {col1} + {col2} + {col3}",
                        'reduced': f"~ {col1} + {col2}",
                        'description': f"Test {col3} effect controlling for {col1} and {col2}"
                    })
    
    # Set the simplest combination as default
    simplest = None
    for combo in combinations:
        if combo['reduced'] == "~ 1":
            simplest = combo
            break
    
    if simplest and simplest in combinations:
        combinations.remove(simplest)
        combinations.insert(0, simplest)
    
    return combinations

# modules/deseq2.py
def get_deseq2_params():
    """Collect DESeq2 parameters - ONLY ask for count matrices if starting step"""
    st.markdown("#### 📁 Input Parameters")
    
    # ALWAYS ask for metadata file (required for DESeq2)
    st.session_state.params['metadata_file'] = st.text_input(
        "**Metadata File (CSV)**",
        help="CSV file with sample metadata and conditions",
        placeholder="/path/to/metadata.csv"
    )
    
    # Only ask for count matrices if this is the starting step (Step 6)
    if st.session_state.start_step == 6:
        st.markdown("#### 📊 Count Matrix Files")
        
        col1, col2 = st.columns(2)
        with col1:
            st.session_state.params['gene_matrix_file'] = st.text_input(
                "**Gene Count Matrix Path**",
                help="Path to gene_count_matrix.csv",
                placeholder="/path/to/gene_count_matrix.csv"
            )
        with col2:
            st.session_state.params['transcript_matrix_file'] = st.text_input(
                "**Transcript Count Matrix Path**", 
                help="Path to transcript_count_matrix.csv",
                placeholder="/path/to/transcript_count_matrix.csv"
            )
    else:
        # For downstream steps, auto-detect from quantification output
        st.info("🔄 Count matrices will be auto-detected from quantification output")
    
    # Rest of the parameters (always ask for these)
    st.markdown("#### ⚙️ DESeq2 Parameters")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.session_state.params['lfc_value'] = st.number_input(
            "**Log Fold Change Threshold**", 
            value=1.0, 
            step=0.5,
            min_value=0.0,
            help="Minimum log fold change for significance"
        )
        
        st.session_state.params['test'] = st.selectbox(
            "**Statistical Test**", 
            ["Wald", "LRT"],
            help="DESeq2 statistical test method"
        )

        st.session_state.params['fitType'] = st.selectbox(
            "**Fit Type**",
            ["parametric", "local", "mean"],
            index=0,
            help="DESeq2 dispersion fit type"
        )
    
    with col2:
        st.session_state.params['altHypothesis'] = st.selectbox(
            "**Alternative Hypothesis**", 
            ["greaterAbs", "less", "greater"],
            help="Alternative hypothesis for testing"
        )
        
        st.session_state.params['alpha'] = st.number_input(
            "**Significance Level (alpha)**",
            value=0.05,
            min_value=0.001,
            max_value=0.1,
            step=0.005,
            help="Adjusted p-value cutoff for significance"
        )

        st.session_state.params['cooksCutoff'] = st.checkbox(
            "**Use Cook's Cutoff**",
            value=True,
            help="Use Cook's distance for outlier detection"
        )
    
    st.session_state.params['analysis_choice'] = st.radio(
        "**Analysis Type**",
        options=["1", "2", "3"],
        format_func=lambda x: {
            "1": "Gene Count Matrix Only",
            "2": "Transcript Count Matrix Only", 
            "3": "Both Gene and Transcript"
        }[x],
        horizontal=True
    )
    
    st.markdown("#### Condition Comparisons")
    
    # Initialize contrasts if not exists
    if 'contrasts' not in st.session_state.params:
        st.session_state.params['contrasts'] = []
    
    # Initialize LRT comparisons if not exists
    if 'lrt_comparisons' not in st.session_state.params:
        st.session_state.params['lrt_comparisons'] = []
    
    # Show different interfaces based on test type
    if st.session_state.params.get('test') == 'LRT':
        st.markdown("##### 📊 LRT Model Combinations")
        
        # Read metadata to get available columns for model combinations
        available_combinations = []
        if st.session_state.params.get('metadata_file') and Path(st.session_state.params['metadata_file']).exists():
            try:
                import pandas as pd
                metadata = pd.read_csv(st.session_state.params['metadata_file'])
                available_columns = metadata.columns.tolist()
                
                # Generate common model combinations based on available columns
                available_combinations = generate_model_combinations(available_columns)
                
                # Display dropdown for model combinations
                selected_combinations = st.multiselect(
                    "**Select full and reduced models**",
                    options=available_combinations,
                    default=[available_combinations[0]] if available_combinations else [],
                    format_func=lambda x: f"{x['full']} vs {x['reduced']}",
                    help="These options are given based on your metadata file. Select one or more model combinations for LRT testing."
                )
                
                # Update session state with selected combinations
                st.session_state.params['lrt_comparisons'] = selected_combinations
                
                if not selected_combinations:
                    st.warning("⚠️ Please select at least one model combination for LRT analysis")
                
            except Exception as e:
                st.error(f"❌ Error reading metadata file: {e}")
        else:
            st.warning("⚠️ Please provide a valid metadata file to generate model combinations")
    
    else:  # Wald test
        
        # Interface for adding contrasts
        col1, col2, col3 = st.columns([2, 2, 1])
        
        with col1:
            new_group1 = st.text_input("Condition 1", key="new_group1", placeholder="Treatment")
        with col2:
            new_group2 = st.text_input("Condition 2", key="new_group2", placeholder="Control")
        with col3:
            if st.button("➕ Add", use_container_width=True):
                if new_group1 and new_group2:
                    contrast_name = f"{new_group1}_vs_{new_group2}"
                    new_contrast = {
                        'name': contrast_name,
                        'group1': new_group1,
                        'group2': new_group2
                    }
                    if new_contrast not in st.session_state.params['contrasts']:
                        st.session_state.params['contrasts'].append(new_contrast)
                        st.success(f"✅ Added: {contrast_name}")
                    else:
                        st.warning("⚠️ Comparison already exists")
                else:
                    st.warning("⚠️ Please enter both conditions")
        
        # Display and manage existing contrasts
        if st.session_state.params['contrasts']:
            st.markdown("**📋 Current Comparisons:**")
            contrasts_to_remove = []
            
            for i, contrast in enumerate(st.session_state.params['contrasts']):
                col1, col2, col3, col4 = st.columns([3, 1, 3, 1])
                with col1:
                    st.write(f"**{contrast['group1']}**")
                with col2:
                    st.write("vs")
                with col3:
                    st.write(f"**{contrast['group2']}**")
                with col4:
                    if st.button("🗑️", key=f"remove_{i}"):
                        contrasts_to_remove.append(i)
            
            # Remove selected contrasts
            for i in sorted(contrasts_to_remove, reverse=True):
                removed = st.session_state.params['contrasts'].pop(i)
                st.info(f"🗑️ Removed: {removed['name']}")
                st.rerun()

def run_deseq2_step(params, output_dir):
    """Run DESeq2 differential expression analysis - WITH PROPER AUTO-DETECTION"""
    try:
        deseq2_dir = output_dir / "deseq2_analysis"
        deseq2_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if R script exists
        if not Path(Config.R_SCRIPT_PATH).exists():
            add_to_log(f"❌ R script not found: {Config.R_SCRIPT_PATH}")
            return False
        
        # Validate inputs
        if not params.get('metadata_file') or not Path(params['metadata_file']).exists():
            add_to_log("❌ Metadata file not found")
            return False
        
        if not params.get('contrasts'):
            add_to_log("❌ No condition comparisons defined")
            return False
        
        add_to_log("📈 Preparing DESeq2 analysis...")
        
        # Prepare parameters for R script
        analysis_params = {
            'mode': 'deseq2',
            'output_dir': str(deseq2_dir),
            'deseq2_params': {
                'metadata_file': params['metadata_file'],
                'lfc_value': str(params.get('lfc_value', 1.0)),
                'test': params.get('test', 'Wald'),
                'altHypothesis': params.get('altHypothesis', 'greaterAbs'),
                'analysis_choice': params.get('analysis_choice', '1'),
                'alpha': params.get('alpha', 0.05),
                'fitType': params.get('fitType', 'parametric'),
                'cooksCutoff': params.get('cooksCutoff', True),
                'contrasts': params.get('contrasts', []),
                'lrt_comparisons': params.get('lrt_comparisons', [])
            }
        }
        
        # Handle count matrix paths - AUTO-DETECT from quantification when not starting step
        if st.session_state.start_step == 6:
            # Starting from DESeq2 step, use provided count matrix paths
            if params.get('gene_matrix_file') and params['gene_matrix_file'] != "auto_detect":
                gene_matrix_path = Path(params['gene_matrix_file'])
                if gene_matrix_path.exists():
                    analysis_params['deseq2_params']['gene_matrix_file'] = str(gene_matrix_path)
                    add_to_log(f"✅ Using provided gene count matrix: {gene_matrix_path}")
            
            if params.get('transcript_matrix_file') and params['transcript_matrix_file'] != "auto_detect":
                transcript_matrix_path = Path(params['transcript_matrix_file'])
                if transcript_matrix_path.exists():
                    analysis_params['deseq2_params']['transcript_matrix_file'] = str(transcript_matrix_path)
                    add_to_log(f"✅ Using provided transcript count matrix: {transcript_matrix_path}")
        else:
            # Auto-detect from quantification output
            count_dir = output_dir / "quantification" / "count_matrix"
            if count_dir.exists():
                gene_matrix = count_dir / "gene_count_matrix.csv"
                transcript_matrix = count_dir / "transcript_count_matrix.csv"
                
                if gene_matrix.exists():
                    analysis_params['deseq2_params']['gene_matrix_file'] = str(gene_matrix)
                    add_to_log(f"✅ Auto-detected gene count matrix: {gene_matrix}")
                
                if transcript_matrix.exists():
                    analysis_params['deseq2_params']['transcript_matrix_file'] = str(transcript_matrix)
                    add_to_log(f"✅ Auto-detected transcript count matrix: {transcript_matrix}")
        
        # Validate that we have at least one count matrix
        has_gene = 'gene_matrix_file' in analysis_params['deseq2_params'] 
        has_transcript = 'transcript_matrix_file' in analysis_params['deseq2_params']
        
        if not has_gene and not has_transcript:
            add_to_log("❌ No count matrix files found! Please provide at least one count matrix.")
            return False
        
        # Save parameters to JSON file
        param_file = deseq2_dir / "deseq2_params.json"
        with open(param_file, 'w') as f:
            json.dump(analysis_params, f, indent=2)
        
        add_to_log(f"✅ Parameters saved to: {param_file}")
        
        # Run R script
        cmd = f"Rscript {Config.R_SCRIPT_PATH} {param_file}"
        
        add_to_log(f"Running: {cmd}")
        
        # Run command
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            add_to_log("✅ DESeq2 analysis completed successfully")
            # Check for DESeq2 output
            deseq2_outputs = list(deseq2_dir.glob("*"))
            if deseq2_outputs:
                add_to_log(f"✅ Generated {len(deseq2_outputs)} output files")
                add_to_log(f"📁 Results saved to: {deseq2_dir}")
                return True
            else:
                add_to_log("❌ No DESeq2 output files generated")
                return False
        else:
            add_to_log(f"❌ DESeq2 analysis failed: {result.stderr}")
            return False
        
    except Exception as e:
        add_to_log(f"❌ DESeq2 analysis failed: {str(e)}")
        return False
    