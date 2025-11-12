# modules/index_genome_module.py
import streamlit as st
import subprocess
from pathlib import Path

def get_index_genome_params():
    """Get parameters for genome indexing"""
    
    # Index type selection
    index_type = st.radio(
        "## 🔧 Select Index Type",
        options=[
            "HFM_index",
            "HGFM_index_transcripts", 
            "HGFM_index_SNPs",
            "HGFM_index_SNPs_transcripts"
        ],
        format_func=lambda x: {
            "HFM_index": "Build HFM index (basic)",
            "HGFM_index_transcripts": "Build HGFM index with transcripts",
            "HGFM_index_SNPs": "Build HGFM index with SNPs", 
            "HGFM_index_SNPs_transcripts": "Build HGFM index with SNPs and transcripts"
        }[x],
        help="Select the type of genome index to build"
    )
    
    st.session_state.params['index_type'] = index_type
    
    # Common parameters
    st.markdown("---")
    st.markdown("### 📁 Input Files")
    
    # Genome FASTA file (required for all)
    genome_fasta = st.text_input(
        "**Genome FASTA File Path** *",
        placeholder="/path/to/genome.fa",
        help="Path to the genome FASTA file"
    )
    if genome_fasta:
        st.session_state.params['genome_fasta'] = genome_fasta
    
    # Output directory
    output_dir = st.text_input(
        "📁 **Output Directory** *",
        placeholder="/path/to/index_output",
        help="Directory where index files will be saved"
    )
    if output_dir:
        st.session_state.params['index_output_dir'] = output_dir
    
    # Index basename
    default_basenames = {
        "HFM_index": "genome",
        "HGFM_index_transcripts": "genome_tran", 
        "HGFM_index_SNPs": "genome_snp",
        "HGFM_index_SNPs_transcripts": "genome_snp_tran"
    }
    
    index_basename = st.text_input(
        "**Index Basename**",
        value=default_basenames[index_type],
        help="Base name for the index files"
    )
    st.session_state.params['index_basename'] = index_basename
    
    # GTF file required for transcript-based indexing
    if index_type in ["HGFM_index_transcripts", "HGFM_index_SNPs_transcripts"]:
        st.markdown("---")
        st.markdown("#### Transcript Files")
        
        # File type selection
        annotation_type = st.radio(
            "**Annotation File Type**",
            options=["GTF", "GFF3"],
            help="Select whether you have GTF or GFF3 annotation file",
            horizontal=True
        )
        st.session_state.params['annotation_type'] = annotation_type
        
        if annotation_type == "GTF":
            annotation_file = st.text_input(
                "**GTF Annotation File** *",
                placeholder="/path/to/annotation.gtf",
                help="GTF file for extracting splice sites and exons"
            )
            if annotation_file:
                st.session_state.params['annotation_file'] = annotation_file
        else:
            annotation_file = st.text_input(
                "**GFF3 Annotation File** *",
                placeholder="/path/to/annotation.gff3",
                help="GFF3 file that will be converted to GTF for extracting splice sites and exons"
            )
            if annotation_file:
                st.session_state.params['annotation_file'] = annotation_file
            
            # Optional: output path for converted GTF
            converted_gtf_path = st.text_input(
                "**Converted GTF Output Path** (optional)",
                placeholder="/path/to/converted_annotation.gtf",
                help="Path where the converted GTF file will be saved. If not provided, will be saved in output directory."
            )
            if converted_gtf_path:
                st.session_state.params['converted_gtf_path'] = converted_gtf_path
        
        # Option to use existing exon/ss files or generate new
        use_existing = st.radio(
            "**Splice site and exon files**",
            options=["generate", "existing"],
            format_func=lambda x: "Generate from annotation" if x == "generate" else "Use existing files",
            horizontal=True
        )
        
        if use_existing == "existing":
            exon_file = st.text_input(
                "**Exon File Path**",
                placeholder="/path/to/genome.exon",
                help="Path to existing exon file"
            )
            if exon_file:
                st.session_state.params['exon_file'] = exon_file
            
            ss_file = st.text_input(
                "**Splice Sites File Path**", 
                placeholder="/path/to/genome.ss",
                help="Path to existing splice sites file"
            )
            if ss_file:
                st.session_state.params['ss_file'] = ss_file
        else:
            st.session_state.params['exon_file'] = "generate"
            st.session_state.params['ss_file'] = "generate"
    
    # SNP files required for SNP-based indexing  
    if index_type in ["HGFM_index_SNPs", "HGFM_index_SNPs_transcripts"]:
        st.markdown("---")
        st.markdown("#### SNP Files")
        
        snp_option = st.radio(
            "**SNP Data Source**",
            options=["existing", "generate"],
            format_func=lambda x: "Use existing SNP files" if x == "existing" else "Generate from VCF",
            horizontal=True
        )
        
        if snp_option == "existing":
            snp_file = st.text_input(
                "**SNP File Path** *",
                placeholder="/path/to/genome.snp",
                help="Path to existing SNP file"
            )
            if snp_file:
                st.session_state.params['snp_file'] = snp_file
            
            haplotype_file = st.text_input(
                "**Haplotype File Path** *",
                placeholder="/path/to/genome.haplotype", 
                help="Path to existing haplotype file"
            )
            if haplotype_file:
                st.session_state.params['haplotype_file'] = haplotype_file
        else:
            # Generate SNP files from VCF
            vcf_file = st.text_input(
                "**VCF File Path** *",
                placeholder="/path/to/variants.vcf",
                help="VCF file containing variants for SNP extraction"
            )
            if vcf_file:
                st.session_state.params['vcf_file'] = vcf_file
            
            snp_output_dir = st.text_input(
                "**SNP Files Output Directory**",
                value=output_dir if output_dir else "",
                placeholder="/path/to/snp_output",
                help="Directory where SNP and haplotype files will be saved"
            )
            if snp_output_dir:
                st.session_state.params['snp_output_dir'] = snp_output_dir

def get_index_advanced_params():
    """Get parameters for genome indexing"""
    st.markdown("---")
    st.markdown("### ⚙️ Parameters")
    
    # Threads
    threads = st.number_input(
        "🧵 **Number of Threads**",
        min_value=2,
        max_value=200,
        value=16,
        step=2,
        help="Number of CPU threads to use for indexing"
    )
    st.session_state.params['index_threads'] = threads
    
    # Large index option
    large_index = st.checkbox(
        "**Use Large Index**",
        value=False,
        help="Use if your genome is larger than 4 billion bases (e.g., human + extra scaffolds, large plants, or combined genome-transcriptome builds)"
    )
    st.session_state.params['large_index'] = large_index
    
    # Offrate
    offrate = st.number_input(
        "**Offrate Value**",
        min_value=1,
        max_value=100,
        value=4,
        help="Controls how dense the suffix array index is; Lower values → denser index → faster alignment but larger index size. Higher values → sparser index → smaller index but slower alignment."
    )
    st.session_state.params['offrate'] = offrate

def validate_index_genome_params():
    """Validate genome indexing parameters"""
    if not st.session_state.params.get('genome_fasta'):
        st.error("❌ Please specify genome FASTA file path")
        return False
    
    if not st.session_state.params.get('index_output_dir'):
        st.error("❌ Please specify output directory")
        return False
    
    index_type = st.session_state.params.get('index_type')
    
    # Validate transcript-based indexing
    if index_type in ["HGFM_index_transcripts", "HGFM_index_SNPs_transcripts"]:
        if not st.session_state.params.get('annotation_file'):
            st.error("❌ Please specify annotation file for transcript-based indexing")
            return False
        
        # If using existing files, check they are provided
        if (st.session_state.params.get('exon_file') != "generate" and 
            not st.session_state.params.get('exon_file')):
            st.error("❌ Please specify exon file path or select 'Generate from annotation'")
            return False
            
        if (st.session_state.params.get('ss_file') != "generate" and 
            not st.session_state.params.get('ss_file')):
            st.error("❌ Please specify splice sites file path or select 'Generate from annotation'")
            return False
    
    # Validate SNP-based indexing
    if index_type in ["HGFM_index_SNPs", "HGFM_index_SNPs_transcripts"]:
        # If using existing SNP files
        if st.session_state.params.get('snp_file') and st.session_state.params.get('haplotype_file'):
            if not st.session_state.params.get('snp_file'):
                st.error("❌ Please specify SNP file path")
                return False
            if not st.session_state.params.get('haplotype_file'):
                st.error("❌ Please specify haplotype file path")
                return False
        # If generating from VCF
        elif st.session_state.params.get('vcf_file'):
            if not st.session_state.params.get('vcf_file'):
                st.error("❌ Please specify VCF file path for SNP generation")
                return False
        else:
            st.error("❌ Please either provide existing SNP files or VCF file to generate them")
            return False
    
    return True

def run_index_genome_step(params, output_dir):
    """Run genome indexing step"""
    try:
        # Create output directory
        index_dir = Path(params['index_output_dir'])
        index_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate transcript files if needed
        if params['index_type'] in ["HGFM_index_transcripts", "HGFM_index_SNPs_transcripts"]:
            if params.get('exon_file') == "generate" or params.get('ss_file') == "generate":
                if not generate_transcript_files(params, index_dir):
                    return False
        
        # Generate SNP files if needed
        if params['index_type'] in ["HGFM_index_SNPs", "HGFM_index_SNPs_transcripts"]:
            if params.get('vcf_file'):
                if not generate_snp_files(params, index_dir):
                    return False
        
        # Build the index
        return build_hisat2_index(params, index_dir)
        
    except Exception as e:
        st.error(f"❌ Error in genome indexing: {str(e)}")
        return False

def convert_gff3_to_gtf(params, output_dir):
    """Convert GFF3 file to GTF format using gffread"""
    try:
        gff3_file = params['annotation_file']
        
        # Determine output GTF path
        if params.get('converted_gtf_path'):
            gtf_file = Path(params['converted_gtf_path'])
        else:
            gtf_file = output_dir / f"{Path(gff3_file).stem}.converted.gtf"
        
        # Create command for gffread conversion
        cmd = f"gffread {gff3_file} -T -o {gtf_file}"
        
        st.info(f"Converting GFF3 to GTF: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            st.error(f"❌ Failed to convert GFF3 to GTF: {result.stderr}")
            return None
        
        st.success(f"✅ Successfully converted GFF3 to GTF: {gtf_file}")
        return str(gtf_file)
        
    except Exception as e:
        st.error(f"❌ Error converting GFF3 to GTF: {str(e)}")
        return None

def generate_transcript_files(params, output_dir):
    """Generate exon and splice site files from annotation file"""
    try:
        annotation_file = params['annotation_file']
        annotation_type = params.get('annotation_type', 'GTF')
        
        # If GFF3 file is provided, convert it to GTF first
        if annotation_type == "GFF3":
            st.info("🔧 GFF3 file detected, converting to GTF format...")
            gtf_file = convert_gff3_to_gtf(params, output_dir)
            if not gtf_file:
                return False
            # Update the annotation file to use the converted GTF
            annotation_file = gtf_file
        else:
            gtf_file = annotation_file
        
        base_name = Path(gtf_file).stem
        
        # Generate splice sites
        ss_file = output_dir / f"{base_name}.ss"
        cmd_ss = f"hisat2_extract_splice_sites.py {gtf_file} > {ss_file}"
        
        st.info(f"🔧 Generating splice sites file: {cmd_ss}")
        result_ss = subprocess.run(cmd_ss, shell=True, capture_output=True, text=True)
        
        if result_ss.returncode != 0:
            st.error(f"❌ Failed to generate splice sites: {result_ss.stderr}")
            return False
        
        # Generate exons
        exon_file = output_dir / f"{base_name}.exon"
        cmd_exon = f"hisat2_extract_exons.py {gtf_file} > {exon_file}"
        
        st.info(f"🔧 Generating exon file: {cmd_exon}")
        result_exon = subprocess.run(cmd_exon, shell=True, capture_output=True, text=True)
        
        if result_exon.returncode != 0:
            st.error(f"❌ Failed to generate exon file: {result_exon.stderr}")
            return False
        
        # Update params with generated file paths
        params['ss_file'] = str(ss_file)
        params['exon_file'] = str(exon_file)
        
        st.success("✅ Successfully generated transcript files")
        return True
        
    except Exception as e:
        st.error(f"❌ Error generating transcript files: {str(e)}")
        return False

def generate_snp_files(params, output_dir):
    """Generate SNP and haplotype files from VCF"""
    try:
        genome_fasta = params['genome_fasta']
        vcf_file = params['vcf_file']
        base_name = Path(genome_fasta).stem
        
        snp_file = output_dir / f"{base_name}.snp"
        haplotype_file = output_dir / f"{base_name}.haplotype"
        
        cmd = f"hisat2_extract_snps_haplotypes_VCF.py {genome_fasta} {vcf_file} {snp_file} {haplotype_file}"
        
        st.info(f"🔧 Generating SNP and haplotype files: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            st.error(f"❌ Failed to generate SNP files: {result.stderr}")
            return False
        
        # Update params with generated file paths
        params['snp_file'] = str(snp_file)
        params['haplotype_file'] = str(haplotype_file)
        
        st.success("✅ Successfully generated SNP and haplotype files")
        return True
        
    except Exception as e:
        st.error(f"❌ Error generating SNP files: {str(e)}")
        return False

def build_hisat2_index(params, output_dir):
    """Build HISAT2 index with specified parameters"""
    try:
        # Base command
        cmd = f"hisat2-build -p {params.get('index_threads', 16)}"
        
        # Add large index option if selected
        if params.get('large_index'):
            cmd += " --large-index"
        
        # Add offrate
        cmd += f" --offrate {params.get('offrate', 4)}"
        
        # Add transcript options for transcript-based indexing
        if params['index_type'] in ["HGFM_index_transcripts", "HGFM_index_SNPs_transcripts"]:
            cmd += f" --exon {params['exon_file']} --ss {params['ss_file']}"
        
        # Add SNP options for SNP-based indexing
        if params['index_type'] in ["HGFM_index_SNPs", "HGFM_index_SNPs_transcripts"]:
            cmd += f" --snp {params['snp_file']} --haplotype {params['haplotype_file']}"
        
        # Add input and output
        cmd += f" {params['genome_fasta']} {output_dir / params['index_basename']}"
        
        st.info(f"🔧 Building HISAT2 index: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            st.success(f"✅ Successfully built genome index: {params['index_basename']}")
            st.info(f"📁 Index files saved in: {output_dir}")
            return True
        else:
            st.error(f"❌ Failed to build index: {result.stderr}")
            return False
            
    except Exception as e:
        st.error(f"❌ Error building HISAT2 index: {str(e)}")
        return False