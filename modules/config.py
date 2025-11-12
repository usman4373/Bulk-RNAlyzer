import os

class Config:
    """Configuration settings for the pipeline"""

    # Get the absolute path of the directory where config.py is located
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))

    # Tool paths
    TRIMMOMATIC_JAR = os.path.join(BASE_DIR, '..', 'Trimmomatic-0.39', 'trimmomatic-0.39.jar')
    PREPDE_SCRIPT = os.path.join(BASE_DIR, 'prepDE.py3')
    R_SCRIPT_PATH = os.path.join(BASE_DIR, 'deseq2.R')

    # Visualization script paths
    HEATMAP_SCRIPT_PATH = os.path.join(BASE_DIR, 'heatmap.R')
    BOXPLOT_SCRIPT_PATH = os.path.join(BASE_DIR, 'boxplot.R')
    COMMON_SCRIPT_PATH = os.path.join(BASE_DIR, 'viz_file_handling.R')
