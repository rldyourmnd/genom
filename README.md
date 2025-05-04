# Transcription Factor Binding Site Analysis

This project analyzes gene promoter regions for transcription factor (TF) binding sites using motifs from the JASPAR database. It creates a comprehensive database of TFs and analyzes gene sequences for potential binding sites.

## Features

- Retrieves the latest transcription factor motifs from JASPAR database via API
- Creates a local SQLite database with detailed TF information
- Analyzes gene sequences for TF binding sites in the -2000 to +500bp region relative to TSS
- Generates visualizations of binding site distributions and patterns
- Creates an interactive HTML report with detailed gene and TF profiles

## Requirements

- Python 3.8+
- Required packages:
  - pyjaspar
  - pandas
  - biopython
  - matplotlib
  - seaborn
  - sqlalchemy

## Installation

1. Clone this repository:
```
git clone https://github.com/yourusername/tf-binding-analysis.git
cd tf-binding-analysis
```

2. Install required packages:
```
pip install pyjaspar pandas biopython matplotlib seaborn sqlalchemy
```

## Usage

### Basic Usage

Run the main script to perform a complete analysis:

```
python main.py
```

This will:
1. Fetch all TF motifs from JASPAR (if not already downloaded)
2. Create a SQLite database with TF information
3. Analyze all gene sequences in the `genes/` directory
4. Generate visualizations and summary reports
5. Create an HTML report in the `html_report/` directory

### Command Line Options

- `--genes-dir`: Directory containing gene FASTA files (default: 'genes')
- `--results-dir`: Directory to store results (default: 'results')
- `--html-dir`: Directory to store HTML visualization (default: 'html_report')
- `--threshold`: Threshold for binding site detection (0.0-1.0) (default: 0.8)
- `--test-mode`: Run in test mode with a subset of TFs
- `--update-db`: Force update of the JASPAR database
- `--skip-analysis`: Skip gene analysis, only generate visualization from existing results

### Examples

Run in test mode (faster, uses only a subset of TFs):
```
python main.py --test-mode
```

Update the JASPAR database and rerun the analysis:
```
python main.py --update-db
```

Use a different binding site detection threshold:
```
python main.py --threshold 0.7
```

Skip analysis and just generate HTML report from existing results:
```
python main.py --skip-analysis
```

## Project Structure

- `jaspar_api.py`: Interface to JASPAR database, fetches TF motifs
- `tf_database_manager.py`: Manages the SQLite database of TF information
- `tf_binding_analyzer.py`: Analyzes gene sequences for TF binding sites
- `create_html_visualization.py`: Generates the HTML report
- `main.py`: Orchestrates the entire analysis workflow
- `genes/`: Directory containing gene FASTA files
- `results/`: Directory for analysis results and visualizations
- `html_report/`: Directory for HTML report files

## Gene Input Format

Place your gene promoter sequence files in the `genes/` directory. Each file should be in FASTA format and named after the gene (e.g., `GAPDH.fa`). The sequence should cover the region from -2000 to +500 relative to the transcription start site (TSS).

Example FASTA format:
```
>GENE_NAME promoter region -2000 to +500
ACTGATCGATCGATCGTAGCTAGCTAGCTACGTAGCTAGCTACGT...
```

## Output Files

- `tf_database.json`: JSON file with all TF information
- `tf_database.sqlite`: SQLite database with structured TF information
- `transcription_factors.csv`: CSV file with all TF information
- `binding_analysis_results.json`: JSON file with binding site analysis results
- `binding_sites_summary.csv`: CSV summary of binding sites per gene per TF
- `binding_heatmap.png`: Heatmap visualization of binding sites
- `binding_site_distribution.png`: Distribution of binding sites relative to TSS
- HTML report files in `html_report/` directory

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

- JASPAR database for providing transcription factor binding profiles
- Biopython for sequence analysis tools  