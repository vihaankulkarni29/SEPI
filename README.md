# SEPI
SEPI (Sensitive Efflux Protein Identifier) 🧬
SEPI is a command-line bioinformatics tool to automate the identification and retrieval of reference protein sequences for specific efflux pumps and their regulators from known antibiotic-sensitive strains of Escherichia coli and Klebsiella pneumoniae.

Overview
The primary function of SEPI is to generate a clean, curated dataset of "wild-type" reference proteins. This dataset can then be used in comparative genomics and proteomics workflows to analyze mutations and variations in antibiotic-resistant isolates. By streamlining this crucial data acquisition step from NCBI, SEPI accelerates research into the mechanisms of antibiotic resistance.

The tool is designed for bioinformaticians, computational biologists, and molecular microbiology researchers studying antibiotic resistance.

Features
Targeted Protein Search: Performs automated searches for a predefined or user-specified list of proteins.

Strain-Specific Retrieval: Prioritizes and retrieves proteins from specific reference strains (e.g., E. coli K-12 MG1655).

Automated NCBI Fetching: Programmatically queries NCBI databases (prioritizing RefSeq) using the E-utilities API.

Robust Fallback Logic: If a strict search yields no results, the tool automatically relaxes its criteria to ensure a sequence is found.

Convenient Output Packaging: Consolidates all retrieved FASTA files into a single .zip archive.

Detailed Accession Reporting: Generates a summary .csv file mapping protein names to their retrieved NCBI accession numbers.

Installation
To get started with SEPI, clone the repository and install the required dependencies.

1. Clone the repository:

Bash

git clone https://github.com/your-username/SEPI.git
cd SEPI
2. Install dependencies:
Make sure you have Python 3.9+ installed. Then, install the required packages using the requirements.txt file.

Bash

pip install -r requirements.txt
Usage
SEPI is run from the command line. The two required arguments are --organism and --email.

Basic Usage
To retrieve all pre-configured proteins for an organism:

For Escherichia coli:

Bash

python sepi.py --organism "ecoli" --email "your.email@example.com"
For Klebsiella pneumoniae:

Bash

python sepi.py --organism "klebsiella" --email "your.email@example.com"
Advanced Usage
To retrieve a specific list of proteins and define a custom output name:

Bash

python sepi.py --organism "ecoli" --proteins "AcrA,AcrB,TolC" --output "Ecoli_AcrAB_TolC" --email "your.email@example.com"
Arguments
--organism: (Required) Target organism. Choices: ecoli, klebsiella.

--email: (Required) Your email address for NCBI API compliance.

--proteins: (Optional) Comma-separated list of protein names. Defaults to "all".

--output: (Optional) The base name for the output files. Defaults to "SEPI_output".

Output
The tool generates two files based on the --output name:

1. ZIP Archive (<output_name>.zip)

A compressed file containing individual FASTA files for each retrieved protein.

Each FASTA file is named using the convention: ProteinName_AccessionNumber.fasta (e.g., AcrB_WP_000068326.1.fasta).

2. CSV Accession Report (<output_name>_accessions.csv)

A two-column table mapping the queried protein name to its retrieved NCBI accession number.

Columns: Protein_Name, Accession_Number
