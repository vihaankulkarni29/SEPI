#!/usr/bin/env python3

"""
SEPI (Sensitive Efflux Protein Identifier)
Version: 1.2 (Final fallback search logic)
Author: Gemini (based on PRD by Vihaan)
Date: August 29, 2025

Description:
A command-line bioinformatics tool to automate the identification and retrieval
of reference protein sequences for specific efflux pumps and their regulators from
known antibiotic-sensitive strains of Escherichia coli and Klebsiella pneumoniae.
"""

import argparse
import logging
import os
import sys
import time
import zipfile
from io import StringIO
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO

# --- Pre-configured Protein Lists (as per PRD section 3.2) ---

PROTEIN_CONFIG = {
    "ecoli": {
        "proteins": [
            # AcrAB-TolC Pump & AcrZ
            "AcrA", "AcrB", "TolC", "AcrZ",
            # Transcriptional Regulators
            "AcrR", "MarA", "MarR", "RamA", "RamR", "SoxS", "Rob", "EnvR",
            # Other Efflux Pumps
            "AcrD", "AcrE", "AcrF", "MdtB", "MdtC"
        ],
        # Strain-specific requirements for higher accuracy
        "strain_specific": {
            "AcrA": "K-12 MG1655",
            "AcrB": "K-12 MG1655"
        },
        "organism_name": "Escherichia coli"
    },
    "klebsiella": {
        "proteins": [
            "OqxA", "OqxB", "EefA", "EefB", "EefC", "KexD", "KexE", "KexF"
        ],
        "strain_specific": {}, # No specific strains mentioned in PRD
        "organism_name": "Klebsiella pneumoniae"
    }
}

# --- Logging Configuration ---

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def search_and_fetch_protein(protein_name: str, organism: str, email: str) -> tuple[str, str] | None:
    """
    Constructs a query, searches NCBI, and fetches the top protein sequence.
    Includes a fallback mechanism to relax search criteria if the initial search fails.
    """
    Entrez.email = email
    config = PROTEIN_CONFIG[organism]
    organism_full_name = config["organism_name"]

    base_query_parts = [
        f'"{organism_full_name}"[Organism]',
        f'"{protein_name}"[Protein Name]',
    ]
    if strain := config["strain_specific"].get(protein_name):
        base_query_parts.insert(1, f'"{strain}"[Strain]')

    # Define a series of queries from most to least strict
    queries = [
        # 1. Strict: RefSeq, complete genome, NOT resistance
        " AND ".join(base_query_parts + ['("complete genome"[Filter])', '(srcdb_refseq[PROP])', 'NOT (resistance OR resistant OR multidrug OR hypothetical)']),
        # 2. Relaxed: RefSeq, complete genome (allows resistance terms)
        " AND ".join(base_query_parts + ['("complete genome"[Filter])', '(srcdb_refseq[PROP])']),
        # 3. Broad: All protein DB, complete genome
        " AND ".join(base_query_parts + ['("complete genome"[Filter])']),
        # 4. Broadest: All protein DB, no genome filter
        " AND ".join(base_query_parts)
    ]

    protein_id = None
    for i, query in enumerate(queries):
        try:
            logging.info(f"Searching for '{protein_name}' (Attempt {i+1})...")
            
            handle = Entrez.esearch(db="protein", term=query, retmax=1)
            record = Entrez.read(handle)
            handle.close()
            time.sleep(0.35)  # Adhere to NCBI's limit

            if record["IdList"]:
                protein_id = record["IdList"][0]
                logging.info(f"Found a candidate for '{protein_name}' on Attempt {i+1}.")
                break # Exit the loop once an ID is found
        except Exception as e:
            logging.error(f"An error occurred during search attempt {i+1} for '{protein_name}': {e}")
            continue # Try next query

    if not protein_id:
        logging.warning(f"All search attempts failed for '{protein_name}'. No entry found.")
        return None

    # Fetch the record using the found ID
    try:
        fetch_handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        fasta_data = fetch_handle.read()
        fetch_handle.close()
        time.sleep(0.35)

        fasta_io = StringIO(fasta_data)
        seq_record = SeqIO.read(fasta_io, "fasta")
        accession = seq_record.id

        logging.info(f"Successfully retrieved '{protein_name}' with Accession: {accession}")
        return accession, fasta_data
    except Exception as e:
        logging.error(f"Found ID {protein_id} but failed to fetch the record for '{protein_name}': {e}")
        return None


def main():
    """Main function to parse arguments and run the SEPI workflow."""
    parser = argparse.ArgumentParser(
        description="SEPI: A tool to fetch sensitive efflux protein reference sequences from NCBI.",
        epilog="Example: python sepi.py --organism ecoli --proteins AcrA,AcrB,TolC --output Ecoli_AcrAB_TolC --email user@example.com"
    )
    parser.add_argument(
        "--organism",
        choices=["ecoli", "klebsiella"],
        required=True,
        help="Specifies the target organism."
    )
    parser.add_argument(
        "--proteins",
        type=str,
        default="all",
        help='A comma-separated list of protein names, or "all" for the pre-configured list. (Default: all)'
    )
    parser.add_argument(
        "--output",
        type=str,
        default="SEPI_output",
        help="The base name for the output files. (Default: SEPI_output)"
    )
    parser.add_argument(
        "--email",
        required=True,
        help="Your email address (required by NCBI for API usage)."
    )

    args = parser.parse_args()

    # Determine the list of proteins to retrieve
    if args.proteins.lower() == "all":
        target_proteins = PROTEIN_CONFIG[args.organism]["proteins"]
    else:
        target_proteins = [p.strip() for p in args.proteins.split(',')]

    if not target_proteins:
        logging.error("No protein names provided. Exiting.")
        sys.exit(1)

    logging.info(f"Starting SEPI for organism: '{args.organism}'")
    logging.info(f"Target proteins: {', '.join(target_proteins)}")
    logging.info(f"Output file base name: '{args.output}'")

    # --- Workflow Execution ---
    results = []
    output_dir = Path(f"{args.output}_fasta_files")
    output_dir.mkdir(exist_ok=True)

    for protein in target_proteins:
        fetched_data = search_and_fetch_protein(protein, args.organism, args.email)
        if fetched_data:
            accession, fasta_sequence = fetched_data
            results.append({"Protein_Name": protein, "Accession_Number": accession})

            # Save individual FASTA file (as per PRD section 6.1)
            file_name = f"{protein}_{accession}.fasta"
            file_path = output_dir / file_name
            with open(file_path, "w") as f_out:
                f_out.write(fasta_sequence)

    if not results:
        logging.warning("No proteins were successfully retrieved. No output files will be generated.")
        # Clean up empty directory
        try:
            os.rmdir(output_dir)
        except OSError:
            pass # Directory might not be empty if user created files in it
        sys.exit(0)

    # --- Generate Output Files (as per PRD section 6) ---
    
    # 1. Create CSV Accession Report
    csv_filename = f"{args.output}_accessions.csv"
    try:
        df = pd.DataFrame(results)
        df.to_csv(csv_filename, index=False)
        logging.info(f"Accession report saved to '{csv_filename}'")
    except Exception as e:
        logging.error(f"Failed to create CSV report: {e}")

    # 2. Create ZIP Archive
    zip_filename = f"{args.output}.zip"
    try:
        with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for fasta_file in output_dir.glob("*.fasta"):
                zipf.write(fasta_file, arcname=fasta_file.name)
                os.remove(fasta_file) # Clean up individual file after adding to zip
        os.rmdir(output_dir) # Clean up the temporary directory
        logging.info(f"FASTA files packaged into '{zip_filename}'")
    except Exception as e:
        logging.error(f"Failed to create ZIP archive: {e}")

    logging.info("SEPI run completed successfully. ✅")


if __name__ == "__main__":
    main()