#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import os
import re
import argparse
from Bio import SeqIO
from pathlib import Path
import datetime


def extract_gene(seq_record, gene_pattern, product_pattern, output_dir, base_filename):
    """
    Attempts to extract the gene sequence from the SeqRecord object if it matches the provided gene pattern (regular expression).
    If the gene name does not match, it will try to match the product name.
    Returns True if the gene was found and extracted, False otherwise.
    """
    try:
        for feat in seq_record.features:
            gene_name = feat.qualifiers.get('gene', [''])[0]
            product_name = feat.qualifiers.get('product', [''])[0]

            # Try matching gene name
            if gene_name and re.search(gene_pattern, gene_name, re.IGNORECASE):
                sequence = feat.extract(seq_record.seq)
                accession = seq_record.annotations['accessions'][0]
                species = seq_record.annotations['organism'].replace(' ', '_')
                output_file = output_dir / f"{base_filename}_{gene_name.upper()}.fasta"

                with open(output_file, 'w') as fasta_output:
                    fasta_output.write(f">{species}_{accession} {gene_name.upper()}\n{sequence}\n")
                return True

            # If gene name doesn't match, try matching product name
            elif product_name and re.search(product_pattern, product_name, re.IGNORECASE):
                sequence = feat.extract(seq_record.seq)
                accession = seq_record.annotations['accessions'][0]
                species = seq_record.annotations['organism'].replace(' ', '_')
                output_file = output_dir / f"{base_filename}_{gene_name.upper()}.fasta"

                with open(output_file, 'w') as fasta_output:
                    fasta_output.write(f">{species}_{accession} {gene_name.upper()}\n{sequence}\n")
                return True

    except Exception as e:
        print(f"Error processing {base_filename}: {str(e)}")

    return False


def main():
    parser = argparse.ArgumentParser(description="Extract COX1 or CYTB gene sequences from GenBank files.")
    parser.add_argument("-i", "--input-dir", required=True, help="Path to the directory containing GenBank files.")
    parser.add_argument("-g", "--gene", choices=["COX1", "CYTB"], default="COX1",
                        help="Specify the gene to extract (COX1 or CYTB). Default is COX1.")

    args = parser.parse_args()

    # Validate input directory
    if not os.path.isdir(args.input_dir):
        print(f"The directory {args.input_dir} does not exist or is not a directory.")
        return

    # Record the start time
    start_time = datetime.datetime.now()

    # Set the output directory within the current working directory
    output_dir = Path.cwd() / f"{args.gene}_sequences"
    # Create the output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize log contents list
    log_contents = []

    # Get all .gb files in the input directory
    genbank_files = [f for f in os.listdir(args.input_dir) if f.lower().endswith('.gb')]
    total_files = len(genbank_files)
    extracted_count = 0

    # Regular expression patterns for gene names, to account for various spellings and typos
    gene_pattern = r'\bCOX 1\b|\bCOX1\b|\bCOXI\b|\bCOX I\b|\bCOI\b' if args.gene == "COX1" else r'CYTB|cytb'
    product_pattern = (
        r'\bcytochrome c oxidase subunit (I{1}|[01])\b|cytochrome oxidase subunit (I{1}|[01])\b|\bCOX 1\b|\bCOX1\b|\bCOXI\b|\bCOX I\b|\bCOI\b'
        if args.gene == "COX1"
        else r'cytb|cytochrome b|CYTB'
    )

    # Process each file
    for genbank_file in genbank_files:
        genbank_path = os.path.join(args.input_dir, genbank_file)
        base_filename = os.path.splitext(genbank_file)[0]

        with open(genbank_path, 'r') as handle:
            for seq_record in SeqIO.parse(handle, "genbank"):
                # Attempt to extract the gene
                if extract_gene(seq_record, gene_pattern, product_pattern, output_dir, base_filename):
                    species = seq_record.annotations['organism']
                    log_contents.append((species, base_filename))
                    extracted_count += 1
                    break
                else:
                    log_contents.append((None, base_filename))

    # Sort log contents by species name
    log_contents.sort(key=lambda x: x[0].lower() if x[0] else "")

    # Log file setup
    log_file_path = Path.cwd() / f"{args.gene}_extraction.log"

    # Write to log file
    with open(log_file_path, 'w') as log_file:
        log_file.write(f"Run started at: {start_time}\n")
        for species, filename in log_contents:
            if species:
                log_file.write(f"{args.gene} sequence extracted for {filename} (Species: {species})\n")
            else:
                log_file.write(f"No {args.gene} gene found in {filename}.gb\n")

        # Record the end time and calculate the duration
        end_time = datetime.datetime.now()
        duration = end_time - start_time
        log_file.write(f"\nRun ended at: {end_time}\n")
        log_file.write(f"Total run time: {duration}\n")
        log_file.write(f"Total files processed: {total_files}\n")
        log_file.write(f"Total {args.gene} sequences extracted: {extracted_count}\n")

    print(f"All {args.gene} sequences have been successfully extracted to {output_dir}")
    print(f"The log has been saved to {log_file_path}")


if __name__ == "__main__":
    main()

