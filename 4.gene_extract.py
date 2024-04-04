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
    Additionally, it will search the 'note' field if neither gene nor product names match.
    Returns True if the gene was found and extracted, False otherwise.
    """
    try:
        for feat in seq_record.features:
            gene_name = feat.qualifiers.get('gene', [''])[0]  # This should be a string
            product_name = feat.qualifiers.get('product', [''])[0]  # This should be a string

            if gene_name and re.search(gene_pattern, gene_name, re.IGNORECASE):
                sequence = feat.extract(seq_record.seq)
                accession = seq_record.annotations['accessions'][0]
                species = seq_record.annotations['organism'].replace(' ', '_')
                gene_name_for_filename = gene_name.replace(' ', '_')
                output_file = output_dir / f"{base_filename}_{gene_name_for_filename.upper()}.fasta"

                with open(output_file, 'w') as fasta_output:
                    fasta_output.write(f">{species}_{accession} {gene_name.upper()}\n{sequence}\n")
                return True

            elif product_name and re.search(product_pattern, product_name, re.IGNORECASE):
                sequence = feat.extract(seq_record.seq)
                accession = seq_record.annotations['accessions'][0]
                species = seq_record.annotations['organism'].replace(' ', '_')
                gene_or_product = gene_name.upper() if gene_name else product_name.upper()
                gene_or_product_for_filename = gene_or_product.replace(' ', '_')
                output_file = output_dir / f"{base_filename}_{gene_or_product_for_filename}.fasta"

                with open(output_file, 'w') as fasta_output:
                    fasta_output.write(f">{species}_{accession} {gene_or_product}\n{sequence}\n")
                return True

            elif 'note' in feat.qualifiers:  # Checking the note field
                note = feat.qualifiers.get('note', [''])[0]
                if note and re.search(product_pattern, note, re.IGNORECASE):
                    sequence = feat.extract(seq_record.seq)
                    accession = seq_record.annotations['accessions'][0]
                    species = seq_record.annotations['organism'].replace(' ', '_')
                    gene_or_product_or_note = gene_name.upper() if gene_name else (product_name.upper() if product_name else note.upper())
                    gene_or_product_or_note_for_filename = re.sub(r'\W+', '_', gene_or_product_or_note).strip('_')
                    output_file = output_dir / f"{base_filename}_{gene_or_product_or_note_for_filename}.fasta"

                    with open(output_file, 'w') as fasta_output:
                        fasta_output.write(f">{species}_{accession} {gene_or_product_or_note}\n{sequence}\n")
                    return True

    except Exception as e:
        print(f"Error processing {base_filename}: {str(e)}")

    return False


def main():
    parser = argparse.ArgumentParser(description="Extract mitochondrial gene sequences from GenBank files.")
    parser.add_argument("-i", "--input-dir", required=True, help="Path to the directory containing GenBank files.")
    parser.add_argument("-g", "--gene", choices=["COX1", "CYTB", "ND1", "ND2", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "12S", "16S"], default="COX1",
                        help="Specify the mitochondrial gene to extract. Default is COX1.")

    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"The directory {args.input_dir} does not exist or is not a directory.")
        return

    start_time = datetime.datetime.now()

    output_dir = Path.cwd() / f"{args.gene}_sequences"
    output_dir.mkdir(parents=True, exist_ok=True)

    log_contents = []

    genbank_files = [f for f in os.listdir(args.input_dir) if f.lower().endswith('.gb')]
    total_files = len(genbank_files)
    extracted_count = 0

    # Update this pattern to include the new genes
    gene_pattern_dict = {
        "COX1": r'\bCOX 1\b|\bCOX1\b|\bCOXI\b|\bCOX I\b|\bCOI\b',
        "CYTB": r'\bCYTB\b',
        "ND1": r'\bND(1|I)\b',
        "ND2": r'\bND2(2|II)\b',
        "COX2": r'\bCOX 2\b|\bCOX2\b|\bCOXII\b|\bCOX II\b|\bCOII\b',
        "ATP8": r'\bATP8\b',
        "ATP6": r'\bATP6\b',
        "COX3": r'\bCOX 3\b|\bCOX3\b|\bCOXIII\b|\bCOX III\b|\bCOIII\b',
        "ND3": r'\bND(3|III)\b',
        "ND4L": r'\bND(4|IV)L\b',
        "ND4": r'\bND(4|IV)(?!L)\b',
        "ND5": r'\bND(5|V)\b',
        "ND6": r'\bND(6|VI)\b',
        "12S": r'\b12S rRNA\b|\b12S ribosomal RNA\b',
        "16S": r'\b16S rRNA\b|\b16S ribosomal RNA\b'
    }

    gene_pattern = gene_pattern_dict.get(args.gene, args.gene)
    # For mitochondrial genes, the product name often includes the gene name, making a separate product pattern unnecessary for most
    product_pattern = {
        "COX1": r'\bcytochrome c oxidase subunit (I{1}|[01])\b|cytochrome oxidase subunit (I{1}|[01])\b|\bCOX 1\b|\bCOX1\b|\bCOXI\b|\bCOX I\b|\bCOI\b',
        "CYTB": r'cytb|cytochrome b|CYTB',
        "ND1": r'\bND(1|I)\b|NADH dehydrogenase subunit (1|I)\b',
        "ND2": r'\bND(2|II)\b|NADH dehydrogenase subunit (2|II)\b',
        "COX2": r'\b(cytochrome c oxidase subunit (II|2)|cytochrome oxidase subunit (II|2)|\bCOX2\b|\bCOXII\b|\bCOX II\b|\bCOII\b)',
        "ATP8": r'\bATP8\b|ATP synthase F0 subunit 8',
        "ATP6": r'\bATP6\b|ATP synthase F0 subunit 6',
        "COX3": r'\b(cytochrome c oxidase subunit (III|3)|cytochrome oxidase subunit (III|3)|\bCOX3\b|\bCOXIII\b|\bCOX III\b|\bCOIII\b)',
        "ND3": r'\bND(3|III)\b|NADH dehydrogenase subunit (3|III)\b',
        "ND4L": r'\bND(4|IV)L\b|\bNADH dehydrogenase subunit (4|IV) ?L\b|\bNADh dehydrogenase subunit IV ?L\b',
        "ND4": r'\bND(4|IV)(?!L)\b|(?i)\bNADH dehydrogenase subunit (4|IV)\b(?! ?L)',
        "ND5": r'\bND(5|V)\b|NADH dehydrogenase subunit (5|V)\b',
        "ND6": r'\bND(6|VI)\b|NADH dehydrogenase subunit (6|VI)\b',
        "12S": r'\b12S rRNA\b|\b12S ribosomal RNA\b',
        "16S": r'\b16S rRNA\b|\b16S ribosomal RNA\b'
    }.get(args.gene, '')

    for genbank_file in genbank_files:
        genbank_path = os.path.join(args.input_dir, genbank_file)
        base_filename = os.path.splitext(genbank_file)[0]

        with open(genbank_path, 'r') as handle:
            for seq_record in SeqIO.parse(handle, "genbank"):
                if extract_gene(seq_record, gene_pattern, product_pattern, output_dir, base_filename):
                    species = seq_record.annotations['organism']
                    log_contents.append((species, base_filename))
                    extracted_count += 1
                    break
                else:
                    log_contents.append((None, base_filename))

    log_contents.sort(key=lambda x: x[0].lower() if x[0] else "")

    log_file_path = Path.cwd() / f"{args.gene}_extraction.log"

    with open(log_file_path, 'w') as log_file:
        log_file.write(f"Run started at: {start_time}\n")
        for species, filename in log_contents:
            if species:
                log_file.write(f"{args.gene} sequence extracted for {filename} (Species: {species})\n")
            else:
                log_file.write(f"No {args.gene} gene found in {filename}.gb\n")

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
