from Bio import SeqIO
import os
import csv
import argparse

def extract_species_and_accession_from_genbank(genbank_path):
    with open(genbank_path, "r") as file:
        for record in SeqIO.parse(file, "genbank"):
            organism = record.annotations.get("organism", "Unknown")
            accession = record.id
            return organism, accession

def count_sequences(folder_path):
    species_sequences = {}

    for filename in os.listdir(folder_path):
        if filename.endswith(".gb"):
            file_path = os.path.join(folder_path, filename)

            # Extract species name and accession number from GenBank file
            species_name, accession_number = extract_species_and_accession_from_genbank(file_path)

            # Initialize or update the dictionary entry for the species
            if species_name not in species_sequences:
                species_sequences[species_name] = {"count": 0, "accessions": []}
            
            species_sequences[species_name]["count"] += 1
            species_sequences[species_name]["accessions"].append(accession_number)

    return species_sequences

def write_to_tsv(result, output_file):
    with open(output_file, "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file, delimiter='\t')
        
        # Write header
        writer.writerow(["Species", "Sequence Count", "Accession Numbers"])

        # Write data
        for species, info in result.items():
            writer.writerow([species, info["count"], ",".join(info["accessions"])])

def parse_arguments():
    parser = argparse.ArgumentParser(description="Count sequences and list accession numbers for each species in GenBank files.")
    parser.add_argument("-i", "--input_folder", required=True, help="Input folder containing GenBank files.")
    parser.add_argument("-o", "--output_file", default="./count_result.tsv", help="Output TSV file path.")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    folder_path = args.input_folder
    output_file = args.output_file
    result = count_sequences(folder_path)

    # Write the result to TSV file
    write_to_tsv(result, output_file)

    print(f"Results written to {output_file}")

