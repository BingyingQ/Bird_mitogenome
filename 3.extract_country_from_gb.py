import os
from Bio import SeqIO
import csv
from datetime import datetime
import argparse

def extract_country_info(genbank_folder, output_file):
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter='\t')
        csv_writer.writerow(['Species', 'ACCESSION', 'Country', 'Source'])

        for filename in os.listdir(genbank_folder):
            if filename.endswith(".gb"):
                filepath = os.path.join(genbank_folder, filename)

                for record in SeqIO.parse(filepath, 'genbank'):
                    species_name = record.annotations.get('organism', '')
                    accession = record.id
                    country = ''
                    source = ''  # Added to track the source of the country information

                    for feature in record.features:
                        if 'country' in feature.qualifiers:
                            country = feature.qualifiers['country'][0]
                            source = 'feature'
                            break

                    if not country:
                        earliest_date = datetime.max
                        for reference in record.annotations['references']:
                            if reference.journal.startswith('Submitted'):
                                date_str = reference.journal.split('(')[-1].split(')')[0].strip()
                                date_format = "%d-%b-%Y"
                                try:
                                    date = datetime.strptime(date_str, date_format)
                                    if date < earliest_date:
                                        earliest_date = date
                                        country = reference.journal.split()[-1]
                                        source = 'journal'
                                except ValueError:
                                    pass

                    if not country:
                        for reference in record.annotations['references']:
                            if reference.journal.startswith('Submitted'):
                                journal_parts = reference.journal.split()
                                if len(journal_parts) > 0:
                                    country = journal_parts[-1]
                                    source = 'journal'
                                break

                    if country and ',' in country:
                        country = country.split(',')[-1].strip()

                    csv_writer.writerow([species_name, accession, country, source])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract country information from GenBank files.')
    parser.add_argument('-i', '--input_folder', default='./downloaded_mtgenomes_gb/', help='Path to the GenBank folder')
    parser.add_argument('-o', '--output_file', default='country_output.tsv', help='Path to the output CSV file')
    args = parser.parse_args()

    extract_country_info(args.input_folder, args.output_file)

