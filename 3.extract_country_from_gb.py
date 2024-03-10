import os
from Bio import SeqIO
import csv
from datetime import datetime
import argparse
import pandas as pd

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

def analyze_china_presence(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t")

    # 提取物种名的前两部分，创建新的列
    df['Species_prefix'] = df['Species'].apply(lambda x: ' '.join(x.split()[:2]))

    # 按前两部分物种名分组
    grouped_species = df.groupby('Species_prefix')

    # 创建结果的DataFrame
    result_df = pd.DataFrame(columns=['Species_prefix', 'China_presence'])

    # 遍历每个物种的所有行
    for species_prefix, group in grouped_species:
        has_china_sequence = 'Yes' if 'China' in group['Country'].values else 'No'
        result_df = result_df.append({'Species_prefix': species_prefix, 'China_presence': has_china_sequence}, ignore_index=True)

    # 输出结果到新的CSV文件
    result_df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract and analyze country information from GenBank files.')
    parser.add_argument('-i', '--input_folder', default='./downloaded_mtgenomes_gb/', help='Path to the GenBank folder')
    parser.add_argument('-o', '--output_file', default='country_output.csv', help='Path to the output CSV file')
    parser.add_argument('-c', '--china_output_file', default='china_or_not_result.csv', help='The CSV file presented whether these species contains sequences from China')
    args = parser.parse_args()

    # 调用提取国家信息的函数
    extract_country_info(args.input_folder, args.output_file)

    # 调用分析China是否存在的函数
    analyze_china_presence(args.output_file, args.china_output_file)
