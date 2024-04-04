import os
import argparse
from Bio import SeqIO

def delete_duplicate_sequences(directory_path):
    # 遍历目录中的所有文件
    all_files = os.listdir(directory_path)

    # 分离出所有以NC开头的文件和其他文件
    nc_files = [f for f in all_files if f.startswith('NC')]
    other_files = [f for f in all_files if not f.startswith('NC')]

    # 读取并存储所有文件的序列
    species_sequences = {}
    file_to_species = {}
    for file in all_files:
        species_name = '_'.join(file.split('_')[:2])
        filepath = os.path.join(directory_path, file)
        file_to_species[file] = species_name
        try:
            for record in SeqIO.parse(filepath, "genbank"):
                seq_str = str(record.seq)
                if species_name in species_sequences:
                    if seq_str in species_sequences[species_name]:
                        # 如果序列已存在，则标记当前文件为重复文件
                        species_sequences[species_name][seq_str].append(file)
                    else:
                        species_sequences[species_name][seq_str] = [file]
                else:
                    species_sequences[species_name] = {seq_str: [file]}
        except Exception as e:
            print(f"Error processing file {file}: {e}. Skipping...")

    # 删除重复的序列文件，只保留一份
    for species, seqs in species_sequences.items():
        for seq, files in seqs.items():
            if len(files) > 1:
                # 除了第一个文件外，其他都删除
                for file_to_delete in files[1:]:
                    filepath = os.path.join(directory_path, file_to_delete)
                    try:
                        os.remove(filepath)
                        print(f"Deleted duplicate: {file_to_delete}")
                    except Exception as e:
                        print(f"Error deleting file {file_to_delete}: {e}. Skipping...")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Delete duplicate sequences from a directory.')
    parser.add_argument('-i', '--input', type=str, help='Input directory path')
    args = parser.parse_args()

    if args.input:
        directory_path = args.input

    delete_duplicate_sequences(directory_path)
