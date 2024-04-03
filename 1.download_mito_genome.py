import argparse
import textwrap
from Bio import Entrez
import time
import datetime
import os

def download_mt_genomes(species_names, email, output_dir, log_file, rettype):
    Entrez.email = email
    start_time = datetime.datetime.now()
    log_file.write(f"Script started at {start_time}\n")

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for species in species_names:
        species_start_time = datetime.datetime.now()
        try:
            # Build the search query
            search_query = f"{species}[orgn] AND mitochondrion[filter] AND complete genome[title]"
            search_handle = Entrez.esearch(db="nucleotide", term=search_query)
            search_results = Entrez.read(search_handle)
            search_handle.close()

            # Check if sequences are found
            if not search_results["IdList"]:
                log_file.write(f"No mitochondrial genomes found for {species} (checked at {species_start_time})\n")
                continue

            downloaded = 0
            # Download and save sequences
            for seq_id in search_results["IdList"]:
                fetch_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype=rettype, retmode="text")
                seq_record = fetch_handle.read()
                fetch_handle.close()
                # Get the ACCESSION number
                accession_number = seq_record.split("ACCESSION")[1].split()[0]
                filename = os.path.join(output_dir, f"{species.replace(' ', '_')}_{accession_number}.gb")
                with open(filename, "w") as file:
                    file.write(seq_record)
                downloaded += 1

            log_file.write(f"Downloaded {downloaded} sequences for {species} (started at {species_start_time})\n")
            # Avoid making requests too frequently
            time.sleep(1)

        except Exception as e:
            log_file.write(f"An error occurred while processing {species}: {e} (at {species_start_time})\n")

    end_time = datetime.datetime.now()
    log_file.write(f"Script finished at {end_time}, total duration: {end_time - start_time}\n")

def load_species_list(filename):
    with open(filename, 'r') as file:
        return [line.strip() for line in file if line.strip()]

def extract_species_info_from_file(file_path):
    downloaded_info = {}
    not_downloaded_species = []

    with open(file_path, 'r') as file:
        for line in file:
            if "Downloaded" in line:
                parts = line.split(' sequences for ')
                species_name = parts[1].split(' (')[0]
                sequence_count = int(parts[0].split()[1])
                downloaded_info[species_name] = sequence_count
            elif "No mitochondrial genomes found" in line:
                species_name = line.split(' for ')[1].split(' (')[0]
                not_downloaded_species.append(species_name)

    return downloaded_info, not_downloaded_species

def save_statistics_to_tsv(downloaded_info, not_downloaded_species, output_tsv_path):
    with open(output_tsv_path, 'w') as tsv_file:
        tsv_file.write("Species\tSequence Count\n")
        for species_name, sequence_count in downloaded_info.items():
            tsv_file.write(f"{species_name}\t{sequence_count}\n")
        
        for species_name in sorted(not_downloaded_species):
            tsv_file.write(f"{species_name}\t0\n")

def compute_statistics(downloaded_info, not_downloaded_species):
    total_downloaded_species = len(downloaded_info)
    total_not_downloaded_species = len(not_downloaded_species)
    total_sequences = sum(downloaded_info.values())

    avg_sequences_per_species = total_sequences / total_downloaded_species if total_downloaded_species > 0 else 0
    most_sequences_species = max(downloaded_info, key=downloaded_info.get) if downloaded_info else None
    least_sequences_species = min(downloaded_info, key=downloaded_info.get) if downloaded_info else None

    sorted_downloaded_species = sorted(downloaded_info.keys())
    sorted_not_downloaded_species = sorted(not_downloaded_species)

    return {
        "Total downloaded species": total_downloaded_species,
        "Total not downloaded species": total_not_downloaded_species,
        "Total sequences downloaded": total_sequences,
        "Average sequences per downloaded species": avg_sequences_per_species,
        "Species with most sequences": most_sequences_species,
        "Species with least sequences": least_sequences_species,
        "Downloaded species (sorted)": sorted_downloaded_species,
        "Not downloaded species (sorted)": sorted_not_downloaded_species
    }

def save_statistics_to_file(statistics, downloaded_info, output_file_path):
    with open(output_file_path, 'w') as file:
        for stat, value in statistics.items():
            if isinstance(value, list):
                file.write(f"{stat}:\n")
                for item in value:
                    file.write(f"  {item}\n")
            else:
                file.write(f"{stat}: {value}\n")

        file.write("\nDownloaded Species Sequence Counts:\n")
        for species in sorted(downloaded_info.keys()):
            file.write(f"  {species}: {downloaded_info[species]} sequences\n")

def main():
    parser = argparse.ArgumentParser(description="Download mitochondrial genomes for a list of species from NCBI and analyze the results.", 
                                     usage='use "python %(prog)s --help" for more information',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-l", "--list", dest="species_filename", required=True, 
                        help= textwrap.dedent('''The file of species list.Each species should be on a separate line, e.g.,
		            Elachura formosa
		            Troglodytes troglodytes
                            Sitta europaea 
                                  ... 
                                               '''))
    parser.add_argument("-e", "--email", dest="email", required=True, help="Your email to identify yourself to NCBI.")
    parser.add_argument("-o", "--output", dest="output_directory", default="downloaded_mtgenomes",
                    help="Directory where downloaded genomes will be saved.")
    parser.add_argument("-t", "--type", dest="rettype", default="gb", help="Type of data to be retrieved from NCBI. gb or fa, default is 'gb'.")
    parser.add_argument("--log", dest="log_filename", default="mito_download.log", help="Name of the log file. Default is 'mito_download.log'.")

    args = parser.parse_args()

    species_list = load_species_list(args.species_filename)
    
    output_file_path = "species_download_stats.txt"
    tsv_output_file_path = "species_download_counts.tsv"
    
    with open(args.log_filename, 'w') as log_file:
        download_mt_genomes(species_list, args.email, args.output_directory, log_file, args.rettype)

    downloaded_info, not_downloaded_species = extract_species_info_from_file(args.log_filename)
    statistics = compute_statistics(downloaded_info, not_downloaded_species)

    # Saving statistics to a file
    save_statistics_to_file(statistics, downloaded_info, output_file_path)
    print(f"Statistics saved to {output_file_path}")

    save_statistics_to_tsv(downloaded_info, not_downloaded_species, tsv_output_file_path)
    print(f"TSV file saved to {tsv_output_file_path}")

    print(f"#####")
    print(f"Mitochondrial genomes download finished. Please check sequences file in [{args.output_directory}], statistics of download result in [{output_file_path}], [{tsv_output_file_path}]")
    print(f"#####")

if __name__ == "__main__":
    main()

