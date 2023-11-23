#!/usr/bin/python3

import os
import sys
import requests
import subprocess
import xml.etree.ElementTree as ET
import json
import time


'''SET WORKING DIRECTORY'''
def set_working_directory():
    script_path = os.path.abspath(__file__)
    
    #getting current directory of the script
    script_dir = os.path.dirname(script_path)

    #change working dir
    os.chdir(script_dir)
    print(f"working directory changed to: {script_dir}")

'''RUN BASH COMMANDS FUNCTION'''
def run_command(command):
    try:
        # if command is a list, convert it to a string
        if isinstance(command, list):
            command = ' '.join(command)
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e.stderr}")
        sys.exit(1)
      
'''RUN RESTFUL API'''  
def run_get(url):
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Error making request: {response.status_code}")
        print(response.text)
        return None
    return response.json()

'''DRAW A TAXONOMIC TREE'''
def get_taxonomic_tree(taxonomic_group):
    # Define the base URL for NCBI E-utilities
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # Perform the search to get the NCBI Taxonomy ID
    search_url = f"{base_url}/esearch.fcgi"
    search_params = {
        "db": "taxonomy",
        "term": taxonomic_group,
        "retmode": "xml"
    }
    search_response = requests.get(search_url, params=search_params)
    search_root = ET.fromstring(search_response.text)
    
    tax_id_element = search_root.find(".//Id")
    if tax_id_element is None:
        print(f"Taxonomic ID for '{taxonomic_group}' not found.")
        return ""

    tax_id = tax_id_element.text

    # Get the lineage using the fetched taxonomic ID
    fetch_url = f"{base_url}/efetch.fcgi"
    fetch_params = {
        "db": "taxonomy",
        "id": tax_id,
        "retmode": "xml"
    }
    fetch_response = requests.get(fetch_url, params=fetch_params)
    lineage_root = ET.fromstring(fetch_response.text)

    # Extract lineage information
    lineage_ids = [taxon.find('TaxId').text for taxon in lineage_root.findall(".//LineageEx/Taxon")]
    lineage_names = [taxon.find('ScientificName').text for taxon in lineage_root.findall(".//LineageEx/Taxon")]

    # Generate a tree view as a string
    tree_str = ""
    for depth, (taxid, name) in enumerate(zip(lineage_ids, lineage_names)):
        tree_str += " " * depth + f"- {name} ('{taxid}')\n"
    
    return tree_str

'''FETCH ALL SEQUENCES AND WRITE IN ON FILE, WITH A LIMIT'''
def fetch_sequences_all(taxonomy, protein_family, limit):
    query = f'"{protein_family}[Protein Name] AND {taxonomy}[Organism]"'
    command = f"esearch -db protein -query {query} -retmax {limit} | efetch -format fasta"
    protein_sequences = run_command(command)

    # Save sequences to a file
    folder_path = f'{taxonomy}/{protein_family}'
    os.makedirs(folder_path, exist_ok=True)
    with open(f'{folder_path}/sequences.fasta', 'w') as file:
        file.write(protein_sequences)

'''GET SPECIES~PROTEIN_IDS LIST'''
def get_species_protein_ids(taxonomy, protein_family, limit):
    # https request
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    esearch_params = {
        'db': 'protein',
        'term': f"{protein_family}[Protein Name] AND {taxonomy}[Organism]",
        'retmax': limit,
        'retmode': 'xml',
    }
    esearch_response = requests.get(esearch_url, params=esearch_params)
    if esearch_response.status_code != 200 or not esearch_response.text:
        print("Error fetching data from ESearch")
        return {}

    # parse esearch results to extract id list
    root = ET.fromstring(esearch_response.text)
    ids = [id_elem.text for id_elem in root.findall('.//Id')]

    if not ids:
        print("No IDs found in ESearch result")
        return {}

    # efetch results
    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    efetch_params = {
        'db': 'protein',
        'id': ','.join(ids),
        'retmode': 'xml',
    }
    efetch_response = requests.get(efetch_url, params=efetch_params)
    if efetch_response.status_code != 200 or not efetch_response.text:
        print("Error fetching data from EFetch")
        return {}

    # parse efetch results，establishing mapping relationship between species and ids
    species_protein_map = {}
    efetch_root = ET.fromstring(efetch_response.text)
    for seq_entry in efetch_root.findall('.//GBSeq'):
        protein_id = seq_entry.find('.//GBSeq_primary-accession').text
        species_name = seq_entry.find('.//GBSeq_organism').text
        species_protein_map[species_name] = protein_id

    return species_protein_map

'''PARSE DICTIONARY OBJ TO FASTA FORMAT'''
def to_fasta_format(seq_dict):
    fasta_format_str = ""
    for species, seq in seq_dict.items():
        # Add the identifier line
        fasta_format_str += f">{species}\n"
        # Add the sequence line
        fasta_format_str += f"{seq}\n"
    return fasta_format_str

'''DOWNLOAD SEQUENCES BASED ON CHOSEN SPECIES'''
def fetch_sequences(species_protein_map, taxonomy, folder_path):
    sequences = {}
    for species, protein_id in species_protein_map.items():
        # fetch the protein sequence using efetch
        protein_sequence = run_command(["efetch", "-db", "protein", "-id", protein_id, "-format", "fasta"])

        # write the sequence to a file
        file_path = os.path.join(folder_path, f"{species}.fasta")
        with open(file_path, 'w') as file:
            file.write(protein_sequence)

        # Save the sequence to a dictionary
        sequences[species] = protein_sequence
    return sequences

'''EMOBOSS CONSERVATION ANALYSIS'''
def run_emboss_conservation(taxonomy, protein_family, sequences):
    # to analyze sequence conservation.
    try:
    # 1st to run clustal omega sequence alignment
        subprocess.run(["clustalo", "-i", f'{taxonomy}/{protein_family}/sequences.fasta', "-o", f"{taxonomy}/{protein_family}/aligned_sequences.fasta", "--force", "--outfmt=fasta"])

    # 2nd using plotcon for conservation analysis and generate png plot
        subprocess.run(["plotcon", "-sequence", f"{taxonomy}/{protein_family}/aligned_sequences.fasta", 
                "-winsize", "4", "-graph", "png", "-goutfile", f"{taxonomy}/{protein_family}/conservation_plot"])
    except Exception as e:
        print(f"Error Msg: {e}")

'''PARSE SEQUENCES FROM FASTA FORMAT FILE'''
def parse_fasta_sequences(filename):
    with open(filename, 'r') as file:
        sequences = []
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                sequence = ''
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    return sequences

'''SIMPLE CONSERVATION ANALYSIS BY CALCULATING SCORE OF EACH POSITION'''
def calculate_conservation(sequences):
    sequence_length = len(sequences[0])
    num_sequences = len(sequences)
    
    conservation_scores = []
    for i in range(sequence_length):
        position = [seq[i] for seq in sequences]
        most_common_residue = max(set(position), key=position.count)
        frequency = position.count(most_common_residue) / num_sequences
        conservation_scores.append(frequency)
    return conservation_scores


'''FUNCTION FOR PROSITE SCANNING'''
def scan_prosite_motifs(sequence_file, folder_path, species):
    # create the output directory if it does not exist
    output_dir = os.path.join(folder_path, "Scan_reports")
    os.makedirs(output_dir, exist_ok=True)

    # define the path for the result file
    result_file = os.path.join(output_dir, f"{species}.txt")
    command = f"patmatmotifs -sequence \"{sequence_file}\" -outfile \"{result_file}\" -full"
    print(f"Scan for {sequence_file} begins.\n")
    result =  run_command(command)
    print(f"Scan result for {species} saved as {result_file}\n")
    return result_file

def scan_all_sequences(selected_species, folder_path):
    results = {}
    for species in selected_species.keys():
        sequence_file = os.path.join(folder_path, f"{species}.fasta")
        results[species] = scan_prosite_motifs(sequence_file, folder_path, species)
        with open(results[species], 'r') as file:
            print(f"Results for {species}:\n")
            for line in file:
                # print line by line
                print(line, end='')
            print("\n")
    # emptyness check
    if not results:
        print("No similar motifs detected.")
        return {}
    else:
        return results

'''FIND HIGH SIMILARITY MOTIFS'''
def find_high_similarity_motifs(folder_path, similarity_threshold):
    motif_data = []

    # List all .txt files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt'):
            file_path = os.path.join(folder_path, filename)

            # Read the file and extract motifs
            with open(file_path, 'r') as file:
                for line in file:
                    if "some_condition_to_identify_motif_lines" in line:
                        # Extract motif and similarity score
                        motif, similarity = extract_motif_and_similarity(line)
                        if similarity >= similarity_threshold:
                            motif_data.append((motif, similarity))

    return motif_data

'''Phylogenetic tree'''  
'''merge all fasta files into one''' 
def merge_fasta_files(input_folder, output_file):
    with open(output_file, 'w') as merged_file:
        for file in os.listdir(input_folder):
            if file.endswith(".fasta") or file.endswith(".fa"):
                file_path = os.path.join(input_folder, file)
                with open(file_path, 'r') as fasta_file:
                    merged_file.write(fasta_file.read() + "\n")

'''convert FASTA format to phylip format.'''
def convert_fasta_to_phylip(input_fasta, output_phylip):
    with open(input_fasta, 'r') as fasta_file, open(output_phylip, 'w') as phylip_file:
        sequences = {}
        current_seq_name = None

        for line in fasta_file:
            line = line.strip()
            if not line:
                continue # ignore null line
            if line.startswith(">"):
                current_seq_name = line[1:].split()[0]  # extract sequence name
                sequences[current_seq_name] = []
            elif current_seq_name:
                sequences[current_seq_name].append(line)
            else:
                print("Error in FASTA file format: Sequence data found before sequence name.")
                return

        if not sequences:
            print("No sequences found in the FASTA file.")
            return

        num_sequences = len(sequences)
        sequence_length = len(''.join(sequences[next(iter(sequences))]))
        phylip_file.write(f"{num_sequences} {sequence_length}\n")

        for seq_name, seq_parts in sequences.items():
            sequence = ''.join(seq_parts)
            phylip_file.write(f"{seq_name[:10].ljust(10)} {sequence}\n")

'''run phyml to generate a phylogenetic tree.'''
def run_phyml(input_phylip, output_tree):
    print("run phyml command")
    command = f"phyml -i {input_phylip} -o {output_tree} -d nt -b 100 -m GTR"
    try:
        result = run_command(command)
        print(f"Phylogenetic tree generated at: {output_tree}")
    except Exception as e:
        print(f"Error generating phylogenetic tree: {e}")

'''SWISS-MODEL'''
def predict_structure(file_path):
    with open(file_path, 'r') as file:
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if sequence:
                    break
            else:
                sequence += line.strip()

    if not sequence:
        return "No sequence found in the FASTA file"

    # headers
    token = "1c816e2026a070e2a4c76b828275e7f41853a370"
    headers={ 'Content-Type': 'application/json', "Authorization": f"Token {token}" }

    # query url
    swiss_model_url = "https://swissmodel.expasy.org/automodel"

    # data
    data = json.dumps({"target_sequences": [sequence], "project_title":"ICA2 project"})
    # create a post
    response = requests.post(swiss_model_url, data = data, headers = headers)

    # check response
    # Obtain the project_id from the response created above
    project_id = response.json()["project_id"]

    # And loop until the project completes
    while True:
        # We wait for some time
        time.sleep(10)

        # Update the status from the server 
        response = requests.get(
            f"https://swissmodel.expasy.org/project/{project_id}/models/summary/", 
            headers={ "Authorization": f"Token {token}" })

        # Update the status
        status = response.json()["status"]

        print('Job status is now', status)

        if status in ["COMPLETED", "FAILED"]:
            break

'''codes for interaction with user'''
def choice_wrap(folder):
    taxonomy = input("Enter the taxonomic group: ")
    protein_family = input("Enter the protein family: ")
    limit = input("Enter the max number of fetched results: ")
    # Initialize sequences as an empty dictionary
    sequences = {}
    selected_species = {}

    # Create directory for the inerested proteins
    folder_path = f'{folder}/{taxonomy}'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    os.makedirs(folder_path, exist_ok=True)

    #species mapping with a limit
    species_protein_map = get_species_protein_ids(taxonomy, protein_family, limit)
    # giving a list of species and corresponding protein ids
    print("\nAvailable species and corresponding protein IDs:")
    for idx, (species, protein_id) in enumerate(species_protein_map.items(), start=1):
        print(f"{idx}. {species} - {protein_id}")

    # giving options for user
    download_choice = input("\nEnter protein id of your choice. Type 'all' for all sequences, or 'select' to choose specific ones: ")
    # download all sequence
    if download_choice.lower() == 'all':
        sequences = fetch_sequences(species_protein_map, taxonomy, folder_path)
        print("You choosed all.")
    # download selected ones
    elif download_choice.lower() == 'select':
        selected_indexes = input("Enter the numbers of species you are interested in (separated by commas): ")
        selected_species = {species: species_protein_map[species] 
                            for idx, species in enumerate(species_protein_map, start=1)
                            if str(idx) in selected_indexes.split(',')}
        sequences = fetch_sequences(selected_species, taxonomy, folder_path)
        print("your choice:")
        print(selected_species)
    # return choices
    return sequences, taxonomy, folder_path

'''MAIN MENU'''
def main_menu():
    print("\nMain Menu:")
    print("0. Draw a Taxonomy Tree")
    print("1. Conservation Analysis")
    print("2. Scan for PROSITE Motifs")
    print("3. Construct Phylogenetic Tree")
    print("4. Structure Prediction with SWISS-MODEL")
    print("5. Exit")
    choice = input("Choose a module by typing a number: ")
    return choice


'''ALL WORKS START HERE'''
def main():
    '''SETTING THE WORKING DIRECTORY AND EXCEPTION CHECK'''
    try:
        set_working_directory()
    except Exception as e:
        print(f"ERROR:can't set working director: {e}")
        sys.exit(1)

    '''CHOOSE FUNCTIONALITY'''
    while True:
        choice = main_menu()

        if choice == '0':
            # input taxonomy
            taxonomy = input("Enter the taxonomic group: ")
            # print a tree with a total number in brackets
            print(get_taxonomic_tree(taxonomy))
            
        '''CONSERVATION ANALYSIS'''
        if choice == '1':
            # input the main info of the protein
            taxonomy = input("Enter the taxonomic group: ")
            protein_family = input("Enter the protein family: ")
            limit = input("Enter the max number of fetched results: ")

            # fetch all the sequence within a limited number
            all_sequences = fetch_sequences_all(taxonomy, protein_family, limit)

            # run emboss conservation analysis
            conservation_data = run_emboss_conservation(taxonomy, protein_family, all_sequences)

            # Conservation analysis
            aligned_sequences_file = f"{taxonomy}/{protein_family}/aligned_sequences.fasta"
            aligned_sequences = parse_fasta_sequences(aligned_sequences_file)
            conservation_scores = calculate_conservation(aligned_sequences)

            # Print conservation scores
            for i, score in enumerate(conservation_scores):
                print(f"Position {i+1}: Conservation Score = {score:.2f}")


        '''SCAN PROTEIN SEQUENCE WITH INTEREST'''
        if choice == '2':
            # download selected sequences in a folder
            selected_species, taxonomy, folder_path = choice_wrap("prosite")
            scan_all_sequences(selected_species, folder_path)
            reports_path = os.path.join(folder_path,"Scan_reports")
            similarity_threshold = 0.8  # threshold for high similarity
            high_similarity_motifs = find_high_similarity_motifs(reports_path, similarity_threshold)

            for motif, similarity in high_similarity_motifs:
                print(f"Motif: {motif}, Similarity: {similarity}")
        
        '''phylogenetic tree'''
        if choice =='3':
            # download selected sequences in a folder
            sequences, taxonomy, folder_path = choice_wrap("phylotree")
            # merge all FASTA files into one
            merged_fasta_file = os.path.join(folder_path, "merged_sequences.fasta")
            merge_fasta_files(folder_path, merged_fasta_file)

            # convert merged FASTA to PHYLIP format
            phylip_file = os.path.join(folder_path, "converted_sequences.phylip")
            convert_fasta_to_phylip(merged_fasta_file, phylip_file)

            # generate the phylogenetic tree using PhyML
            output_tree = os.path.join(folder_path, "phylogenetic_tree")
            run_phyml(phylip_file, output_tree)

            print(f"Phylogenetic tree generated at: {output_tree}")
            
        
        '''Predict Structure'''
        if choice =='4':
            sequences, taxonomy, folder_path = choice_wrap("swiss-model")
            merged_fasta_file = os.path.join(folder_path, "merged_sequences.fasta")
            merge_fasta_files(folder_path,merged_fasta_file)
            print(predict_structure(merged_fasta_file))

        '''EXIT'''
        if choice =='5':
            print("Goodbye!")
            break

if __name__ == '__main__':
    #set working dir as cwd or a configured one
    working_directory = os.getcwd() if len(sys.argv) == 1 else sys.argv[1]
    input_folder = os.path.join(working_directory, "input_files")
    output_folder = os.path.join(working_directory, "output_files")
    main()
