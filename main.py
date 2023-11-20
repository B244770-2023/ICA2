#!/usr/bin/python3

import os
import io
import sys
import requests
import subprocess
from Bio.ExPASy import ScanProsite
from Bio.Blast import NCBIWWW
from bs4 import BeautifulSoup

'''SET WORKING DIRECTORY'''
def set_working_directory():
    script_path = os.path.abspath(__file__)
    
    #getting current directory of the script
    script_dir = os.path.dirname(script_path)

    #change working dir
    os.chdir(script_dir)
    print(f"working directory changed to: {script_dir}")

'''RUN COMMANDS FUNCTION'''
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
        
'''DRAW A TAXONOMIC TREE'''
def get_taxonomic_tree(taxonomic_group):
    # Search for the taxonomic group to get the NCBI Taxonomy ID
    search_xml = run_command(f"esearch -db taxonomy -query '{taxonomic_group}' | efetch -format xml")
    root = ET.fromstring(search_xml)
    tax_id = root.find(".//Id").text

    # get the lineage
    lineage_xml = run_command(f"efetch -db taxonomy -id {tax_id} -format xml")
    lineage_root = ET.fromstring(lineage_xml)

    # extract lineage information
    lineage_ids = [taxon.find('TaxId').text for taxon in lineage_root.findall(".//LineageEx/Taxon")]
    lineage_names = [taxon.find('ScientificName').text for taxon in lineage_root.findall(".//LineageEx/Taxon")]

    # generate a tree view as a string
    tree_str = ""
    for depth, (taxid, name) in enumerate(zip(lineage_ids, lineage_names)):
        tree_str += " " * depth + f"- {name} ('{taxid}')\n"
    
    return tree_str

'''FETCH ALL SEQUENCES AND WRITE IN ON FILE, WITH A LIMIT'''
def fetch_sequences_all(taxonomy, protein_family, limit):
    query = f"{protein_family}[Protein Name] AND {taxonomy}[Organism]"
    protein_sequences = run_command(["esearch", "-db", "protein", "-query", query, "|", "efetch", "-format", "fasta"])

    # save sequences to a file
    folder_path = f'{taxonomy}/{protein_family}'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    with open(f'{folder_path}/sequences.fasta', 'w') as file:
        file.write(protein_sequences)

'''GET SPECIES~PROTEIN_IDS LIST'''
def get_species_protein_ids(taxonomy, protein_family, limit):
    # fetch all sequence IDs
    query = f"{protein_family}[Protein Name] AND {taxonomy}[Organism]"
    ids = run_command(f"esearch -db protein -query '{query}' -retmax {limit} | efetch -format docsum")

    # parse XML to extract protein IDs
    root = ET.fromstring(ids)
    protein_ids = [docsum.find('Id').text for docsum in root.findall('DocumentSummary')]

    species_protein_map = {}
    # fetch names of each protein's species
    for protein_id in protein_ids:
        xml_data = run_command(f"efetch -db protein -id {protein_id} -format gb")
        xml_root = ET.fromstring(xml_data)
        species_name = xml_root.find(".//GBSeq_organism").text
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
def fetch_sequences(species_protein_map, taxonomy):
    
    # Create directory for the inerested proteins
    folder_path = f'interested/{taxonomy}'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    sequences = {}
    # fetch handle
    for species, protein_id in species_protein_map.items():
        # Fetch the protein sequence using efetch
        protein_sequence = run_command(["efetch", "-db", "protein", "-id", protein_id, "-format", "fasta"])

        # Write the sequence to a file
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

'''FUNCTION FOR PROSITE SCANNING'''
def scan_prosite(sequences):
    #scan sequences with motifs from prosite db
    matched_motifs = {}
    for sequence in sequences:
        if sequence: # empty check
            search_handle = ScanProsite.scan(seq=sequence, output='html')
            scan_results = ScanProsite.read(search_handle)
            motifs_results[species] = [result['signature_ac'] for result in scan_results]
    return matched_motifs

'''BLAST'''
def do_blast(sequences):
    for sequence in sequences:
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence.seq)
    return result_handle.read()

'''SWISS-MODEL'''
def predict_structure(protein_sequence):
    # query url
    swiss_model_url = "https://swissmodel.expasy.org/interactive"

    # create a post
    response = requests.post(swiss_model_url, data={"seq": protein_sequence})

    # check response
    if response.status_code == 200:
        # get response url using soup
        soup = BeautifulSoup(response.text, "html.parser")
        result_link = soup.find("a", class_="models-result-link")
        if result_link:
            result_url = result_link.get("href")
            return result_url
    else:
        return "Failed to request SWISS-MODEL API"

'''FETCH CODE WRAP'''
'''codes for fetch prefered protein sequences and return in dictionary'''
def fetch_wrap():
    taxonomy = 'aves'
    # taxonomy = input("Enter the taxonomic group: ")
    protein_family = 'glucose-6-phosphatase'
    # protein_family = input("Enter the protein family: ")
    limit = 5
    # limit = input("Enter the max number of fetched results: ")
    # Initialize sequences as an empty dictionary
    sequences = {}
    #species mapping
    species_protein_map = get_species_protein_ids(taxonomy, protein_family, limit)
    # giving a list of species and corresponding protein ids
    print("\nAvailable species and corresponding protein IDs:")
    for idx, (species, protein_id) in enumerate(species_protein_map.items(), start=1):
        print(f"{idx}. {species} - {protein_id}")

    # giving options for user
    download_choice = input("\nEnter protein id of your choice. Type 'all' for all sequences, or 'select' to choose specific ones: ")
    # download all sequence
    if download_choice.lower() == 'all':
        sequences = fetch_sequences(species_protein_map, taxonomy)
    # download selected ones
    elif download_choice.lower() == 'select':
        selected_indexes = input("Enter the numbers of species you are interested in (separated by commas): ")
        selected_species = {species: species_protein_map[species] 
                            for idx, species in enumerate(species_protein_map, start=1)
                            if str(idx) in selected_indexes.split(',')}
        # fetch the sequence
        sequences = fetch_sequences(selected_species, taxonomy)
    return sequences

'''MAIN MENU'''
def main_menu():
    print("0. Draw a Taxonomy Tree")
    print("1. Conservation Analysis")
    print("2. Scan for PROSITE Motifs")
    print("3. Perform Blast Analysis")
    print("4. Structure Prediction with SWISS-MODEL")
    print("5. Exit")
    print("6. Exit")
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
            # input the main info of the protein
            taxonomy = input("Enter the taxonomic group: ")
            protein_family = input("Enter the protein family: ")
            limit = input("Enter the max number of fetched results: ")

            get_taxonomic_tree(taxonomy)
            
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

        '''SCAN PROTEIN SEQUENCE WITH INTEREST'''
        if choice == '2':
            sequences = fetch_wrap()
            # scan motifs
            matched_motifs = scan_prosite(to_fasta_format(sequences))
            # display
            for species, motifs in matched_motifs.items():
                print(f"\n{species} has the following motifs:")
                for motif in motifs:
                    print(motif)
        
        '''BLAST'''
        if choice =='3':
            do_blast(fetch_wrap())
        
        '''Predict Structure'''
        if choice =='4':
            print(predict_structure(fetch_wrap()))

        '''BLAST'''
        if choice =='5':
            do_blast(fetch_wrap)

        '''BLAST'''
        if choice =='6':
            do_blast(fetch_wrap)

if __name__ == '__main__':
    #set working dir as cwd or a configured one
    working_directory = os.getcwd() if len(sys.argv) == 1 else sys.argv[1]
    input_folder = os.path.join(working_directory, "input_files")
    output_folder = os.path.join(working_directory, "output_files")
    main()
