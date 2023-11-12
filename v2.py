from Bio import Entrez

def get_protein_sequences(protein_family, taxonomic_group, email):
    Entrez.email = email
    # Search for the protein family within the specified taxonomic group
    search_handle = Entrez.esearch(db="protein", term=f"{protein_family}[Protein Name] AND {taxonomic_group}[Organism]")
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    # Fetch the sequence data for the obtained IDs
    protein_ids = search_results['IdList']
    fetch_handle = Entrez.efetch(db="protein", id=",".join(protein_ids), rettype="fasta", retmode="text")
    protein_sequences = fetch_handle.read()
    fetch_handle.close()
    
    return protein_sequences

# User inputs
protein_family_input = input("Enter the protein family: ")
taxonomic_group_input = input("Enter the taxonomic group: ")
email_input = 'your@email.com'  # User should provide their actual email here

# Retrieve and print the sequences
sequences = get_protein_sequences(protein_family_input, taxonomic_group_input, email_input)
print(sequences)

