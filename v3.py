from Bio import Entrez
from ete3 import NCBITaxa

def get_taxonomic_tree(taxonomic_group, email):
    # Initialize the NCBITaxa
    ncbi = NCBITaxa()
    # Search for the taxonomic group to get the NCBI Taxonomy ID
    Entrez.email = email
    search_handle = Entrez.esearch(db="taxonomy", term=taxonomic_group)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    tax_id_list = search_results['IdList']
    # Retrieve the taxonomic tree for the first ID in the list
    tree = ncbi.get_descendant_taxa(tax_id_list[0], collapse_subspecies=False, return_tree=True)
    return tree

def get_protein_sequences(protein_family, taxonomic_group, email):
    Entrez.email = email
    query = f"{protein_family}[Protein Name] AND {taxonomic_group}[Organism]"
    search_handle = Entrez.esearch(db="protein", term=query, retmax=10)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    protein_ids = search_results['IdList']
    fetch_handle = Entrez.efetch(db="protein", id=",".join(protein_ids), rettype="fasta", retmode="text")
    protein_sequences = fetch_handle.read()
    fetch_handle.close()
    return protein_sequences

# User inputs
email_input = input("Please enter your email: ")
protein_family_input = input("Enter the protein family: ")
taxonomic_group_input = input("Enter the taxonomic group: ")

# Retrieve and print the taxonomic tree
tree = get_taxonomic_tree(taxonomic_group_input, email_input)
print(tree)  # This will print the tree object, but you might want to visualize it using tree.render() or similar.

# Retrieve and print the protein sequences
sequences = get_protein_sequences(protein_family_input, taxonomic_group_input, email_input)
print(sequences)

