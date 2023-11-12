
from Bio import Entrez

def get_taxonomic_lineage(taxonomic_group, email):
    Entrez.email = email
    # Search for the taxonomic group to get the NCBI Taxonomy ID
    search_handle = Entrez.esearch(db="taxonomy", term=taxonomic_group)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    tax_id = search_results['IdList'][0]
    # Fetch the detailed taxonomic information
    fetch_handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
    tax_data = Entrez.read(fetch_handle)
    fetch_handle.close()
    # Parse the lineage information
    lineage = []
    for taxon in tax_data[0]['LineageEx']:
        lineage.append((taxon['ScientificName'], taxon['TaxId']))
    return lineage

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

# Retrieve and print the taxonomic lineage
lineage = get_taxonomic_lineage(taxonomic_group_input, email_input)
for name, tax_id in lineage:
    print(f"{name} (ID: {tax_id})")

# Retrieve and print the protein sequences
sequences = get_protein_sequences(protein_family_input, taxonomic_group_input, email_input)
print(sequences)

