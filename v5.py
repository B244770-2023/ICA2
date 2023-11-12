from Bio import Entrez
from ete3 import Tree, NCBITaxa

def search_taxon(taxonomic_group, email):
    Entrez.email = email
    handle = Entrez.esearch(db='taxonomy', term=taxonomic_group)
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']

def build_tree(taxon_id):
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxon_id)
    names = ncbi.get_taxid_translator(lineage)
    tree_structure = ncbi.get_topology(lineage)
    return tree_structure

def visualize_tree(tree_structure):
    # Here you can use the ETE Toolkit's rendering capabilities
    print(tree_structure.get_ascii(attributes=["sci_name", "rank"]))

# Main code
user_email = input("Enter your email: ")
taxonomic_group_input = input("Enter the taxonomic group: ")
taxon_ids = search_taxon(taxonomic_group_input, user_email)

for taxon_id in taxon_ids:
    tree = build_tree(taxon_id)
    visualize_tree(tree)

