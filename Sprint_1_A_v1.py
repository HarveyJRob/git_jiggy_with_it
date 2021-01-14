
"""
Add all modules that need to be imported
"""

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import re, json

"""
Uncomment the relevant line below to select whether you would like the
Entrez.email variable to be hard-coded or pass in by the user.
If you select hard-coded, replace the email address with the one you wish to be used.
We should probably add some type of check on the user input?
"""

Entrez.email = "harvey.j.rob@gmail.com"
#Entrez.email = input("Please enter your email address:")

"""
Function which returns the Entrez gene id when provided with a valid HGNC gene symbol
"""

def get_entrez_gene_id(gene_symbol):
    term = '%s[gene] "Homo sapiens"[orgn]' % gene_symbol
    handle = Entrez.esearch(db='gene', term=term, retmode='xml')
    record = Entrez.read(handle)
    handle.close()
    return (record['IdList'][0])

"""
Uncomment if you want to hard code the gene_id to be passed to entrez
"""

#entrez_id = (get_entrez_gene_id('BRCA1'))

"""
OR Uncomment if you want the user to pass in their chosen gene_symbol from the command line
Similar to email address is there some sort of checking or error handling we should do on their input?
"""

gene_symbol = input("Please enter gene of interest (e.g. BRCA1):")
entrez_id = (get_entrez_gene_id(gene_symbol))

"""
Uncomment if you want to give the option of hiding the sequence information
"""

display_seq = input("Would you like to see sequence information (Y/N)?")

if display_seq == "Y" or display_seq == "y":
    display_seq_len = int(input("Enter a sequence length limit (e.g. 20) If blank full sequence will be displayed:"))
else:
    display_seq_len = 0

"""
Function which returns the valid transcript reference sequences when provided with an Entrez gene ID
"""

def get_transcript_list(entrez_id, gene_table=False):
    handle = Entrez.efetch(db="gene", id=entrez_id, rettype="gene_table", retmode="text")
    record = handle.read()
    handle.close()
    if gene_table is True:
        return record
    coding = re.findall(r"NM_\d+.\d+", record)
    noncoding = re.findall(r"NR_\d+.\d+", record)
    tx_list = coding + noncoding
    tx_list = list(dict.fromkeys(tx_list))
    return tx_list

"""
Pass the entrez_id into the get_transcript_list function.
Save the returned list into the variable transcript_list.
"""

transcript_list = get_transcript_list(str(entrez_id))

"""
Uncomment any of the options below if you would like them printed to the terminal for testing purposes
"""

#print("User Email:", Entrez.email)
#print("Gene Symbol:", gene_symbol)
#print("Gene ID:", entrez_id)
#print("Display Seq Len:", display_seq_len)
#print(transcript_list)


"""
Create an empty dictionary called final to be populated by the code below
"""

final = {}

"""
Run the following code for each transcript in transcript_list
"""

for x in transcript_list:
    handle = Entrez.efetch(db="nucleotide", id=x, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "gb")
    handle.close()
    final[x] = {'id': record.id, 'description': record.description}

    my_cdna_obj = record.seq
    my_cdna = str(record.seq[:display_seq_len])
    #my_cdna_short = (my_cdna[:10] + '..') if len(my_cdna) > 10 else my_cdna - to delete
    my_rna_obj = record.seq.transcribe()
    my_rna = str(my_rna_obj.lower()[:display_seq_len])
    #my_rna_short = (my_rna[:10] + '..') if len(my_rna) > 10 else my_rna - to delete
    final[x]['cdna'] = my_cdna
    final[x]['rna'] = my_rna

    final[x]['gene'] = {}
    gene_symbol = record.features[5].qualifiers['gene'][0]
    final[x]['gene']['gene_symbol'] = gene_symbol
    final[x]['type'] = record.features[5].type.lower()

    """
    Run the following code only for coding transcripts
    """

    if record.features[5].type == "CDS":
        hgnc_id = (re.search(r'HGNC:HGNC:\d+', str(record.features[5].qualifiers['db_xref'])))[0].split(':')[-1]
        hgnc_id = 'HGNC:%s' % hgnc_id
        final[x]['gene']['hgnc_id'] = hgnc_id
        entrez_id = (re.search(r'GeneID:\d+', str(record.features[5].qualifiers['db_xref'])))[0].split(':')[-1]
        final[x]['gene']['entrez_id'] = entrez_id
        final[x]['translation'] = {}
        cds_info = re.search(r"\d+:\d+", str(record.features[5].location))[0]
        final[x]['translation']['cds_start'], final[x]['translation']['cds_end'] = cds_info.split(":")
        final[x]['translation']['id'] = record.features[5].qualifiers['protein_id'][0]
        final[x]['translation']['description'] = record.features[5].qualifiers['product'][0]
        final[x]['translation']['prot'] = record.features[5].qualifiers['translation'][0][:display_seq_len]
    else:
        continue

"""
Convert our python dictionary into a JSON string using json.dumps()
Use four indents to format the JSON string and make it easier to read
"""

final_json = json.dumps(final, indent=4)

"""
Print out the final result
"""

print(final_json)


