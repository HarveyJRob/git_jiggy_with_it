
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC
import re, json

Entrez.email = "harvey.j.rob@gmail.com"

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



def create_final_dict(gene_symbol, sequence_display):

    if sequence_display is None:
        display_seq_len = 0
    else:
        display_seq_len = int(sequence_display)

    # get_entrez_gene_id returns entrez_id for given gene.
    entrez_id = (get_entrez_gene_id(gene_symbol))

    # get_transcript_list with the entrez id returned by get_entrez_gene_id
    transcript_list = get_transcript_list(str(entrez_id))

    #Create an empty dictionary called final to be populated by the code below
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
    return json.dumps(final, indent=4)



# call function create_dictionary_for_nonmodel_transcript - uncomment to test
#create_final_dict("BRCA1", "10")
