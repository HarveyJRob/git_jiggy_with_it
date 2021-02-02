
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import re, json, time
import logging, sys, requests
import logging.handlers as handlers

"""
Logging Information 
"""

logging.basicConfig(level=logging.INFO, filename='usage.log', filemode='a', format='%(name)s - %(levelname)s - %(message)s')

logger = logging.getLogger('entrez_request')
logHandler = handlers.RotatingFileHandler('entrez_request.log', maxBytes=500000, backupCount=2)
logHandler.setLevel(logging.ERROR)
logger.addHandler(logHandler)

#logging.debug('this is a debug message')
#logging.info('this is an info message')
#logging.warning('This will get logged to a file')
#logging.error('this is an error error error')
#logging.critical('this is a total critial issue my man')


# Define your custom error
#class MyZeroDivisionError(Exception):
#    pass

# Use an exception block to capture a known error
#try:
#    print(1/0)
#except ZeroDivisionError as e:
    # Once the Exception in caught, we can raise our custom exception with a better value
#    error = 'The warning "%s" is actually trying to tell you that it is illegal to divide something by Zero!' % str(e)
#    raise MyZeroDivisionError(error)



try:
    response = requests.get(input("Enter a URL e.g. https://www.bbc.co.uk: "))
except requests.exceptions.ConnectionError as e:
    logger.exception("Exception occurred", exc_info=True)
    print(-1, 'Connection Error')
else:
    print(response.status_code, response.content)



"""
Function to collection required information from the user information
"""

def request_user_input():
    user_email = input("Please enter your email address:") or "robert.harvey-3@postgrad.manchester.ac.uk"
    gene_symbol = input("Please enter gene of interest (e.g. BRCA1):")
    display_seq = input("Would you like to see sequence information (Y/N)?")

    if display_seq == "Y" or display_seq == "y":
        try:
            display_seq_len = int(input("Enter a sequence length limit (e.g. 20) or leave blank to see full sequence:"))
        except ValueError:
            logger.error(ValueError)
            display_seq_len = 1000000
            print("Full Sequence will be displayed")
    else:
        display_seq_len = 0

    return {'input1':user_email, 'input2':gene_symbol, 'input3':display_seq_len}

"""
Run request_user_input() and save the dictionary output as user_input
Next save the value of each key:value pair to its own variable
"""

user_input = request_user_input()

#Should really do some sort of verification that user has entered an email address
user_email = user_input["input1"]

# Always tell NCBI who you are
Entrez.email = user_email

gene_symbol = user_input["input2"]

display_seq_len = user_input["input3"]

# Always tell NCBI who you are
Entrez.email = user_email


"""
Function which returns the Entrez gene id when provided with a valid HGNC gene symbol
"""

def get_entrez_gene_id(gene_symbol):

    #Create search term
    term = '%s[gene] "Homo sapiens"[orgn]' % gene_symbol

    #Request data from Entrez using Bio.Entrez
    handle = Entrez.esearch(db='gene',
                            term=term,
                            retmode='xml'
                            )

    #The data is returned in xml format
    #Convert it into a python dict using Entrez.read
    record = Entrez.read(handle)

    #Close the handle
    handle.close()

    # Return the gene id which is stored as a list element keyworded 'IdList'
    # There should only ever be 1 gene ID for each gene!
    return (record['IdList'][0])


"""
Function which returns the valid transcript reference sequences when provided with an Entrez gene ID
"""

def get_transcript_list(entrez_id, gene_table=False):

    # Request data from Entrez using Bio.Entrez
    handle = Entrez.efetch(db="gene",
                           id=entrez_id,
                           rettype="gene_table",
                           retmode="text"
                           )

    # We have returned text so that it can be easily manipulated.
    # We read the text as if it were a file
    record = handle.read()

    # Close the handle
    handle.close()

    # Optional request for the gene_table rather than the transcript list
    if gene_table is True:
        return record

    # Search the text for all coding transcripts using re.findall - returns a list
    coding = re.findall(r"NM_\d+.\d+", record)
    # Search the text for all non-coding transcripts using re.findall - returns a list
    noncoding = re.findall(r"NR_\d+.\d+", record)

    # Combine the list
    tx_list = coding + noncoding

    # Remove duplicate items - see * below
    tx_list = list(dict.fromkeys(tx_list))

    # Return the transcript list
    return tx_list


"""
Uncomment any of the options below if you would like them printed to the terminal for testing purposes
"""

print("User Email:", Entrez.email)
print("Gene Symbol:", gene_symbol)
print("Length of sequence to display:", display_seq_len)

"""
Function creates a dictionary for each transcript in transcript list and appends it to final_list
"""

def build_transcript_list(gene_symbol, display_seq_len):

    print("Building transcript list.... please wait")

    # Create an empty list to store the dictionaries created for each transcript
    final_list = []

    display_seq_len = int(display_seq_len)

    # Pass the gene_symbol into the get_entrez_gene_id function
    # Save the result to entrez_id
    entrez_id = (get_entrez_gene_id(gene_symbol))

    # Pass the entrez_id into the get_transcript_list function.
    # Save the returned list into the variable transcript_list.
    transcript_list = get_transcript_list(str(entrez_id))

    # Run the following code for each transcript in transcript_list
    for transcript in transcript_list:

        # Fetch the sequence record from RefSeq.
        # Note db is the database and gb is a Genbank record
        handle = Entrez.efetch(db="nucleotide", id=transcript, rettype="gb", retmode="text")

        # Parse the handle into a SeqRecord object called record
        record = SeqIO.read(handle, "gb")

        # Close the handle
        handle.close()

        # Assign record.seq to an object
        my_cdna_obj = record.seq

        # Extract the cDNA sequence
        # Limit the length of the variable based on display_seq_len
        my_cdna = str(record.seq[:display_seq_len])
        # my_cdna_short = (my_cdna[:10] + '..') if len(my_cdna) > 10 else my_cdna - to delete

        # Re-create the Seq object in RNA format
        my_rna_obj = record.seq.transcribe()

        # Extract the RNA sequence, make it lower case
        # Limit the length of the variable based on display_seq_len
        my_rna = str(my_rna_obj.lower()[:display_seq_len])

        # SeqRecord.features is a list of SeqFeature objects - find the index of the CDS info
        element = 0
        is_cds = False
        for feature in record.features:
            if feature.type == 'CDS':  # Note: a non-coding transcript type == ncRNA
                is_cds = True
                break
            else:
                element = element + 1

        # Extract data into variables
        if is_cds == True:
            my_hgnc_id = \
                (re.search(r'HGNC:HGNC:\d+', str(record.features[element].qualifiers['db_xref'])))[0].split(':')[-1]
            hgnc_id = 'HGNC:%s' % my_hgnc_id
            my_type = str(record.features[element].type)
            my_cds_digits = re.findall(r'\d+', str(record.features[element].location))
            my_cds_start = my_cds_digits[0]
            my_cds_end = my_cds_digits[1]
            my_id = str(record.features[element].qualifiers['protein_id'][0])
            my_translation_description = str(record.features[element].qualifiers['product'][0])
            my_prot = str(record.features[element].qualifiers['translation'][0])

        # Add the data to the dictionary
        t_dict = {'id': record.id, 'description': record.description}
        t_dict['cdna'] = my_cdna
        t_dict['rna'] = my_rna
        t_dict['gene'] = {}
        t_dict['gene']['gene_symbol'] = gene_symbol
        t_dict['gene']['entrez_id'] = entrez_id

        if is_cds == True:
            t_dict['gene']['hgnc_id'] = hgnc_id
            t_dict['type'] = my_type  # Note: non-coding RNAs have the type ncrna (ncRNA)
            t_dict['translation'] = {}
            t_dict['translation']['cds_start'] = my_cds_start  # Note: ncRNA don't have protein information so set to None
            t_dict['translation']['cds_end'] = my_cds_end
            t_dict['translation']['id'] = my_id
            t_dict['translation']['description'] = my_translation_description
            t_dict['translation']['prot'] = my_prot[:display_seq_len]

        # append each dictionary to the list
        final_list.append(t_dict)

    # Convert our python dictionary into a JSON string using json.dumps()
    # Use four indents to format the JSON string and make it easier to read
    final_json = json.dumps(final_list, indent=4, sort_keys=True)

    return final_json


time_start = time.time()

final_list = build_transcript_list(gene_symbol, display_seq_len)

time_end = time.time()

time_diff = time_end - time_start

logger.info("Response time: {}, User email: {}, Gene_Symbol: {}, Display_seq_len: {}"
            .format(time_diff, user_email, gene_symbol, display_seq_len))


"""
Print out the final result
"""

print(final_list)

