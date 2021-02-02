"""
Simple rest interface for VariantRecoder, VariantValidator & Entrez APIs
Built using Flask Flask-RESTPlus and Swagger UI
Run the program using $python Sprint_2_API_v1.py
Then open a web browser - http://127.0.0.1:8000/

"""

# Import modules
from flask import Flask
from flask_restplus import Api, Resource
import Sprint_2_Functions_v1
import requests, sys


# Define the application as a Flask app with the name defined by __name__ (i.e. the name of the current module)
# Most tutorials define application as "app", but I have had issues with this when it comes to deployment,
# so application is recommended
application = Flask(__name__)

# Define the API as api
api = Api(app=application)


# Define a name-space to be read Swagger UI which is built in to Flask-RESTPlus
# The first variable is the path of the namespace the second variable describes the space


vrvv_space = api.namespace('VariantRecoder_VariantValidator', description='Variant Recoder into Variant Validator API')
@vrvv_space.route("/<string:ensembl_variant>")
class VRVV_Class(Resource):
    def get(self, ensembl_variant):

        # Protocol and base URL for ensembl API
        vr_server = "http://rest.ensembl.org"

        # Path or route to the namespace we want
        vr_ext = "variant_recoder/human"

        # The Ensembl variant entered by the user in this/our API
        vr_variant = ensembl_variant

        # Create the URL for the Variant Recoder request
        vr_url = '/'.join([vr_server, vr_ext, vr_variant])

        # Variable to store the response of requests.get
        # If we can't to offer XML and JSON then we would have to change the headers
        vr_r = requests.get(vr_url + "?", headers={"Content-Type": "application/json"})

        # Returns an HTTPError object if an error has occurred during the process
        # Raises a SystemExit exception
        if not vr_r.ok:
            vr_r.raise_for_status()
            sys.exit()

        # Save the response to a new variable.
        # I think the .json() converts the JSON object to it's equivalent in python
        # In this case it seems to be a list
        vr_decoded = vr_r.json()
        #print(repr(decoded))

        #Save the value of the "hgvsg" key into a variable to be passed to Variant Validator
        variant_description = vr_decoded[0]["hgvsg"][0]

        # These are other variables the Variant Validator API requires.
        # What's the best way to pass them to manage them?
        genome_build = "GRCh37"
        select_transcripts = "all"

        # Protocol and base URL for Variant Validator API
        vv_server = "http://rest.variantvalidator.org"

        # Path or route to the namespace we want
        vv_ext = "VariantValidator/variantvalidator"

        # Create the URL for the Variant Validator request
        vv_url = '/'.join([vv_server, vv_ext, genome_build, variant_description, select_transcripts])

        # Variable to store the response of requests.get
        # If we can't to offer XML and JSON then we would have to change the headers
        vv_r = requests.get(vv_url + "?", headers={"Content-Type": "application/json"})

        #validation = requests.get(vv_url)

        # Save the response to a new variable.
        # I think the .json() converts the JSON object to it's equivalent in python
        # In this case it seems to be a dictionary
        vv_content = vv_r.json()

        return vv_content


ensembl_lookup_space = api.namespace('Ensembl_Lookup', description='Ensembl Lookup APIs')
@ensembl_lookup_space.route("/<string:genome_build>/<string:variant>")
class Ensembl_Lookup_Class(Resource):
    def get(self, genome_build, variant):

        user_input = {}
        user_input["genome_build"] = genome_build
        user_input["full_input"] = variant

        variant_split_1 = variant.split(":")

        if len(variant_split_1) == 1:
            user_input["seq_id"] = str(variant_split_1[0])
        elif len(variant_split_1) == 2:
            user_input["seq_id"] = str(variant_split_1[0])
            user_input["variant"] = str(variant_split_1[1])
        else:
            return "There is a problem with the variant input"

        variant_split_2 = variant_split_1[0].split(".")

        if len(variant_split_2) == 1:
            user_input["ascession"] = variant_split_2[0]
        elif len(variant_split_2)== 2:
            user_input["ascession"] = variant_split_2[0]
            user_input["version"] = variant_split_2[1]
        else:
            return "There is a problem with the variant input"



        # Protocol and base URL for ensembl API
        if genome_build == "grch37":
            el_server = "https://grch37.rest.ensembl.org"
        else:
            el_server = "http://rest.ensembl.org"

        # Path or route to the namespace we want
        el_ext = "lookup/id"

        # Create the URL for the Variant Recoder request
        el_url = '/'.join([el_server, el_ext, user_input["ascession"]])

        # Variable to store the response of requests.get
        # If we can't to offer XML and JSON then we would have to change the headers
        el_r = requests.get(el_url + "?", headers={"Content-Type": "application/json"})

        # Returns an HTTPError object if an error has occurred during the process
        # Raises a SystemExit exception
        if not el_r.ok:
            el_r.raise_for_status()
            sys.exit()

        # Save the response to a new variable.
        # I think the .json() converts the JSON object to it's equivalent in python
        # In this case it seems to be a list
        el_content = el_r.json()

        combined_dicts = {"user_input": user_input, "el_content": el_content}

        #print(repr(el_content))
        return combined_dicts


vr_space = api.namespace('VariantRecoder', description='Variant Recoder APIs')
@vr_space.route("/<string:ensembl_variant>")
class VR_Class(Resource):
    def get(self, ensembl_variant):

        # Protocol and base URL for ensembl API
        vr_server = "http://rest.ensembl.org"

        # Path or route to the namespace we want
        vr_ext = "variant_recoder/human"

        # The Ensembl variant entered by the user in this/our API
        vr_variant = ensembl_variant

        # Create the URL for the Variant Recoder request
        vr_url = '/'.join([vr_server, vr_ext, vr_variant])

        # Variable to store the response of requests.get
        # If we can't to offer XML and JSON then we would have to change the headers
        vr_r = requests.get(vr_url + "?", headers={"Content-Type": "application/json"})

        # Returns an HTTPError object if an error has occurred during the process
        # Raises a SystemExit exception
        if not vr_r.ok:
            vr_r.raise_for_status()
            sys.exit()

        # Save the response to a new variable.
        # I think the .json() converts the JSON object to it's equivalent in python
        # In this case it seems to be a list
        vr_content = vr_r.json()
        print(repr(vr_content))
        return vr_content


vv_space = api.namespace('VariantValidator', description='VariantValidator APIs')
@vv_space.route("/variantvalidator/<string:genome_build>/<string:variant_description>/<string:select_transcripts>")
class VV_Class(Resource):
    def get(self, genome_build, variant_description, select_transcripts):

        # Protocol and base URL for Variant Validator API
        vv_server = "http://rest.variantvalidator.org"

        # Path or route to the namespace we want
        vv_ext = "VariantValidator/variantvalidator"

        # Create the URL for the Variant Validator request
        vv_url = '/'.join([vv_server, vv_ext, genome_build, variant_description, select_transcripts])
        print(vv_url)

        # Variable to store the response of requests.get
        # If we can't to offer XML and JSON then we would have to change the headers
        vv_r = requests.get(vv_url + "?", headers={"Content-Type": "application/json"})
        #vv_r = requests.get(vv_url)

        # Save the response to a new variable.
        # I think the .json() converts the JSON object to it's equivalent in python
        # In this case it seems to be a dictionary
        vv_content = vv_r.json()

        return vv_content


entrez_space = api.namespace('Entrez', description='Request to Entrez APIs and format output')
@entrez_space.route("/<string:gene_symbol>/<string:display_seq_len>")
class Entrez_Class(Resource):
    def get(self, gene_symbol, display_seq_len):
        final_json = Sprint_2_Functions_v1.build_transcript_list(gene_symbol, display_seq_len)
        return final_json



# Allows app to be run in debug mode
if __name__ == '__main__':
    application.debug = True  # Enable debugging mode
    application.run(host="127.0.0.1", port=8000)  # Specify a host and port fot the app