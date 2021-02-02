"""
Simple rest interface for VariantValidator built using Flask Flask-RESTPlus and Swagger UI
"""

# Import modules
from flask import Flask
from flask_restplus import Api, Resource
import Sprint_1_A_v3_RH
import requests, sys



# Define the application as a Flask app with the name defined by __name__ (i.e. the name of the current module)
# Most tutorials define application as "app", but I have had issues with this when it comes to deployment,
# so application is recommended
application = Flask(__name__)

# Define the API as api
api = Api(app=application)


# Define a name-space to be read Swagger UI which is built in to Flask-RESTPlus
# The first variable is the path of the namespace the second variable describes the space
hello_space = api.namespace('hello', description='Simple API that returns a greeting')

@hello_space.route("/")
class HelloClass(Resource):
    def get(self):
        return {
            "greeting": "Hello World"
        }


vrecoder_space = api.namespace('VariantRecoder', description='Variant Recoder APIs')
@vrecoder_space.route("/<string:ensembl_variant>")
class VrecoderClass(Resource):
    def get(self, ensembl_variant):
        server = "http://rest.ensembl.org"
        ext1 = "/variant_recoder/human/"
        variant = ensembl_variant

        r = requests.get(server + ext1 + variant + "?", headers={"Content-Type": "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()
        print(repr(decoded))
        return decoded


ensembl_space = api.namespace('Ensembl', description='Ensemble APIs')
@ensembl_space.route("/<string:ensembl_variant>")
class EnsemblClass(Resource):
    def get(self, ensembl_transcript):
        server = "http://rest.ensembl.org"
        ext1 = "/lookup/id/"
        variant = ensembl_transcript

        r = requests.get(server + ext1 + variant + "?", headers={"Content-Type": "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()
        print(repr(decoded))
        return decoded



vv_space = api.namespace('VariantValidator', description='VariantValidator APIs')
@vv_space.route("/variantvalidator/<string:genome_build>/<string:variant_description>/<string:select_transcripts>")
class VariantValidatorClass(Resource):
    def get(self, genome_build, variant_description, select_transcripts):
        # Make a request to the current VariantValidator rest-API
        url = '/'.join(['http://rest.variantvalidator.org/VariantValidator/variantvalidator', genome_build, variant_description,
                        select_transcripts])
        print(url)
        validation = requests.get(url)
        print(validation)
        content = validation.json()
        return content


entrez_space = api.namespace('Entrez', description='Request to Entrez APIs and format output')

@entrez_space.route("/<string:gene_symbol>/<string:display_seq_len>")
class EntrezClass(Resource):
    def get(self, gene_symbol, display_seq_len):
        final_json = Sprint_1_A_v3_RH.build_transcript_list(gene_symbol, display_seq_len)
        return final_json


# Allows app to be run in debug mode
if __name__ == '__main__':
    application.debug = True  # Enable debugging mode
    application.run(host="127.0.0.1", port=8000)  # Specify a host and port fot the app