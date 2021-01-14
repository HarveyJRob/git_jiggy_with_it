"""
Simple rest interface for VariantValidator built using Flask Flask-RESTPlus and Swagger UI
"""

# Import modules
from flask import Flask
from flask_restplus import Api, Resource
import Sprint_1_A_v2

# Define the application as a Flask app with the name defined by __name__ (i.e. the name of the current module)
# Most tutorials define application as "app", but I have had issues with this when it comes to deployment,
# so application is recommended
application = Flask(__name__)
#application.config['JSONIFY_PRETTYPRINT_REGULAR'] = True

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

name_space = api.namespace('name', description='Return a name provided by the user')

@name_space.route("/<string:name>")
class NameClass(Resource):
    def get(self, name):
        return {
            "My name is": name
        }

vv_space = api.namespace('VariantValidator', description='VariantValidator APIs')

@vv_space.route("/variantvalidator/<string:genome_build>/<string:variant_description>/<string:select_transcripts>")
class VariantValidatorClass(Resource):
    def get(self, genome_build, variant_description, select_transcripts):
        # Make a request to the current VariantValidator rest-API
        url = '/'.join(['http://rest.variantvalidator.org/variantvalidator', genome_build, variant_description,
                        select_transcripts])
        print(url)
        validation = requests.get(url)
        print(validation)
        content = validation.json()
        return content

entrez_space = api.namespace('Entrez', description='Request to Entrez APIs and format output')

@entrez_space.route("/<string:gene_symbol>/<string:sequence_display>")
class EntrezClass(Resource):
    def get(self, gene_symbol, sequence_display):
        content = Sprint_1_A_v2.create_final_dict(gene_symbol, sequence_display)
        return content


# Allows app to be run in debug mode
if __name__ == '__main__':
    application.debug = True  # Enable debugging mode
    application.run(host="127.0.0.1", port=5000)  # Specify a host and port fot the app