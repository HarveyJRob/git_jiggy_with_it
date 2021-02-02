"""
Simple rest interface that interacts with Ensembl, VariantValidator & Entrez APIs
Built using Flask Flask-RESTPlus and Swagger UI
Run the program using $python Sprint_2_API_v2.py
Then open a web browser - http://127.0.0.1:8000/
"""

# Import modules
from flask import Flask, make_response, request
from flask_restplus import Api, Resource, reqparse
import requests
import sys
from requests.exceptions import ConnectionError
from dicttoxml import dicttoxml

import logging
import logging.handlers as handlers
import time


"""
Logging Information
"""

logging.basicConfig(level=logging.ERROR,
                    filename='usage.log',
                    filemode='a',
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p')

logger = logging.getLogger()

"""
logger = logging.getLogger(__name___)
logger.setLevel(logging.INFO)


logHandler = handlers.RotatingFileHandler('sprint2.log',
                                          maxBytes=500000,
                                          backupCount=2)
logHandler.setLevel(logging.INFO)

logger.addHandler(logHandler)

"""

# Define the application as a Flask app with the name defined by __name__ (i.e. the name of the current module)
# Most tutorials define application as "app", but I have had issues with this when it comes to deployment,
# so application is recommended
application = Flask(__name__)

# By default, show all endpoints (collapsed)
application.config.SWAGGER_UI_DOC_EXPANSION = 'list'

# Define the API as api
api = Api(app=application,
          version="2.0",
          title="Rest_API_Sprint2",
          description="### REST API for Sprint 2 of Intro to Programming<br>"
                      "Interacts with Ensembl, VariantValidator & Entrez APIs")



# Create a RequestParser object to identify specific content-type requests in HTTP URLs
# The requestparser allows us to specify arguements passed via a URL, in this case, ....?content-type=application/json
parser = reqparse.RequestParser()
parser.add_argument('content-type',
                    type=str,
                    help='***Select the response format***',
                    choices=['application/json', 'application/xml'])


"""
Register custom exceptions
class RemoteConnectionError(Exception):
    code=504
"""

"""
Representations
 - Adds a response-type into the "Response content type" drop-down menu displayed in Swagger
 - When selected, the APP will return the correct response-header and content type
 - The default for flask-restplus is aspplication/json
"""

# Add additional representations using the @api.representation decorator
# Requires the module make_response from flask and dicttoxml
@api.representation('application/xml')
def xml(data, code, headers):
    data = dicttoxml(data)
    resp = make_response(data, code)
    resp.headers['Content-Type'] = 'application/xml'
    return resp

@api.representation('application/json')
def json(data, code, headers):
    resp = make_response(data, code)
    resp.headers['Content-Type'] = 'application/json'
    return resp


"""
A function to process the user input
Some basic error checking is done (this needs to be improved)
Data is stored and returned in a python dict called user_input
"""

def crunch_input(genome_build, variant, output_format):

    # Create an empty dict and two empty lists
    user_input = {}
    user_input["input_error"] = []
    variant_split_1 = []
    variant_split_2 = []

    # Force output_format to be lower case
    # Check which output_format the user entered and save it to user_input
    if str(output_format).lower() == "c":
        user_input["output_format"] = "compact"
    elif str(output_format).lower() == "f":
        user_input["output_format"] = "full"

    # Any problems save an input_error message
    else:
        user_input["input_error"].append("PLEASE ENTER A VALID OUTPUT FORMAT")

    # Force genome_build to be lower case
    # Check which genome build the user entered and save it to user_input
    if str(genome_build).lower() == "grch37" or str(genome_build).lower() == "grch38":
        user_input["genome_build"] = str(genome_build).lower()

    # Any problems save an input_error message
    else:
        user_input["input_error"].append("PLEASE ENTER A VALID GENOME BUILD")

    # Save un-processed user_input for reference
    user_input["variant_unchecked"] = str(variant)

    # Force variant to be upper case - not sure why but I got errors forcing it to be lower case
    # Check if variant string entered by the user contains a ":"
    # If so split the string and save it to a new list
    variant_split_1 = str(variant).split(":")

    # If variant string doesn't contain a ":" save an error message
    if len(variant_split_1) == 1:
            user_input["input_error"].append("PLEASE ENTER A VARIANT NOT JUST A TRANSCRIPT")

    # If variant string contain more than one ":" save an error message
    if len(variant_split_1) >= 3:
            user_input["input_error"].append("THE VARIANT CONTAINS TOO MANY -> : ")

    # If the variant string contains just one ":"
    if len(variant_split_1) == 2:

        transcript_bit = str(variant_split_1[0]).upper()
        variant_bit = str(variant_split_1[1])

        # Save the string after ":" to user_input
        user_input["variant"] = variant_bit

        # If the string before ":" begins with "ENST" save it to user_input
        if transcript_bit[0:4] == "ENST":
            user_input["ensembl_prefix"] = "ENST"

        # If not save an error message
        else:
            user_input["input_error"].append("INVALID ENSEMBL PREFIX")

        variant_split_2 = transcript_bit.split(".")
        stable_id = str(variant_split_2[0])

        if len(variant_split_2) == 1:

            if len(stable_id) == 15 and stable_id[4:].isdecimal():
                user_input["ascession_num"] = stable_id[4:]
                user_input["version"] = "NOT ENTERED"
            else:
                user_input["input_error"].append("INVALID ASCESSION NUMBER")

        if len(variant_split_2) >= 3:
            user_input["input_error"].append("THE STABLE ID CONTAINS TOO MANY -> . ")

        if len(variant_split_2) == 2:

            version_bit = str(variant_split_2[1])

            if len(stable_id) == 15 and stable_id[4:].isdecimal():
                user_input["ascession_num"] = stable_id[4:]
            else:
                user_input["input_error"].append("INVALID ASCESSION NUMBER")

            if version_bit.isdecimal() and len(version_bit) < 3:
                user_input["version"] = version_bit
            else:
                user_input["input_error"].append("INVALID VERSION NUMBER")

    if user_input["input_error"] != []:
        return user_input

    else:
        user_input["stable_id_no_version"] = str(user_input["ensembl_prefix"]) + str(user_input["ascession_num"])

        if user_input["version"] != "NOT ENTERED":
            user_input["stable_id_with_version"] = user_input["ensembl_prefix"] + user_input["ascession_num"] + \
                                               user_input["version"]
        else:
            user_input["stable_id_with_version"] = user_input["stable_id_no_version"]

        user_input["variant_checked"] = user_input["stable_id_with_version"] + ":" + user_input["variant"]

    # Return the user_input dict
    return user_input


"""
A function to select/format what is returned to the user
Useful whilst testing so we can try out different output formats
"""

def crunch_output(user_input, lookup_content, recoder_content, vv_content):

    if user_input["output_format"] == "compact":

        # The information requested by Pete
        compact_output = {"Variant you provided": user_input["variant_unchecked"],
                          "Genome Build you selected": user_input["genome_build"],
                          "Corrected Variant": recoder_content[0]["input"],
                          "HGVS Genomic": recoder_content[0]["hgvsg"][0]}

        return compact_output

    else:

        # All the info returned from the various APIs. Useful to learn what each one is doing
        full_output = {"user_input": user_input,
                       "lookup_api_content": lookup_content,
                       "recoder_api_content": recoder_content,
                       "vv_api_content": vv_content}

        return full_output


"""
A function to interact with the various ensembl rest APIs
"""

def ensembl_api(build, ext, user_input):

    # Switch the server address based on the genome build
    if build == "grch37":
        server = "https://grch37.rest.ensembl.org"
    else:
        server = "http://rest.ensembl.org"

    # Create the URL for the Variant Recoder request
    url = '/'.join([server, ext, user_input])

    # Variable to store the response of requests.get
    # If we want to offer XML and JSON then we would have to change the headers
    r = requests.get(url + "?", headers={"Content-Type": "application/json"})

    # Returns an HTTPError object if an error has occurred during the process
    # Raises a SystemExit exception
    #if not r.ok:
    #    r.raise_for_status()
    #    sys.exit()

    # Save the response to a new variable.
    # I think the .json() converts the JSON object to it's equivalent in python
    ensembl_content = r.json()

    return ensembl_content


"""
A function to interact with the variant validator rest api
"""

def variantvalidator_api(build, user_input, select_transcripts):

    # Protocol and base URL for Variant Validator API
    server = "http://rest.variantvalidator.org"

    # Path or route to the namespace we want
    ext = "VariantValidator/variantvalidator"

    # Create the URL for the Variant Recoder request
    url = '/'.join([server, ext, build, user_input, select_transcripts])
    # Variable to store the response of requests.get
    # If we want to offer XML and JSON then we would have to change the headers
    #try:
    r = requests.get(url + "?", headers={"Content-Type": "application/json"})
    #except ConnectionError:
    #    raise RemoteConnectionError("Variant Validator server currently unavailable")

    # Save the response to a new variable.
    # I think the .json() converts the JSON object to it's equivalent in python
    vv_content = r.json()

    return vv_content


"""
Define a name-space to be read by Swagger UI which is built in to Flask-RESTPlus
"""

sprint2_space = api.namespace('Sprint2', description='Simple API for Sprint 2')

@sprint2_space.route("/<string:genome_build>/<string:variant>/<string:output_format>")

@sprint2_space.param("variant", "***Accepted Format - case insensitive:***\n"
                                       ">   ENST00000366667.4:c.803T>C - ascession.version:variant\n"
                                       ">   ENST00000366667:c.803T>C - ascession:variant\n")

@sprint2_space.param("genome_build", "***Accepted Values - case insensitive:***\n"
                                ">   GRCh37\n"
                                ">   GRCh38\n")

@sprint2_space.param("output_format", "***Accepted Values - case insensitive:***\n"
                                ">   C = Compact (Key info only)\n"
                                ">   F = Full (All info)\n")


class Sprint2_Class(Resource):
    # Add documentation about the parser
    @api.expect(parser, validate=True)
    def get(self, genome_build, variant, output_format):

        # start a timer
        time_start = time.time()

        # Run the crunch_input() to check and store all input provided by the user
        user_input = crunch_input(genome_build, variant, output_format)

        # Checks to see if the input_error key is an empty list (no errors)
        # If so, this function returns the user_input and exits
        if user_input["input_error"] != []:
            time_end = time.time()
            time_diff = time_end - time_start
            user_input["query_time"] = time_diff
            logger.warning("Response time: {}, Genome Build: {}, Variant: {}, Output Format: {}"
                        .format(time_diff, genome_build, variant, output_format))
            return user_input

        # Run ensemble_api() hitting the lookup/id endpoint for the relevant genome build.
        # The final argument is the stable id w/o a version number extracted from crunching the users input.
        lookup_content = ensembl_api(user_input["genome_build"], "lookup/id", user_input["stable_id_no_version"])

        # Checks to see if a key:value pair for error exists in lookup_content
        # This would have been created if stable ID couldn't be found
        # If so this function returns the user_input and exits
        if "error" in lookup_content.keys():
            user_input["LookupID_error"] = lookup_content["error"]

        # Run ensembl_api() hitting the variant_recoder endpoint for the relevant genome build.
        # The final argument is a concat of id & version from lookup_content and variant from user_input.
        # It might be more readable to create a variable to store this info and pass the variable into the function?

        recoder_content = {}

        try:
            recoder_content = ensembl_api(user_input["genome_build"],
                                          "variant_recoder/human",
                                          str(lookup_content["id"]) + "." +
                                          str(lookup_content["version"]) + ":" +
                                          str(user_input["variant"]))
        except KeyError as RecoderErr:
            logger.error(LookupError)
            user_input["Recoder_error"] = "Missing ID or Version from lookup content"

        # Run variantvalidator_api()
        # The genome_build is the one provided by the user
        # The variant is the HGVS Genomic info returned from variant recoder
        # At the moment the select_transcripts variable is set to all

        vv_content = {}

        try:
            vv_content = variantvalidator_api(user_input["genome_build"],
                                              recoder_content[0]["hgvsg"][0],
                                              "all")
        except KeyError as VVErr:
            logger.error(VVErr)
            user_input["VV_error"] = "Missing HGVSG from recoder content"

        # Send all dicts to the crunch_output function to format based on the users preference
        crunched_output = crunch_output(user_input, lookup_content, recoder_content, vv_content)

        # Stop the timer
        time_end = time.time()

        # Calculate how long the timer ran for
        time_diff = time_end - time_start

        # Write a log message for the user query and how long it too to complete
        logger.info("Response time: {}, Genome Build: {}, Variant: {}, Output Format: {}"
                    .format(time_diff, genome_build, variant, output_format))

        # If the programme has got this far, return crunched output
        return crunched_output


"""
        # Collect Arguements
        args = parser.parse_args()

        # Overides the default response route so that the standard HTML URL can return any specified format
        if args['content-type'] == 'application/json':
            # example: http://127.0.0.1:5000.....bob?content-type=application/json
            return json(crunched_output, 200, None)
        # example: http://127.0.0.1:5000.....?content-type=application/xml
        elif args['content-type'] == 'application/xml':
            return xml(crunched_output, 200, None)
        else:
            # Return the api default output
            return crunched_output
"""


"""
Error handlers

# Simple function that creates an error message that we will log
def log_exception(type):
    # We want to know the arguments passed and the path so we can replicate the error
    params = dict(request.args)
    params['path'] = request.path
    # Create the message and log
    message = '%s occurred at %s with params=%s' % (type, time.ctime(), params)
    logger.exception(message, exc_info=True)


@application.errorhandler(RemoteConnectionError)
def remote_connection_error_handler(e):
    # Add the Exception to the log ensuring that exc_info is True so that a traceback is also logged
    log_exception('RemoteConnectionError')

    # Collect Arguments
    args = parser.parse_args()
    if args['content-type'] != 'application/xml':
        return json({'message': str(e)},
                                504,
                                None)
    else:
        return xml({'message': str(e)},
                   504,
                   None)


@application.errorhandler(404)
def not_found_error_handler():
    # Collect Arguments
    args = parser.parse_args()
    if args['content-type'] != 'application/xml':
        return json({'message': 'Requested Endpoint not found'},
                                404,
                                None)
    else:
        return xml({'message': 'Requested Endpoint not found'},
                   404,
                   None)


@application.errorhandler(500)
def default_error_handler():
    # Add the Exception to the log ensuring that exc_info is True so that a traceback is also logged
    log_exception('RemoteConnectionError')

    # Collect Arguments
    args = parser.parse_args()
    if args['content-type'] != 'application/xml':
        return json({'message': 'unhandled error: contact Sprint 2 systems Administrator'},
                                500,
                                None)
    else:
        return xml({'message': 'unhandled error: contact Sprint 2 Systems Administrator'},
                   500,
                   None)
"""


# Allows app to be run in debug mode
if __name__ == '__main__':
    application.debug = True  # Enable debugging mode
    application.run(host="127.0.0.1", port=8000)  # Specify a host and port fot the app