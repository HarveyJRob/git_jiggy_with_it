"""
Simple rest interface that interacts with Ensembl, VariantValidator & Entrez APIs
Built using Flask Flask-RESTPlus and Swagger UI
Run the program using $python Sprint_2_API_v2.py
Then open a web browser - http://127.0.0.1:8000/
"""

# Import modules
from flask import Flask
from flask_restplus import Api, Resource
import requests, sys


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


"""
A function to process the user input
Some basic error checking is done (this needs to be improved)
Data is stored and returned in a python dict called user_input
"""

def crunch_input(genome_build, variant, output_format):

    # Create an empty user_input dict
    user_input = {}

    if str(output_format).lower() == "c":
        user_input["output_format"] = "compact"
    elif str(output_format).lower() == "f":
        user_input["output_format"] = "full"
    else:
        user_input["input_error"] = "PLEASE ENTER A VALID OUTPUT FORMAT"

    # Force genome_build to be lower case
    # Check which genome build the user entered and save it to user_input
    if str(genome_build).lower() == "grch37" or str(genome_build).lower() == "grch38":
        user_input["genome_build"] = str(genome_build).lower()

    # Any problems save an input_error message
    else:
        user_input["input_error"] = "PLEASE ENTER A VALID GENOME BUILD"

    # Save un-processed user_input for reference
    user_input["full_user_input"] = str(variant)

    # Force variant to be upper case - not sure why but I got errors forcing it to be lower case
    # Check if variant string entered by the user contains a ":"
    # If so split the string and save it to a new list
    variant_split_1 = str(variant).upper().split(":")

    # If variant string doesn't contain a ":" save an error message
    if len(variant_split_1) == 1:
        user_input["input_error"] = "PLEASE ENTER A VARIANT NOT JUST A TRANSCRIPT"

    # Else if string contains ":" check it begins with "ENST", is 15 chars long or contains a "." at char 16
    # If so save everything before ":" to ascession key and after ":" to variant key
    elif len(variant_split_1) == 2 \
            and (len(variant_split_1[0]) == 15 or str(variant_split_1[0][15] == ".")) \
            and str(variant_split_1[0][:4]) == "ENST":
        user_input["ascession"] = str(variant_split_1[0])
        user_input["variant"] = str(variant_split_1[1])

    # Otherwise save an input_error message
    else:
        user_input["input_error"] = "THERE IS A PROBLEM WITH THE VARIANT FORMAT (1)"

    # Check if the first item in the list variant_split_1 contains a :.:.
    # If so split the string and save it to a new list
    variant_split_2 = variant_split_1[0].split(".")

    # If variant string doesn't contain a ".", is 15 chars long and begins with "ENST"
    # Save it to ascession key and set version key to "NOT PROVIDED"
    if len(variant_split_2) == 1 \
            and len(variant_split_2[0]) == 15 \
            and str(variant_split_2[0][:4]) == "ENST":
        user_input["ascession"] = variant_split_2[0]
        user_input["version"] = "NOT PROVIDED"

    # If variant contains a "." is 15 chars long and begins with "ENST"
    # Save everything before the "." to ascession and after the "." to version
    elif len(variant_split_2) == 2 \
            and len(variant_split_2[0]) == 15 \
            and str(variant_split_2[0][:4]) == "ENST":
        user_input["ascession"] = variant_split_2[0]
        user_input["version"] = variant_split_2[1]

    # Otherwise save an input_error message
    else:
        user_input["input_error"] = "THERE IS A PROBLEM WITH THE VARIANT FORMAT (2)"

    # Return the user_input dict
    return user_input


"""
A function to select/format what is returned to the user
Useful whilst testing so we can try out different output formats
"""

def crunch_output(user_input, lookup_content, recoder_content, vv_content):

    if user_input["output_format"] == "compact":

        # The information requested by Pete
        compact_output = {"Variant you provided": user_input["full_user_input"],
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
    if not r.ok:
        r.raise_for_status()
        sys.exit()

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
    r = requests.get(url + "?", headers={"Content-Type": "application/json"})

    # Returns an HTTPError object if an error has occurred during the process
    # Raises a SystemExit exception
    if not r.ok:
        r.raise_for_status()
        sys.exit()

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
    def get(self, genome_build, variant, output_format):

        # Run the crunch_input() to check and store all input provided by the user
        user_input = crunch_input(genome_build, variant, output_format)


        # Checks to see if a key:value pair for input_error exists
        # This would have been created whilst crunching the user input
        # If so this function returns the user_input and exits
        if "input_error" in user_input.keys():
            return user_input


        # Run ensemble_api() hitting the lookup/id endpoint for the relevant genome build.
        # The final argument is the ascession number extracted from crunching the users input.
        lookup_content = ensembl_api(user_input["genome_build"], "lookup/id", user_input["ascession"])


        """
        TO DELETE
        This code can be used to run lookup/id even if a variant isn't entered.
        At the moment program would exits before reaching this code if a variant isn't entered
        
        if "variant" not in user_input.keys():
            combined_dicts = {"user_input": user_input,
                              "lookup_api_content": lookup_content}
            return combined_dicts
        """


        # Run ensembl_api() hitting the variant_recoder endpoint for the relevant genome build.
        # The final argument is a concat of id & version from lookup_content and variant from user_input.
        # It might be more readable to create a variable to store this info and pass the variable into the function?
        recoder_content = ensembl_api(user_input["genome_build"],
                                      "variant_recoder/human",
                                      str(lookup_content["id"]) + "." +
                                      str(lookup_content["version"]) + ":" +
                                      str(user_input["variant"]))


        # Run variantvalidator_api()
        # The genome_build is the one provided by the user
        # The variant is the HGVS Genomic info returned from variant recoder
        # At the moment the select_transcripts variable is set to all
        vv_content = variantvalidator_api(user_input["genome_build"], recoder_content[0]["hgvsg"][0], "all")


        # Send all dicts to the crunch_output function to format based on the users preference
        crunched_output = crunch_output(user_input, lookup_content, recoder_content, vv_content)

        # return to the user the dictionary created in the step above
        return crunched_output


# Allows app to be run in debug mode
if __name__ == '__main__':
    application.debug = True  # Enable debugging mode
    application.run(host="127.0.0.1", port=8000)  # Specify a host and port fot the app