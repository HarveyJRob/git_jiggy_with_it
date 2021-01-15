from urllib.parse import quote

import requests
import json

raw_txt = input("Type your name: ")
quote_txt = quote(raw_txt)

base_url = 'http://127.0.0.1:8000/'


def make_request(base_url, api_function):
    # Tell the User the full URL of their call to the rest API
    url = '%s%s' % (base_url, api_function)
    print("Querying rest API with URL: " + url)

    # Make the request and pass to a response object that the function returns
    response = requests.get(url)
    return response

response = make_request(base_url, 'name/{0}'.format(quote_txt))

print(response.status_code)

print(response.headers)

body = response.json()

print(body)
