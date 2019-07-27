import json
import sys

print(sys.argv)

with open('client/vlSpec.json', 'r') as vlspec:
    vldata = json.load(vlspec)

with open('data.json', 'r') as datafile:
    data = json.load(datafile)

values = {}
values['values'] = data

vldata['width'] = 700

domain = {}
domain['domain'] = [sys.argv[1], sys.argv[2]]
vldata['encoding']['x']['scale'] = domain

vldata['data'] = values
with open('static_vega_lite.json', 'w') as newvlspec:
    json.dump(vldata, newvlspec, indent=4)

print('created static_vega_lite.json')