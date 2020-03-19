import json
import sys

print(sys.argv)

with open('client/vegaSpecs.json', 'r') as vspec:
    vdata = json.load(vspec)

with open('data.json', 'r') as datafile:
    data = json.load(datafile)

values = {}
values['values'] = data
values['name'] = 'fasta'

vdata['width'] = 700

domain = {}
domain = [sys.argv[1], sys.argv[2]]
vdata['scales'][0]['domain'] = domain

vdata['data'][1] = values
with open('static_vega_lite.json', 'w') as newvspec:
    json.dump(vdata, newvspec, indent=4)

print('created static_vega.json')