import json
import sys

with open('static/vegaSpecs.json', 'r') as vspec:
    vdata = json.load(vspec)

with open('data.json', 'r') as datafile:
    data = json.load(datafile)

values = {}
values['values'] = data
values['name'] = 'fasta'

for key in values['values']:
    if key['marker_type'] in ["A", "G", "T", "C"]:
        key['base'] = key['marker_type']
    elif key['marker_type'] in ["Deletion", "Insertion", "Match", "Pairing"]:
        key['typ'] = key['marker_type']

    if key['marker_type'] == "Insertion":
        key['inserts'] = key['bases']

vdata['width'] = 700

domain = {}
domain = [sys.argv[1], sys.argv[2]]
vdata['scales'][0]['domain'] = domain

vdata['data'][1] = values

json.dump(vdata, sys.stdout, indent=4)