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

for key in values['values']:
    if key['marker_type'] == "A" or key['marker_type'] == "G" or key['marker_type'] == "T" or key['marker_type'] == "C":
        key['base'] = key['marker_type']
    elif key['marker_type'] == "Deletion" or key['marker_type'] == "Insertion" or key['marker_type'] == "Match" or key['marker_type'] == "Pairing":
        key['typ'] = key['marker_type']

    if key['marker_type'] == "Insertion":
        key['inserts'] = key['bases']

#TODO: fix order that matches go last due to displaying insertions etc. over matches on overlapping reads

vdata['width'] = 700

domain = {}
domain = [sys.argv[1], sys.argv[2]]
vdata['scales'][0]['domain'] = domain

vdata['data'][1] = values
with open('static_vega.json', 'w') as newvspec:
    json.dump(vdata, newvspec, indent=4)

print('created static_vega.json')