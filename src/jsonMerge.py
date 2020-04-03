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

def compare(item1, item2):
    if item1 == "Match":
        return -1
    elif item2 == "Match":
        return 1
    elif item1 in ["A", "G", "T", "C"]:
        return -1
    elif item2 in ["A", "G", "T", "C"]:
        return 1
    else:
        return 0


values['values'] = sorted(values['values'], cmp= lambda i1, i2: compare(i1['marker_type'], i2['marker_type']))

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
with open('static_vega.json', 'w') as newvspec:
    json.dump(vdata, newvspec, indent=4)

print('created static_vega.json')