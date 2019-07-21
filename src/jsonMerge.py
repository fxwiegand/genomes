import json

with open('client/vlSpec.json', 'r') as vlspec:
    vldata = json.load(vlspec)

with open('data.json', 'r') as datafile:
    data = json.load(datafile)

values = {}
values['values'] = data

vldata['width'] = 700

vldata['data'] = values
with open('static_vega.json', 'w') as newvlspec:
    json.dump(vldata, newvlspec, indent=4)

print('created static_vega.json')