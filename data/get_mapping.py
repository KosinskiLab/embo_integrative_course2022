import json
import sys

with open(sys.argv[1]) as f:
    j = json.load(f)

sub_chain = {}
for s in j['subunits']:
    sub_chain[s['name']] = s['chainIds'][0]
print(sub_chain)

for d in j['data']:
    if d['type'] == 'sequences':
        for seqn, subunits in d['mapping'].items():
            subunit = subunits[0]
            print(sub_chain[subunit]+', '+seqn)
