#!/usr/bin/env python
import uproot

with uproot.open('/tmp/13TeV__2018__magnet_down_data__Z_candidates.root:DecayTree') as t:
    momenta = t.arrays(["mup_PT","mup_ETA"],library="np")
    print(momenta)

