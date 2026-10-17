#!/usr/bin/env python

import uproot

DATADIR="/storage/epp2/phshgg/Public/ew_analyses/run3_tuples/11/"

with uproot.open(f"{DATADIR}/d13600GeV_24c3_Up_Turbo_Stream.root:U1S/DecayTree") as t:
    momenta = t.arrays(["mup_PX","mup_PY","mup_PZ"],library="np")
    print(momenta)

