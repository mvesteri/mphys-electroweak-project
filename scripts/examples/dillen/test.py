#!/usr/bin/env python

import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

DATADIR="/storage/epp2/phshgg/Public/ew_analyses/run3_tuples/11/"

with uproot.open(f"{DATADIR}/d13600GeV_24c3_Up_Turbo_Stream.root:U1S/DecayTree") as t:

    momenta = t.arrays(library = 'pd')
    print(momenta.columns.tolist())