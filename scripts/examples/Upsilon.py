#!/usr/bin/env python

import uproot
import numpy as np
from matplotlib import pyplot as plt

def main():
    DATADIR="/storage/epp2/phshgg/Public/ew_analyses/run3_tuples/11/"
    MUON_MASS = 105.7
    with uproot.open(f"{DATADIR}/d13600GeV_24c3_Up_Turbo_Stream.root:U1S/DecayTree") as t:
        data = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX","mum_PY","mum_PZ"],library="np")
    
    mup_P,mum_P = np.array([data["mup_PX"],data["mup_PY"],data["mup_PZ"]]),np.array([data["mum_PX"],data["mum_PY"],data["mum_PZ"]])
    mup_E,mum_E = np.sqrt(mup_P[0]**2+mup_P[1]**2+mup_P[2]**2+MUON_MASS**2),np.sqrt(mum_P[0]**2+mum_P[1]**2+mum_P[2]**2+MUON_MASS**2)
    tot_E = mup_E + mum_E
    tot_PX = mup_P[0] + mum_P[0]
    tot_PY = mup_P[1] + mum_P[1]
    tot_PZ = mup_P[2] + mum_P[2]
    tot_P = np.sqrt(tot_PX**2+tot_PY**2+tot_PZ**2)
    mass = np.sqrt(tot_E**2 - tot_P**2)
    print(mass)
    dataHist,bins,_ = plt.hist(mass,bins=100,range=(9100,9600),histtype='step',label="Upsilon mass",density=True)
    plt.legend()
    plt.xlabel("Mass / MeV")
    plt.ylabel("Frequency")
    plt.title("Upsilon invariant mass")
    plt.savefig("Upsilon_mass.pdf")
    plt.clf()

if __name__ == '__main__':
    main()