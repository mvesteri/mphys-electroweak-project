#!/usr/bin/env python
import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
DATADIR="/storage/epp2/phshgg/Public/ew_analyses/run3_tuples/11/"

with uproot.open(f"{DATADIR}/d13600GeV_24c3_Up_Full_Stream.root:Z/DecayTree") as t:
    momenta = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX" ,"mum_PY" ,"mum_PZ"],library="np")
    #print(momenta)

print("l")
def mass(momenta):
    M_mass=105.658
    mpx=momenta["mup_PX"]
    mpy=momenta["mup_PY"] 
    mpz=momenta["mup_PZ"] 
    mmx=momenta["mum_PX"] 
    mmy=momenta["mum_PY"] 
    mmz=momenta["mum_PZ"]
    m=[]
    Ep=[]
    Em=[]
    
    EP=np.sqrt(M_mass**2 +mpx**2 +mpy**2+mpz**2)
        
    EN=np.sqrt(M_mass**2 +mmx**2 +mmy**2+mmz**2)
    v=np.sqrt((EP+EN)**2-((mpx+mmx)**2+(mpy+mmy)**2+(mpz+mmz)**2))
    m.append(v)
    return m
points=mass(momenta)

def plot(points):

    plt.hist(points, bins=300, histtype="step",range=(35000,150000),label="Z mass",density=True)
    plt.title("Z invariant mass")
    plt.xlabel("Mass_Mev")
    plt.ylabel("Frequency Density")
    
     #binn is array of bin edges
    #centers=0.5*(binn[1:]+binn[:-1])
    #initial_guess = [max(counts), -0.01,10, np.mean(points), 10] # np.std(points)
    #popt, pcov = curve_fit(fiteq, centers, counts, p0=initial_guess,bounds=([-np.inf,-np.inf,-np.inf,-np.inf,0],[np.inf,0,np.inf,np.inf,100]))
    #plt.plot(centers, fiteq(centers, *popt), 'r-', label='Gaussian+ exponential fit')
    plt.savefig("Zgraph.pdf")
    plt.clf()
    print("gr")

plot(points)