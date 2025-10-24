#!/usr/bin/env python
import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

DATADIR="/storage/epp2/phshgg/Public/ew_analyses/run3_tuples/11/"

with uproot.open(f"{DATADIR}/d13600GeV_24c3_Up_Full_Stream.root:Z/DecayTree") as t:
    momenta = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX" ,"mum_PY" ,"mum_PZ"],library="np")
    #print(momenta)


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


DATADIRsim="/storage/epp2/phshgg/Public/ew_analyses/run3_tuples/11"

with uproot.open(f"{DATADIRsim}/d13600GeV_block8_Up_Z_Sim10f_42112001_Full.root:Z/DecayTree") as t:
    momentasim = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX" ,"mum_PY" ,"mum_PZ"],library="np")
    #print(momenta)

print("l")
def masssim(momentasim):
    M_mass=105.658
    mpx=momentasim["mup_PX"]
    mpy=momentasim["mup_PY"] 
    mpz=momentasim["mup_PZ"] 
    mmx=momentasim["mum_PX"] 
    mmy=momentasim["mum_PY"] 
    mmz=momentasim["mum_PZ"]
    m=[]
    Ep=[]
    Em=[]
    
    EP=np.sqrt(M_mass**2 +mpx**2 +mpy**2+mpz**2)
        
    EN=np.sqrt(M_mass**2 +mmx**2 +mmy**2+mmz**2)
    v=np.sqrt((EP+EN)**2-((mpx+mmx)**2+(mpy+mmy)**2+(mpz+mmz)**2))
    m.append(v)
    return m
pointsim=masssim(momentasim)

def plot(pointsim,points):

    plt.hist(pointsim, bins=300, histtype="step",range=(35000,150000),label="Z mass sim",density=True, linewidth=1.5,color='blue')
    plt.hist(points, bins=300, histtype="step",range=(35000,150000),label="Z mass data",density=True, linewidth=1.5,color='red')
    plt.title("Z invariant mass")
    plt.xlabel("Mass_Mev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
     #binn is array of bin edges
    #centers=0.5*(binn[1:]+binn[:-1])
    #initial_guess = [max(counts), -0.01,10, np.mean(points), 10] # np.std(points)
    #popt, pcov = curve_fit(fiteq, centers, counts, p0=initial_guess,bounds=([-np.inf,-np.inf,-np.inf,-np.inf,0],[np.inf,0,np.inf,np.inf,100]))
    #plt.plot(centers, fiteq(centers, *popt), 'r-', label='Gaussian+ exponential fit')
    plt.savefig("Zgraph_ds.pdf")
    plt.clf()
    print("gr")

plot(pointsim,points)