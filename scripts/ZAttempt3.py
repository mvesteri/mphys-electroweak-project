#!/usr/bin/env python
import uproot
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
DATAIR="/storage/epp2/phshgg/Public/MPhysProject_2025_2026/tuples/1/"
with uproot.open(f"{DATAIR}/DecayTree__Z__Z__d13600GeV_24c4.root:DecayTree") as t:

    #momentasim = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX" ,"mum_PY" ,"mum_PZ"],library="np")
    datam=t.arrays(["mum_eta","mum_phi","mum_pt","mup_eta" ,"mup_phi" ,"mup_pt","true_boson_mass"],library="np")  #i ahve ztrue and i ahve z reconstructed
    
    tmass=datam["true_boson_mass"]

def conaeq(reconmass):
    mpe=reconmass["mup_eta"]
    mppt=reconmass["mup_pt"] 
    mpphi=reconmass["mup_phi"] 
    mme=reconmass["mum_eta"] 
    mmpt=reconmass["mum_pt"] 
    mmphi=reconmass["mum_phi"]
    mpx=mppt*np.cos(mpphi)
    mmx=mmpt*np.cos(mmphi)
    mpy=mppt*np.sin(mpphi)
    mmy=mmpt*np.sin(mmphi)
    mpz=mppt*np.sinh(mpe)
    mmz=mmpt*np.sinh(mme)

    m=[]
    Ep=[]
    Em=[]
    M_mass=105.658*10**(-3)
    EP=np.sqrt(M_mass**2 +mpx**2 +mpy**2+mpz**2)   
    EN=np.sqrt(M_mass**2 +mmx**2 +mmy**2+mmz**2)
    m=np.sqrt((EP+EN)**2-((mpx+mmx)**2+(mpy+mmy)**2+(mpz+mmz)**2))
    
    return m

def fiteq(x,a,b,m,w,scale):
    gamma=np.sqrt(m**2*(m**2+w**2))
    k=(2*np.sqrt(2)*m*gamma*w)/(np.pi*np.sqrt(m**2+gamma))
    f= ((k/((x**2-m**2)**2+(m**2)*w**2))*scale  + a*np.exp(b*x))  # i can just do this manually
    return f

def plot(tmass,datam):
    remassHist,rebinn,_=plt.hist(conaeq(datam), bins=150, histtype="step",range=(70,110),label="Z-reconstructed",density=True,linewidth=1)
    trmassHist,binn,_=plt.hist(tmass, bins=150, histtype="step",range=(70,110),label="Z-true",density=True,linewidth=1) #for true mass
    centers=0.5*(binn[1:]+binn[:-1])
    w=3
    m=91.1876
    fitParam,_ = curve_fit(fiteq,centers,trmassHist,p0=[5,-0.7,m,w,0.4],bounds=([0.0,-1.0,60,0,0],[100.0,0,120,10,1]),maxfev=10000)
    print(fitParam)
    model = fiteq(centers,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    plt.plot(centers,model)
    plt.title("Z-true fit+reconstructed")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
    plt.savefig("Z-true fit +reconstructed.pdf")
    plt.clf()

    plt.figure()  #now making a grpah with a histogram with woeghtings applied
    ratio=fiteq(tmass,fitParam[0],fitParam[1],85,fitParam[3],fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    massHistweights,binnweights,_weights=plt.hist(conaeq(datam), bins=150, histtype="step",range=(70,100),weights=ratio,label="Z-reconstructed weighted",density=True,linewidth=1) #for reconstructed mass with wieghts
    remassHist,rebinn,_=plt.hist(conaeq(datam), bins=150, histtype="step",range=(70,110),label="Z-reconstructed",density=True,linewidth=1)
    plt.plot(centers,model)
    plt.title("Z-reconstucted weighted target 85")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
    plt.savefig("Z-reconstucted weighted.pdf")
    plt.clf()

    plt.figure()
    remassHist,rebinn,_=plt.hist(conaeq(datam), bins=150, histtype="step",range=(70,110),label="Z-reconstructed",density=True,linewidth=1)
    ratio90=fiteq(tmass,fitParam[0],fitParam[1],90,fitParam[3],fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    massHistweights90,binnweights90,_weights90=plt.hist(conaeq(datam), bins=150, histtype="step",range=(70,110),weights=ratio90,label="Z-reconstructed weighted 90",density=True,linewidth=1) #for reconstructed mass with wieghts
    ratio91=fiteq(tmass,fitParam[0],fitParam[1],91,fitParam[3],fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    massHistweights91,binnweights91,_weights91=plt.hist(conaeq(datam), bins=150, histtype="step",range=(70,110),weights=ratio91,label="Z-reconstructed weighted 91",density=True,linewidth=1)
    ratio92=fiteq(tmass,fitParam[0],fitParam[1],92,fitParam[3],fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    massHistweights92,binnweights92,_weights92=plt.hist(conaeq(datam), bins=150, histtype="step",range=(70,110),weights=ratio92,label="Z-reconstructed weighted 92",density=True,linewidth=1)
    plt.title("Z-reconstucted weights")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
    plt.savefig("weights graph.pdf")
    plt.clf()
plot(tmass,datam)