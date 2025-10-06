#!/usr/bin/env python
from ROOT import gROOT
gROOT.SetBatch(True) #prevents ROOT trying to open a GUI
import ROOT

def main():

    tree = ROOT.TChain('DecayTree')
    tree.Add('/storage/epp2/phshgg/Public/MPhysProject_2023_2024/13TeV__2018__magnet_down_data__Z_candidates.root')

    tree.Show(0)
    
    nbins = 100
    xmin = 15.e3
    xmax = 60.e3

    hist = ROOT.TH1D('hist_name','A nice title for the histogram',nbins,xmin,xmax)

    tree.Draw('mup_PT>>'+hist.GetName())
    
    canv = ROOT.TCanvas()
    
    hist.Draw()
    hist.SetXTitle("Positive muon transverse momentum [MeV]")
    hist.SetYTitle("Counts")
    
    canv.SaveAs('transient/plot.pdf')

if __name__ == '__main__':
    main()
