#!/usr/bin/env python
from ROOT import gROOT
gROOT.SetBatch(True) #prevents ROOT trying to open a GUI
import ROOT
import os
scriptname = os.path.basename(__file__)

def main():

    tree = ROOT.TChain('tree')
    tree.Add('/storage/epp2/phshgg/Public/MPhysProject_2023_2024/data__13TeV__2016__magnet_down__BuToJpsiK.root')

    tree.Show(0)
    
    nbins = 100
    xmin = 5.0e3
    xmax = 5.5e3
    
    hist = ROOT.TH1D('hist_name','',nbins,xmin,xmax)
    hist.SetXTitle("B mass [MeV]")
    hist.SetYTitle("Counts")
    
    hist_plus = hist.Clone('plus')
    hist_minus = hist.Clone('minus')

    for name,op in zip(('plus','minus'),('>','<')):
        tree.Draw(f'Hb_M>>{name}',f'Hb_ID {op} 0','goff')
    
    canv = ROOT.TCanvas()

    for h,colour,opt in zip((hist_plus,hist_minus),(ROOT.kOrange,ROOT.kAzure),('','SAME')):
        h.Scale(1./h.Integral())
        h.SetLineColor(colour)
        h.Draw(opt)

    canv.BuildLegend()
    
    canv.SaveAs(f'transient/{scriptname}.pdf')

if __name__ == '__main__':
    main()
