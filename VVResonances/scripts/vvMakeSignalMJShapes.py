#!/usr/bin/env python

import ROOT
from array import array
from CMGTools.VVResonances.python.plotting.TreePlotter import TreePlotter
from CMGTools.VVResonances.python.plotting.MergedPlotter import MergedPlotter
from CMGTools.VVResonances.python.plotting.StackPlotter import StackPlotter
from CMGTools.VVResonances.python.statistics.Fitter import Fitter
from math import log
import os, sys, re, optparse,pickle,shutil,json
ROOT.gROOT.SetBatch(True)

print " **************    vvMakeSignalMJShapes   start   **************"

def returnString(func):
    st='0'
    for i in range(0,func.GetNpar()):
        st=st+"+("+str(func.GetParameter(i))+")"+("*MH"*i)
    return st    


parser = optparse.OptionParser()
parser.add_option("-s","--sample",dest="sample",default='',help="Type of sample")
parser.add_option("-c","--cut",dest="cut",help="Cut to apply for shape",default='')
parser.add_option("-o","--output",dest="output",help="Output JSON",default='')
parser.add_option("-V","--MVV",dest="mvv",help="mVV variable",default='')
parser.add_option("-m","--min",dest="mini",type=float,help="min MJJ",default=40)
parser.add_option("-M","--max",dest="maxi",type=float,help="max MJJ",default=160)
parser.add_option("-e","--exp",dest="doExp",type=int,help="useExponential",default=0)
parser.add_option("-f","--fix",dest="fixPars",help="Fixed parameters",default="1")
parser.add_option("-r","--minMX",dest="minMX",type=float, help="smallest Mx to fit ",default=1000.0)
parser.add_option("-R","--maxMX",dest="maxMX",type=float, help="largest Mx to fit " ,default=7000.0)
parser.add_option("-t","--triggerweight",dest="triggerW",action="store_true",help="Use trigger weights",default=False)

(options,args) = parser.parse_args()
#define output dictionary

isVH = False
samples={}

for filename in os.listdir(args[0]):
    print filename
    if not (filename.find(options.sample)!=-1):
        continue

    if filename.find('hbb'):
     isVH=True
     print "INFO: fitting VH sample with double peak"
     
    fnameParts=filename.split('.')
    fname=fnameParts[0]
    print fname
    ext=fnameParts[1]
    print ext
    if ext.find("root") ==-1:
        continue
        
    mass = float(fname.split('_')[-1])
    print mass
    if mass < options.minMX or mass > options.maxMX: continue	
    samples[mass] = fname
    
    print 'found',filename,'mass',str(mass) 
print options.mvv
leg = options.mvv.split('_')[1]
print leg
graphs={'mean':ROOT.TGraphErrors(),'sigma':ROOT.TGraphErrors(),'alpha':ROOT.TGraphErrors(),'n':ROOT.TGraphErrors(),'f':ROOT.TGraphErrors(),'slope':ROOT.TGraphErrors(),'alpha2':ROOT.TGraphErrors(),'n2':ROOT.TGraphErrors(),
        'meanH':ROOT.TGraphErrors(),'sigmaH':ROOT.TGraphErrors(),'alphaH':ROOT.TGraphErrors(),'nH':ROOT.TGraphErrors(),'fH':ROOT.TGraphErrors(),'slopeH':ROOT.TGraphErrors(),'alpha2H':ROOT.TGraphErrors(),'n2H':ROOT.TGraphErrors() }

#Now we have the samples: Sort the masses and run the fits
N=0
print "*** start cycle on masses in vvMakeSignalMJShapes***"
for mass in sorted(samples.keys()):
    print "in vvMakeSignalMJShapes mass "+str(mass)
    print 'fitting',str(mass) 
    print args[0]
    print samples[mass]
    plotter=TreePlotter(args[0]+'/'+samples[mass]+'.root','AnalysisTree')
    plotter.addCorrectionFactor('genWeight','tree')
    plotter.addCorrectionFactor('puWeight','tree')
    if options.triggerW:
     plotter.addCorrectionFactor('jj_triggerWeight','tree')
     print "Using triggerweight"
       
        
    fitter=Fitter(['x'])
    print " *** called Fitter from statistics ***"
    if isVH:
        print "*** fitter.jetDoublePeakVH ***"
        fitter.jetDoublePeakVH('model','x')
    else:
        print "*** fitter.jetResonanceNOEXP ***"
        fitter.jetResonanceNOEXP('model','x')
    
    if options.fixPars!="1":
        fixedPars =options.fixPars.split(',')
        for par in fixedPars:
            parVal = par.split(':')
	    if len(parVal) > 1:
             fitter.w.var(parVal[0]).setVal(float(parVal[1]))
             fitter.w.var(parVal[0]).setConstant(1)

    print "*** preparing hist ***"
#    fitter.w.var("MH").setVal(mass)
    histo = plotter.drawTH1(options.mvv,options.cut,"1",int((options.maxi-options.mini)/4),options.mini,options.maxi)

    fitter.importBinnedData(histo,['x'],'data')
    print "*****  fitter.fit('model','data',[ROOT.RooFit.SumW2Error(0)]) *****"
    fitter.fit('model','data',[ROOT.RooFit.SumW2Error(0)])
    print "*****  fitter.fit('model','data',[ROOT.RooFit.SumW2Error(0),ROOT.RooFit.Minos(1)]) *****"
    fitter.fit('model','data',[ROOT.RooFit.SumW2Error(0),ROOT.RooFit.Minos(1)])
    print "*****  fitter.projection *****"
    fitter.projection("model","data","x","debugJ"+leg+"_"+options.output+"_"+str(mass)+".png")

    print "*****  preparing graphs and ouptput *****"
    for var,graph in graphs.iteritems():
        value,error=fitter.fetch(var)
        graph.SetPoint(N,mass,value)
        graph.SetPointError(N,0.0,error)
                
    N=N+1
    fitter.delete()
    print "**** in vvMakeSignalMJShapes done with mass "+str(mass)

print "*** end cycle on masses ***"        
F=ROOT.TFile(options.output,"RECREATE")
F.cd()
for name,graph in graphs.iteritems():
    graph.Write(name)
F.Close()
            
print " **************    vvMakeSignalMJShapes   end   **************"
