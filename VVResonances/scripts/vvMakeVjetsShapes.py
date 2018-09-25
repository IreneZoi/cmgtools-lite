#!/usr/bin/env python

import ROOT
from array import array
from CMGTools.VVResonances.plotting.TreePlotter import TreePlotter
from CMGTools.VVResonances.plotting.MergedPlotter import MergedPlotter
from CMGTools.VVResonances.plotting.StackPlotter import StackPlotter
from CMGTools.VVResonances.statistics.Fitter import Fitter
from math import log
import os, sys, re, optparse,pickle,shutil,json
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)



parser = optparse.OptionParser()
parser.add_option("-s","--sample",dest="sample",default='',help="Type of sample")
parser.add_option("-c","--cut",dest="cut",help="Cut to apply for shape",default='')
parser.add_option("-o","--output",dest="output",help="Output JSON",default='')
parser.add_option("-m","--min",dest="mini",type=float,help="min MJJ",default=40)
parser.add_option("-M","--max",dest="maxi",type=float,help="max MJJ",default=160)
parser.add_option("--store",dest="store",type=str,help="store fitted parameters in this file",default="")
parser.add_option("--corrFactorW",dest="corrFactorW",type=float,help="add correction factor xsec",default=1.)
parser.add_option("--corrFactorZ",dest="corrFactorZ",type=float,help="add correction factor xsec",default=41.34/581.8)
parser.add_option("-f","--fix",dest="fixPars",help="Fixed parameters",default="1")
parser.add_option("--minMVV","--minMVV",dest="minMVV",type=float,help="mVV variable",default=1)
parser.add_option("--maxMVV","--maxMVV",dest="maxMVV",type=float, help="mVV variable",default=1)
parser.add_option("--binsMVV",dest="binsMVV",help="use special binning",default="")
parser.add_option("-t","--triggerweight",dest="triggerW",action="store_true",help="Use trigger weights",default=False)

(options,args) = parser.parse_args()
samples={}

def getBinning(binsMVV,minx,maxx,bins):
    l=[]
    if binsMVV=="":
        for i in range(0,bins+1):
            l.append(minx + i* (maxx - minx)/bins)
    else:
        s = binsMVV.split(",")
        for w in s:
            l.append(int(w))
    return l
    
label = options.output.split(".root")[0]
t  = label.split("_")
el=""
for words in t:
    if words.find("HP")!=-1 or words.find("LP")!=-1:
        continue
    el+=words+"_"
label = el

samplenames = options.sample.split(",")
for filename in os.listdir(args[0]):
    for samplename in samplenames:
        if not (filename.find(samplename)!=-1):
            continue


        fnameParts=filename.split('.')
        fname=fnameParts[0]
        ext=fnameParts[1]
        if ext.find("root") ==-1:
            continue
    
        name = fname.split('_')[0]
    
        samples[name] = fname

        print 'found',filename

sigmas=[]
params={}
legs=["l1","l2"]

plotters=[]
for name in samples.keys():
    plotters.append(TreePlotter(args[0]+'/'+samples[name]+'.root','tree'))
    plotters[-1].setupFromFile(args[0]+'/'+samples[name]+'.pck')
    plotters[-1].addCorrectionFactor('xsec','tree')
    plotters[-1].addCorrectionFactor('genWeight','tree')
    plotters[-1].addCorrectionFactor('puWeight','tree')
    if options.triggerW: plotters[-1].addCorrectionFactor('triggerWeight','tree')	

    corrFactor = options.corrFactorW
    if samples[name].find('Z') != -1:
        corrFactor = options.corrFactorZ
    if samples[name].find('W') != -1:
        corrFactor = options.corrFactorW
    plotters[-1].addCorrectionFactor(corrFactor,'flat')
        
plotter=MergedPlotter(plotters)
print 'Fitting Mjet:' 

for leg in legs:
 tmp=[]
 tmp_nonres=[]
 fitter=Fitter(['x'])
 fitter.jetResonanceVjets('model','x')

 if options.fixPars!="1":
     fixedPars =options.fixPars.split(',')
     if len(fixedPars) > 1:
      print "   - Fix parameters: ", fixedPars
      for par in fixedPars:
       if par=="c_0" or par =="c_1" or par=="c_2": continue
       parVal = par.split(':')
       fitter.w.var(parVal[0]).setVal(float(parVal[1]))
       fitter.w.var(parVal[0]).setConstant(1)

 
 histo = plotter.drawTH1("jj_"+leg+"_softDrop_mass",options.cut+"*(jj_"+leg+"_mergedVTruth==1)*(jj_"+leg+"_softDrop_mass>55&&jj_"+leg+"_softDrop_mass<215)","1",80,55,215)
 
 ### fit distribution with a gaussian function and fix parameters of doubleCB ###
 gauss  = ROOT.TF1("gauss" ,"gaus",74,94)  
 histo.Fit(gauss,"R")
 mean = gauss.GetParameter(1)
 sigma = gauss.GetParameter(2)

 print "set paramters of double CB constant aground the ones from gaussian fit"
 fitter.w.var("mean").setVal(mean)
 fitter.w.var("mean").setConstant(1)
 fitter.w.var("sigma").setVal(sigma)
 fitter.w.var("sigma").setConstant(1)
 #######################################################
 
 fitter.importBinnedData(histo,['x'],'data')
 fitter.fit('model','data',[ROOT.RooFit.SumW2Error(1),ROOT.RooFit.Save(1),ROOT.RooFit.Range(55,120)])
 fitter.projection("model","data","x","debugJ"+leg+"_"+options.output+"_Res.png")
 params[label+"_Res_"+leg]={"mean": {"val": fitter.w.var("mean").getVal(), "err": fitter.w.var("mean").getError()}, "sigma": {"val": fitter.w.var("sigma").getVal(), "err": fitter.w.var("sigma").getError()}, "alpha":{ "val": fitter.w.var("alpha").getVal(), "err": fitter.w.var("alpha").getError()},"alpha2":{"val": fitter.w.var("alpha2").getVal(),"err": fitter.w.var("alpha2").getError()},"n":{ "val": fitter.w.var("n").getVal(), "err": fitter.w.var("n").getError()},"n2": {"val": fitter.w.var("n2").getVal(), "err": fitter.w.var("n2").getError()}}

   
        
if options.store!="":
    f=open(options.store,"w")
    for par in params:
        f.write(str(par)+ " = " +str(params[par])+"\n")
