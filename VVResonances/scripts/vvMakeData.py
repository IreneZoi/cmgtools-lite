#!/usr/bin/env python

import ROOT
from array import array
from CMGTools.VVResonances.python.plotting.TreePlotter import TreePlotter
from CMGTools.VVResonances.python.plotting.MergedPlotter import MergedPlotter
from math import log
import os, sys, re, optparse,pickle,shutil,json

print "******* vvMakeData.py start *******"
parser = optparse.OptionParser()
parser.add_option("-s","--samples",dest="samples",default='',help="Type of sample")
parser.add_option("-c","--cut",dest="cut",help="Cut to apply for yield",default='')
parser.add_option("-o","--output",dest="output",help="Output ROOT",default='')
parser.add_option("-v","--vars",dest="vars",help="variables seprataed by comma",default='')
parser.add_option("-b","--bins",dest="bins",help="bins per dimension separated by comma",default='')
parser.add_option("-m","--min",dest="mins",help="minimum separated by comma",default='')
parser.add_option("-M","--max",dest="maxes",help="maximum separated by comma",default='')
parser.add_option("-d","--isData",dest="data",type=int,help="isData",default=1)
parser.add_option("-f","--factors",dest="factors",type=str,help="factor",default="1")
parser.add_option("-n","--name",dest="name",help="name",default="histo")
parser.add_option("--binsMVV",dest="binsMVV",help="use special binning",default="")
parser.add_option("-t","--triggerweight",dest="triggerW",action="store_true",help="Use trigger weights",default=False)

(options,args) = parser.parse_args()
#define output dictionary

samples={}

def getBinning(binsMVV):
    l=[]
    if binsMVV!="":
        s = binsMVV.split(",")
        for w in s:
            l.append(int(w))
    return l

sampleTypes=options.samples.split(',')

dataPlotters=[]

for filename in os.listdir(args[0]):
    for sampleType in sampleTypes:
        if filename.find(sampleType)!=-1:
            fnameParts=filename.split('.')
            fname=fnameParts[0]
            ext=fnameParts[1]
            if ext.find("root") ==-1:
                continue
            dataPlotters.append(TreePlotter(args[0]+'/'+fname+'.root','AnalysisTree'))
            if options.data==0 or options.data==2:
                dataPlotters[-1].setupFromFile(args[0]+'/'+fname+'.pck')
                dataPlotters[-1].addCorrectionFactor('xsec','tree')
                dataPlotters[-1].addCorrectionFactor('genWeight','tree')
                dataPlotters[-1].addCorrectionFactor('puWeight','tree')
                if options.triggerW:
                    print "Using triggerweights! Are you sure?"
                    dataPlotters[-1].addCorrectionFactor('triggerWeight','tree')	
            corrFactors = options.factors.split(',')
	    for c in corrFactors:
	     if len(c.split(':')) < 2: continue
	     if c.split(':')[0] in fname:
	      print "Add correction factor:",fname,c.split(':')[1]
	      dataPlotters[-1].addCorrectionFactor(float(c.split(':')[1]),'flat')



        
    

data=MergedPlotter(dataPlotters)


pvars=options.vars.split(',')
pmins=options.mins.split(',')
pmaxes=options.maxes.split(',')
pbins=options.bins.split(',')

if len(pvars)==1:
    histo=data.drawTH1(pvars[0],options.cut,"1",int(pbins[0]),float(pmins[0]),float(pmaxes[0]))


if len(pvars)==2:
    histo=data.drawTH2(pvars[1]+":"+pvars[0],options.cut,"1",int(pbins[0]),float(pmins[0]),float(pmaxes[0]),int(pbins[1]),float(pmins[1]),float(pmaxes[1]))

if len(pvars)==3:
    binning = getBinning(options.binsMVV)
    print "z ",pvars[2], "y ", pvars[1], "x ", pvars[0]
    print "binsz ",int(pbins[2])," minz ", float(pmins[2])," maxz ", float(pmaxes[2])
    print "binsy ",int(pbins[1])," miny ", float(pmins[1])," maxy ", float(pmaxes[1])
    print "binsx ",int(pbins[0])," minx ", float(pmins[0])," maxx ", float(pmaxes[0])
    if options.binsMVV=="":
        histo=data.drawTH3(pvars[2]+":"+pvars[1]+":"+pvars[0],options.cut,"1",int(pbins[0]),float(pmins[0]),float(pmaxes[0]),int(pbins[1]),float(pmins[1]),float(pmaxes[1]),int(pbins[2]),float(pmins[2]),float(pmaxes[2]))
    else:
        binsx=[]
        for i in range(0,int(pbins[0])+1):
            binsx.append(float(pmins[0])+i*(float(pmaxes[0])-float(pmins[0]))/int(pbins[0]))
        histo=data.drawTH3Binned(pvars[2]+":"+pvars[1]+":"+pvars[0],options.cut,"1",array('f',binsx),array('f',binsx),array('f',binning))

F=ROOT.TFile(options.output,"UPDATE")
F.cd()
histo.Write(options.name)
F.Close()

print "******* vvMakeData.py end *******"

