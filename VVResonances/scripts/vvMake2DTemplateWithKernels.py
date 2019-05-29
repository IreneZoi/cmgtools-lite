#!/usr/bin/env python
import ROOT
from array import array
from CMGTools.VVResonances.python.statistics.Fitter import Fitter
from math import log,exp,sqrt
import os, sys, re, optparse,pickle,shutil,json
import json
import copy
from CMGTools.VVResonances.python.plotting.tdrstyle import *
setTDRStyle()
from CMGTools.VVResonances.python.plotting.TreePlotter import TreePlotter
from CMGTools.VVResonances.python.plotting.MergedPlotter import MergedPlotter
ROOT.gSystem.Load("libCMGToolsVVResonances")
parser = optparse.OptionParser()
parser.add_option("-o","--output",dest="output",help="Output",default='')
parser.add_option("-r","--res",dest="res",help="res",default='')
parser.add_option("-s","--samples",dest="samples",default='',help="Type of sample")
parser.add_option("-c","--cut",dest="cut",help="Cut to apply for yield in gen sample",default='')
parser.add_option("-v","--vars",dest="vars",help="variable for x",default='')
parser.add_option("-b","--binsx",dest="binsx",type=int,help="bins",default=1)
parser.add_option("-B","--binsy",dest="binsy",type=int,help="conditional bins split by comma",default=1)
parser.add_option("-x","--minx",dest="minx",type=float,help="bins",default=0)
parser.add_option("-X","--maxx",dest="maxx",type=float,help="conditional bins split by comma",default=1)
parser.add_option("-y","--miny",dest="miny",type=float,help="bins",default=0)
parser.add_option("-Y","--maxy",dest="maxy",type=float,help="conditional bins split by comma",default=1)
parser.add_option("-w","--weights",dest="weights",help="additional weights",default='')
parser.add_option("-u","--usegenmass",dest="usegenmass",action="store_true",help="use gen mass for det resolution",default=False)
parser.add_option("-e","--firstEv",dest="firstEv",type=int,help="first event",default=0)
parser.add_option("-E","--lastEv",dest="lastEv",type=int,help="last event",default=-1)
parser.add_option("--binsMVV",dest="binsMVV",help="use special binning",default="")
parser.add_option("-t","--triggerweight",dest="triggerW",action="store_true",help="Use trigger weights",default=False)

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



def unequalScale(histo,name,alpha,power=1):
    newHistoU =copy.deepcopy(histo) 
    newHistoU.SetName(name+"Up")
    newHistoD =copy.deepcopy(histo) 
    newHistoD.SetName(name+"Down")
    maxFactor = max(pow(histo.GetXaxis().GetXmax(),power),pow(histo.GetXaxis().GetXmin(),power))
    for i in range(1,histo.GetNbinsX()+1):
        x= histo.GetXaxis().GetBinCenter(i)
        for j in range(1,histo.GetNbinsY()+1):
            nominal=histo.GetBinContent(i,j)
            factor = 1+alpha*pow(x,power) 
            newHistoU.SetBinContent(i,j,nominal*factor)
            newHistoD.SetBinContent(i,j,nominal/factor)
    if newHistoU.Integral()>0.0:        
        newHistoU.Scale(1.0/newHistoU.Integral())        
    if newHistoD.Integral()>0.0:        
        newHistoD.Scale(1.0/newHistoD.Integral())        
    return newHistoU,newHistoD 
    
def mirror(histo,histoNominal,name):
    newHisto =copy.deepcopy(histoNominal) 
    newHisto.SetName(name)
    intNominal=histoNominal.Integral()
    intUp = histo.Integral()
    for i in range(1,histo.GetNbinsX()+1):
        for j in range(1,histo.GetNbinsY()+1):
            up=histo.GetBinContent(i,j)/intUp
            nominal=histoNominal.GetBinContent(i,j)/intNominal
            newHisto.SetBinContent(i,j,histoNominal.GetBinContent(i,j)*nominal/up)
    return newHisto       
	
def expandHisto(histo,options):
    histogram=ROOT.TH2F(histo.GetName(),"histo",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))
    for i in range(1,histo.GetNbinsX()+1):
        proje = histo.ProjectionY("q",i,i)
        graph=ROOT.TGraph(proje)
        for j in range(1,histogram.GetNbinsY()+1):
            x=histogram.GetYaxis().GetBinCenter(j)
            bin=histogram.GetBin(i,j)
            histogram.SetBinContent(bin,graph.Eval(x,0,"S"))
    return histogram
      
def conditional(hist):
    for i in range(1,hist.GetNbinsY()+1):
        proj=hist.ProjectionX("q",i,i)
        integral=proj.Integral()
        if integral==0.0:
            print 'SLICE WITH NO EVENTS!!!!!!!!',hist.GetName()
            continue
        for j in range(1,hist.GetNbinsX()+1):
            hist.SetBinContent(j,i,hist.GetBinContent(j,i)/integral)


(options,args) = parser.parse_args()

weights_ = options.weights.split(',')

random=ROOT.TRandom3(101082)

sampleTypes=options.samples.split(',')
print "Creating datasets for samples: " ,sampleTypes

dataPlotters=[]
dataPlottersNW=[]

for filename in os.listdir(args[0]):
 for sampleType in sampleTypes:
  if filename.find(sampleType)!=-1:
   fnameParts=filename.split('.')
   fname=fnameParts[0]
   ext=fnameParts[1]
   if ext.find("root") ==-1:
       continue
   dataPlotters.append(TreePlotter(args[0]+'/'+fname+'.root','AnalysisTree'))
   dataPlotters[-1].setupFromFile(args[0]+'/'+fname+'.pck')
   dataPlotters[-1].addCorrectionFactor('xsec','tree')
   dataPlotters[-1].addCorrectionFactor('genWeight','tree')
   dataPlotters[-1].addCorrectionFactor('puWeight','tree')
   if options.triggerW: dataPlotters[-1].addCorrectionFactor('triggerWeight','tree')
   for w in weights_:
    if w != '': dataPlotters[-1].addCorrectionFactor(w,'branch')
   dataPlotters[-1].filename=fname
   dataPlottersNW.append(TreePlotter(args[0]+'/'+fname+'.root','AnalysisTree'))
   dataPlottersNW[-1].addCorrectionFactor('puWeight','tree')
   dataPlottersNW[-1].addCorrectionFactor('genWeight','tree')
   if options.triggerW: dataPlottersNW[-1].addCorrectionFactor('triggerWeight','tree')
   for w in weights_: 
    if w != '': dataPlottersNW[-1].addCorrectionFactor(w,'branch')
   dataPlottersNW[-1].filename=fname

data=MergedPlotter(dataPlotters)

fcorr=ROOT.TFile(options.res)
scale_x=fcorr.Get("scaleyHisto")
scale_y=fcorr.Get("scalexHisto")
res_x=fcorr.Get("resyHisto")
res_y=fcorr.Get("resxHisto")

variables=options.vars.split(',')
leg = options.vars.split(',')[0].split('_')[1]



binsx=[]
for i in range(0,options.binsx+1):
    binsx.append(options.minx+i*(options.maxx-options.minx)/options.binsx)

binsy = getBinning(options.binsMVV,options.miny,options.maxy,options.binsy)
print binsy

scaleUp = ROOT.TH1F(scale_x)
scaleUp.SetName("scaleUp")
scaleDown = ROOT.TH1F(scale_x)
scaleDown.SetName("scaleDown")
for i in range(1,scale_x.GetNbinsX()+1):
    scaleUp.SetBinContent(i,scale_x.GetBinContent(i)+0.09)
    scaleDown.SetBinContent(i,scale_x.GetBinContent(i)-0.09)
    
histogram = ROOT.TH2F("histo_nominal","histo_nominal",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))
histogram_scale_up = ROOT.TH2F("histo_nominal_ScaleUp","histo_nominal_ScaleUp",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))
histogram_scale_down = ROOT.TH2F("histo_nominal_ScaleDown","histo_nominal_ScaleDown",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))

histogram_altshapeUp = ROOT.TH2F("histo_altshapeUp","histo_altshapeUp",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))
histogram_altshape_scale_up = ROOT.TH2F("histo_altshape_ScaleUp","histo_altshape_ScaleUp",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))
histogram_altshape_scale_down = ROOT.TH2F("histo_altshape_ScaleDown","histo_altshape_ScaleDown",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))

histogram_altshape2 = ROOT.TH2F("histo_altshape2","histo_altshape2",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))

#systematics
histograms=[
    histogram,
    #histogram_scale_up,
    #histogram_scale_down,
    #histogram_altshapeUp,
    #histogram_altshape_scale_up,
    #histogram_altshape_scale_down,    
    #histogram_altshape2
]

#ok lets populate!

maxEvents = -1
varsDataSet = 'jj_%s_gen_pt,%s,%s'%(leg,variables[1],variables[0])

for plotter,plotterNW in zip(dataPlotters,dataPlottersNW):
 
 #Nominal histogram Pythia8
 if plotter.filename.find(sampleTypes[0].replace('.root','')) != -1:
  print "Preparing nominal histogram for sampletype " ,sampleTypes[0]
  print "filename: ", plotter.filename, " preparing central values histo"
 
  #y:x
  histI2D=plotter.drawTH2Binned("jj_LV_mass:jj_%s_softDrop_mass"%(leg),options.cut,"1",array('f',binsx),array('f',binsy),"Softdrop mass","M_{JJ} mass","GeV","GeV","COLZ" )

  print " - Creating dataset - "
  dataset=plotterNW.makeDataSet(varsDataSet,options.cut,options.firstEv,options.lastEv)

  print " - Creating 2D gaussian template - "
  histTMP=ROOT.TH2F("histoTMP","histo",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))
  if not(options.usegenmass): 
   datamaker=ROOT.cmg.GaussianSumTemplateMaker(dataset,variables[0],variables[1],'jj_%s_gen_pt'%(leg),scale_x,scale_y,res_x,res_y,histTMP)
  else: datamaker=ROOT.cmg.GaussianSumTemplateMaker(dataset,variables[0],variables[1],'jj_%s_gen_softDrop_mass'%(leg),scale_x,scale_y,res_x,res_y,histTMP)

  if histTMP.Integral()>0:
   histTMP.Scale(histI2D.Integral()/histTMP.Integral())
   histogram.Add(histTMP)
   # mjet_mvv_nominal.Add(histI2D)
   
  histTMP.Delete()    
  histI2D.Delete()	
    
 if len(sampleTypes)<2: continue 
 elif plotter.filename.find(sampleTypes[1].replace('.root','')) != -1: #alternative shape Herwig
  print "Preparing alternative shapes for sampletype " ,sampleTypes[1]
  print "filename: ", plotter.filename, " preparing alternate shape histo"

  histI2D=plotter.drawTH2Binned("jj_LV_mass:jj_%s_softDrop_mass"%(leg),options.cut,"1",array('f',binsx),array('f',binsy),"M_{qV} mass","GeV","Softdrop mass","GeV","COLZ" )

  print " - Creating dataset - "
  dataset=plotterNW.makeDataSet(varsDataSet,options.cut,options.firstEv,options.lastEv)

  print " - Creating 2D gaussian template - "
  histTMP=ROOT.TH2F("histoTMP","histo",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))
  if not(options.usegenmass): 
   datamaker=ROOT.cmg.GaussianSumTemplateMaker(dataset,variables[0],variables[1],'jj_%s_gen_pt'%(leg),scale_x,scale_y,res_x,res_y,histTMP)
  else: datamaker=ROOT.cmg.GaussianSumTemplateMaker(dataset,variables[0],variables[1],'jj_%s_gen_softDrop_mass'%(leg),scale_x,scale_y,res_x,res_y,histTMP)

  if histTMP.Integral()>0:
    histTMP.Scale(histI2D.Integral()/histTMP.Integral())
    histogram_altshapeUp.Add(histTMP)
    # mjet_mvv_altshapeUp.Add(histI2D)
   
  histTMP.Delete()

  print " - Creating 2D gaussian template scale up - "
  histTMP=ROOT.TH2F("histoTMP","histo",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))
  if not(options.usegenmass): 
   datamaker=ROOT.cmg.GaussianSumTemplateMaker(dataset,variables[0],variables[1],'jj_%s_gen_pt'%(leg),scaleUp,scale_y,res_x,res_y,histTMP)
  else: datamaker=ROOT.cmg.GaussianSumTemplateMaker(dataset,variables[0],variables[1],'jj_%s_gen_softDrop_mass'%(leg),scaleUp,scale_y,res_x,res_y,histTMP)

  if histTMP.Integral()>0:
    histTMP.Scale(histI2D.Integral()/histTMP.Integral())
    histogram_altshape_scale_up.Add(histTMP)
    
  histTMP.Delete()

  print " - Creating 2D gaussian template scale down - "
  histTMP=ROOT.TH2F("histoTMP","histo",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))
  if not(options.usegenmass): 
   datamaker=ROOT.cmg.GaussianSumTemplateMaker(dataset,variables[0],variables[1],'jj_%s_gen_pt'%(leg),scaleDown,scale_y,res_x,res_y,histTMP)
  else: datamaker=ROOT.cmg.GaussianSumTemplateMaker(dataset,variables[0],variables[1],'jj_%s_gen_softDrop_mass'%(leg),scaleDown,scale_y,res_x,res_y,histTMP)

  if histTMP.Integral()>0:
    histTMP.Scale(histI2D.Integral()/histTMP.Integral())
    histogram_altshape_scale_down.Add(histTMP)
    
  histI2D.Delete()
  histTMP.Delete()
    
 if len(sampleTypes)<3: continue 
 elif plotter.filename.find(sampleTypes[2].replace('.root','')) != -1: #alternative shape Pythia8+Madgraph (not used for syst but only for cross checks)
  print "Preparing alternative shapes for sampletype " ,sampleTypes[2]
  print "filename: ", plotter.filename, " preparing alternate shape histo"

  histI2D=plotter.drawTH2("jj_LV_mass:jj_%s_softDrop_mass"%(leg),options.cut,"1",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy),"M_{qV} mass","GeV","Softdrop mass","GeV","COLZ" )
  histTMP=ROOT.TH2F("histoTMP","histo",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))

  print " - Creating dataset - "
  dataset=plotterNW.makeDataSet(varsDataSet,options.cut,options.firstEv,options.lastEv)

  print " - Creating 2D gaussian template - "
  if not(options.usegenmass): 
   datamaker=ROOT.cmg.GaussianSumTemplateMaker(dataset,variables[0],variables[1],'jj_%s_gen_pt'%(leg),scale_x,scale_y,res_x,res_y,histTMP)
  else: datamaker=ROOT.cmg.GaussianSumTemplateMaker(dataset,variables[0],variables[1],'jj_%s_gen_softDrop_mass'%(leg),scale_x,scale_y,res_x,res_y,histTMP)

  if histTMP.Integral()>0:
    histTMP.Scale(histI2D.Integral()/histTMP.Integral())
    histogram_altshape2.Add(histTMP) 
    # mjet_mvv_altshape2.Add(histI2D)
    
  histI2D.Delete()
  histTMP.Delete()

f=ROOT.TFile(options.output,"RECREATE")
print "Finished producing histograms! Saving to" ,options.output
finalHistograms={}
f.cd()
for hist in histograms:
 print "Working on histogram " ,hist.GetName()
 hist.Write(hist.GetName()+"_coarse")
 print "Creating conditional histogram for ",hist.GetName()
 conditional(hist)
 print "Expanding for " ,hist.GetName()
 expanded=expandHisto(hist,options)
 conditional(expanded)
 expanded.Write()
 finalHistograms[hist.GetName()]=expanded

# ##Mirror Herwig shape
#histogram_altshapeDown=mirror(finalHistograms['histo_altshapeUp'],finalHistograms['histo_nominal'],"histo_altshapeDown")
#conditional(histogram_altshapeDown)
#histogram_altshapeDown.Write()

print options.maxx,options.minx
alpha=1.5/(options.maxx)
histogram_pt_down,histogram_pt_up=unequalScale(finalHistograms['histo_nominal'],"histo_nominal_PT",alpha)
conditional(histogram_pt_down)
histogram_pt_down.Write()
conditional(histogram_pt_up)
histogram_pt_up.Write()

alpha=1.5*(options.minx)
h1,h2=unequalScale(finalHistograms['histo_nominal'],"histo_nominal_OPT",alpha,-1)
conditional(h1)
h1.Write()
conditional(h2)
h2.Write()
        		
f.Close()




