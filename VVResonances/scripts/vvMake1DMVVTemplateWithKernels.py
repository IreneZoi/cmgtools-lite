#!/usr/bin/env python
import ROOT
from array import array
from CMGTools.VVResonances.python.statistics.Fitter import Fitter
from math import log,exp,sqrt
import os, sys, re, optparse,pickle,shutil,json
import copy
import json
from CMGTools.VVResonances.python.plotting.tdrstyle import *
setTDRStyle()
from CMGTools.VVResonances.python.plotting.TreePlotter import TreePlotter
from CMGTools.VVResonances.python.plotting.MergedPlotter import MergedPlotter
ROOT.gSystem.Load("libCMGToolsVVResonances")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

parser = optparse.OptionParser()
parser.add_option("-o","--output",dest="output",help="Output",default='')
parser.add_option("-r","--res",dest="res",help="res",default='')
parser.add_option("-H","--resHisto",dest="resHisto",help="res",default='')
parser.add_option("-s","--samples",dest="samples",default='',help="Type of sample")
parser.add_option("-c","--cut",dest="cut",help="Cut to apply for yield in gen sample",default='')
parser.add_option("-v","--var",dest="var",help="variable for x",default='')
parser.add_option("-b","--bins",dest="binsx",type=int,help="bins",default=1)
parser.add_option("-x","--minx",dest="minx",type=float,help="bins",default=0)
parser.add_option("-X","--maxx",dest="maxx",type=float,help="conditional bins split by comma",default=1)
parser.add_option("-w","--weights",dest="weights",help="additional weights",default='')
parser.add_option("-u","--usegenmass",dest="usegenmass",action="store_true",help="use gen mass for det resolution",default=False)
parser.add_option("-e","--firstEv",dest="firstEv",type=int,help="first event",default=0)
parser.add_option("-E","--lastEv",dest="lastEv",type=int,help="last event",default=-1)
parser.add_option("--binsMVV",dest="binsMVV",help="use special binning",default="")
parser.add_option("-t","--triggerweight",dest="triggerW",action="store_true",help="Use trigger weights",default=False)
parser.add_option("--corrFactorW",dest="corrFactorW",type=float,help="add correction factor xsec",default=1.)
parser.add_option("--corrFactorZ",dest="corrFactorZ",type=float,help="add correction factor xsec",default=1.)


(options,args) = parser.parse_args()


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
    for i in range(1,histo.GetNbinsX()+1):
        x= histo.GetXaxis().GetBinCenter(i)
        nominal=histo.GetBinContent(i)
        factor = 1+alpha*pow(x,power) 
        newHistoU.SetBinContent(i,nominal*factor)
        newHistoD.SetBinContent(i,nominal/factor)
    return newHistoU,newHistoD 
    
def mirror(histo,histoNominal,name):
    newHisto =copy.deepcopy(histoNominal) 
    newHisto.SetName(name)
    intNominal=histoNominal.Integral()
    intUp = histo.Integral()
    for i in range(1,histo.GetNbinsX()+1):
        up=histo.GetBinContent(i)/intUp
        nominal=histoNominal.GetBinContent(i)/intNominal
        newHisto.SetBinContent(i,histoNominal.GetBinContent(i)*nominal/up)
    return newHisto      

  
def smoothTail1D(proj):
    if proj.Integral() == 0:
        print "histogram has zero integral "+proj.GetName()
        return 0
    scale = proj.Integral()
    proj.Scale(1.0/scale)
    
    
    beginFitX = 2100#1500
    endX = 2800
    if period == "2016" or options.output.find("HPHP")!=-1:
        beginFitX=1100
        endX = 1500
    expo=ROOT.TF1("expo","[0]*(1-x/13000.)^[1]/(x/13000)^[2]",2000,8000)
    expo.SetParameters(0,16.,2.)
    expo.SetParLimits(2,1.,20.)
    proj.Fit(expo,"LLMR","",beginFitX,8000)
    beginsmooth = False
    print proj.GetNbinsX()+1
    for j in range(1,proj.GetNbinsX()+1):
        x=proj.GetXaxis().GetBinCenter(j)
        if x>beginFitX:
            if beginsmooth==False:
                if x< endX: #2100: 
                   if abs(proj.GetBinContent(j) - expo.Eval(x)) < 0.00009:# and abs(expo.Derivative(x)- (hist.GetBinContent(j):
                    print beginFitX
                    print "begin smoothing at " +str(x)
                    beginsmooth = True 
               #if abs(proj.GetBinContent(j) - expo.Eval(x)) < 0.00001:# and abs(expo.Derivative(x)- (hist.GetBinContent(j):
                   #print beginFitX
                   #print "begin smoothing at " +str(x)

                   #beginsmooth = True 
                else: beginsmooth = True
            if beginsmooth:
                proj.SetBinContent(j,expo.Eval(x))
    proj.Scale(scale)
    return 1

weights_ = options.weights.split(',')

random=ROOT.TRandom3(101082)

sampleTypes=options.samples.split(',')
period = "2016"
if options.samples.find("HT800")!=-1:
    period = "2017"

stack = ROOT.THStack("stack","")

print "Creating datasets for samples: " ,sampleTypes
dataPlotters=[]
dataPlottersNW=[]
for filename in os.listdir(args[0]):
    for sampleType in sampleTypes:
        if filename.find(sampleType)!=-1:
            fnameParts=filename.split('.')
            fname=fnameParts[0]
            ext=fnameParts[1]
            if ext.find("root") ==-1: continue
            dataPlotters.append(TreePlotter(args[0]+'/'+fname+'.root','AnalysisTree'))
            dataPlotters[-1].setupFromFile(args[0]+'/'+fname+'.pck')
            dataPlotters[-1].addCorrectionFactor('xsec','tree')
            dataPlotters[-1].addCorrectionFactor('genWeight','tree')
            dataPlotters[-1].addCorrectionFactor('puWeight','tree')
            if options.triggerW:
              dataPlotters[-1].addCorrectionFactor('triggerWeight','tree')
              print "Using trigger weights from tree"
            for w in weights_:
	     if w != '': dataPlotters[-1].addCorrectionFactor(w,'branch')
	    corrFactor = 1
            if filename.find('Z') != -1:
                corrFactor = options.corrFactorZ
                print "add correction factor for Z+jets sample"
            if filename.find('W') != -1:
                corrFactor = options.corrFactorW
                print "add correction factor for W+jets sample"
            dataPlotters[-1].addCorrectionFactor(corrFactor,'flat') 

            dataPlotters[-1].filename=fname
            dataPlottersNW.append(TreePlotter(args[0]+'/'+fname+'.root','AnalysisTree'))
            dataPlottersNW[-1].addCorrectionFactor('puWeight','tree')
            dataPlottersNW[-1].addCorrectionFactor('genWeight','tree')
            if options.triggerW: dataPlottersNW[-1].addCorrectionFactor('triggerWeight','tree')
            dataPlottersNW[-1].addCorrectionFactor(corrFactor,'flat')
            for w in weights_: 
             if w != '': dataPlottersNW[-1].addCorrectionFactor(w,'branch')
             if options.triggerW: dataPlottersNW[-1].addCorrectionFactor('triggerWeight','tree')
            dataPlottersNW[-1].addCorrectionFactor(corrFactor,'flat')
            dataPlottersNW[-1].filename=fname
      
data=MergedPlotter(dataPlotters)

fcorr=ROOT.TFile(options.res)

scale = fcorr.Get("scale"+options.resHisto+"Histo")
res   = fcorr.Get("res"  +options.resHisto+"Histo")

binning = getBinning(options.binsMVV,options.minx,options.maxx,options.binsx)
print binning


#distribution of mjet from simulation --> use to validate kernel
mvv_nominal=ROOT.TH1F("mvv_nominal","mvv_nominal",options.binsx,array('f',binning))
mvv_nominal.Sumw2()

mvv_altshapeUp=ROOT.TH1F("mvv_altshapeUp","mvv_altshapeUp",options.binsx,array('f',binning))
mvv_altshapeUp.Sumw2()

mvv_altshape2=ROOT.TH1F("mvv_altshape2","mvv_altshape2",options.binsx,array('f',binning))
mvv_altshape2.Sumw2()

histogram_nominal=ROOT.TH1F("histo_nominal","histo_nominal",options.binsx,array('f',binning))
histogram_nominal.Sumw2()

histogram_altshapeUp=ROOT.TH1F("histo_altshapeUp","histo_altshapeUp",options.binsx,array('f',binning))
histogram_altshapeUp.Sumw2()

histogram_altshape2=ROOT.TH1F("histo_altshape2","histo_altshape2",options.binsx,array('f',binning))
histogram_altshape2.Sumw2()

histograms=[
    histogram_nominal,
    histogram_altshapeUp,
    histogram_altshape2,
    mvv_nominal,
    mvv_altshapeUp,
    mvv_altshape2
  ]

maxEvents = -1
#ok lets populate!
for plotter,plotterNW in zip(dataPlotters,dataPlottersNW):

 #Nominal histogram Pythia8
 c=0
 if plotter.filename.find(sampleTypes[0].replace('.root','')) != -1: 
   print "Preparing nominal histogram for sampletype " ,sampleTypes[0]
   print "filename: ", plotter.filename, " preparing central values histo"
   histI2=plotter.drawTH1Binned('jj_LV_mass',options.cut,"1",array('f',binning))
   canv = ROOT.TCanvas("c1","c1",800,600)
   dataset=plotterNW.makeDataSet('jj_gen_partialMass,jj_l1_gen_pt,jj_l1_gen_softDrop_mass',options.cut,options.firstEv,options.lastEv)     
   
   histTMP=ROOT.TH1F("histoTMP","histo",options.binsx,array('f',binning)) 
   if not(options.usegenmass): 
    datamaker=ROOT.cmg.GaussianSumTemplateMaker1D(dataset,options.var,'jj_l1_gen_pt',scale,res,histTMP)
   else: datamaker=ROOT.cmg.GaussianSumTemplateMaker1D(dataset,options.var,'jj_l1_gen_softDrop_mass',scale,res,histTMP) 

   if histTMP.Integral()>0:
    histTMP.Scale(histI2.Integral()/histTMP.Integral())
    histogram_nominal.Add(histTMP)
    mvv_nominal.Add(histI2)

   histI2.Delete()
   histTMP.Delete()


 if len(sampleTypes)<2: continue
 elif plotter.filename.find(sampleTypes[1].replace('.root','')) != -1: #alternative shape Herwig
   print "Preparing alternative shapes for sampletype " ,sampleTypes[1]
   print "filename: ", plotter.filename, " preparing alternate shape histo"
   histI2=plotter.drawTH1Binned('jj_LV_mass',options.cut,"1",array('f',binning))
 
   dataset=plotterNW.makeDataSet('jj_gen_partialMass,jj_l1_gen_pt,jj_l1_gen_softDrop_mass',options.cut,options.firstEv,options.lastEv)     
   
   histTMP=ROOT.TH1F("histoTMP","histo",options.binsx,array('f',binning))  
   if not(options.usegenmass): 
    datamaker=ROOT.cmg.GaussianSumTemplateMaker1D(dataset,options.var,'jj_l1_gen_pt',scale,res,histTMP)
   else:        
        datamaker=ROOT.cmg.GaussianSumTemplateMaker1D(dataset,options.var,'jj_l1_gen_softDrop_mass',scale,res,histTMP) 

   if histTMP.Integral()>0:
    histTMP.Scale(histI2.Integral()/histTMP.Integral())
    histogram_altshapeUp.Add(histTMP)
    mvv_altshapeUp.Add(histI2)
    
   histI2.Delete()
   histTMP.Delete()

   histogram_altshapeUp.SetLineColor(ROOT.kBlue)
   histogram_altshapeUp.SetFillColorAlpha(ROOT.kBlue, 0.6)
   stack.Add(histogram_altshapeUp)

 if len(sampleTypes)<3: continue
 elif plotter.filename.find(sampleTypes[2].replace('.root','')) != -1: #alternative shape Pythia8+Madgraph (not used for syst but only for cross checks)
   print "Preparing alternative shapes for sampletype " ,sampleTypes[2]
   print "filename: ", plotter.filename, " preparing alternate shape histo"
   
   histI2=plotter.drawTH1Binned('jj_LV_mass',options.cut,"1",array('f',binning))

   dataset=plotterNW.makeDataSet('jj_gen_partialMass,jj_l1_gen_pt,jj_l1_gen_softDrop_mass',options.cut,options.firstEv,options.lastEv)     
   
   histTMP=ROOT.TH1F("histoTMP","histo",options.binsx,array('f',binning))
   if not(options.usegenmass): 
    datamaker=ROOT.cmg.GaussianSumTemplateMaker1D(dataset,options.var,'jj_l1_gen_pt',scale,res,histTMP)
   else: datamaker=ROOT.cmg.GaussianSumTemplateMaker1D(dataset,options.var,'jj_l1_gen_softDrop_mass',scale,res,histTMP) 

   if histTMP.Integral()>0:
    histTMP.Scale(histI2.Integral()/histTMP.Integral())
    histogram_altshape2.Add(histTMP)
    mvv_altshape2.Add(histI2)
    
   histI2.Delete()
   histTMP.Delete()
   histogram_altshape2.SetLineColor(ROOT.kGreen)
   histogram_altshape2.SetFillColorAlpha(ROOT.kGreen, 0.6)
   stack.Add(histogram_altshape2)


print " ********** ALL DONE, now save in output file ", options.output
f=ROOT.TFile(options.output,"RECREATE")
f.cd()
finalHistograms={}
if (options.output).find("Jets")!=-1:
    histograms[0].Add(histograms[1])
    histograms[0].Add(histograms[2])
    
    histograms[3].Add(histograms[4])
    histograms[3].Add(histograms[5])
    print "add all the histograms "
scale = histograms[0].Integral()
scale2 = histograms[3].Integral()
for hist in histograms:
    finalHistograms[hist.GetName()]=hist

if (options.output).find("Jets")!=-1:
    if "histo_altshapeUp" in finalHistograms.keys():    
        finalHistograms["histo_nominal"].Add(finalHistograms["histo_altshapeUp"])
    if "histo_altshape2" in finalHistograms.keys():    
        finalHistograms["histo_nominal"].Add(finalHistograms["histo_altshape2"])
        
    if "mvv_altshapeUp" in finalHistograms.keys():    
        finalHistograms["mvv_nominal"].Add(finalHistograms["mvv_altshapeUp"])
    if "mvv_altshape2" in finalHistograms.keys():    
        finalHistograms["mvv_nominal"].Add(finalHistograms["mvv_altshape2"])
    print "add the histograms for W+jets, Z+jets and ttbar before smoothing the tails"
for hist in finalHistograms.itervalues():
 # hist.Write(hist.GetName()+"_raw")
 if (options.output).find("Jets")!=-1 and hist.GetName()=="histo_nominal":
     print "smooth tails of 1D histogram for vjets background of histo "+hist.GetName()
     if hist.Integral() > 0:
        smoothTail1D(hist)
        if hist.GetName().find("histogram_nominal")!=-1:
            hist.Scale(scale)
        #if hist.GetName().find("mvv_nominal")!=-1:
        #    hist.Scale(scale2)

 hist.Write(hist.GetName())
 finalHistograms[hist.GetName()]=hist

########################################################

options.minx,options.maxx
alpha=1.5/float(options.maxx)
histogram_pt_down,histogram_pt_up=unequalScale(finalHistograms["histo_nominal"],"histo_nominal_PT",alpha)
histogram_pt_down.Write()
histogram_pt_up.Write()

alpha=1.5*float(options.minx)
histogram_opt_down,histogram_opt_up=unequalScale(finalHistograms["histo_nominal"],"histo_nominal_OPT",alpha,-1)
histogram_opt_down.Write()
histogram_opt_up.Write()

alpha=float(options.maxx)*float(options.maxx)
histogram_pt2_down,histogram_pt2_up=unequalScale(finalHistograms["histo_nominal"],"histo_nominal_PT2",alpha,2)
histogram_pt2_down.Write()
histogram_pt2_up.Write()

alpha=float(options.minx)*float(options.minx)
histogram_opt2_down,histogram_opt2_up=unequalScale(finalHistograms["histo_nominal"],"histo_nominal_OPT2",alpha,-2)
histogram_opt2_down.Write()
histogram_opt2_up.Write() 

#################################
c = ROOT.TCanvas("c","C",600,400)
c.SetRightMargin(0.11)
c.SetLeftMargin(0.11)
c.SetTopMargin(0.11)
finalHistograms["histo_nominal"].SetLineColor(ROOT.kBlue)
finalHistograms["histo_nominal"].GetYaxis().SetTitle("arbitrary scale")
finalHistograms["histo_nominal"].GetYaxis().SetTitleOffset(1.5)
finalHistograms["histo_nominal"].GetXaxis().SetTitle("dijet mass")
sf = finalHistograms["histo_nominal"].Integral()
histogram_pt_up     .Scale(sf/histogram_pt_up.Integral())
histogram_pt_down   .Scale(sf/histogram_pt_down.Integral())
histogram_opt_up    .Scale(sf/histogram_opt_up.Integral())
histogram_opt_down  .Scale(sf/histogram_opt_down.Integral())
finalHistograms["histo_nominal"].Draw("hist")
#stack.Draw("histsame")
histogram_pt_up.SetLineColor(ROOT.kRed)
histogram_pt_up.SetLineWidth(2)
histogram_pt_up.Draw("histsame")
histogram_pt_down.SetLineColor(ROOT.kRed)
histogram_pt_down.SetLineWidth(2)
histogram_pt_down.Draw("histsame")
histogram_opt_up.SetLineColor(ROOT.kGreen)
histogram_opt_up.SetLineWidth(2)
histogram_opt_up.Draw("histsame")
histogram_opt_down.SetLineColor(ROOT.kGreen)
histogram_opt_down.SetLineWidth(2)
histogram_opt_down.Draw("histsame")
text = ROOT.TLatex()
text.DrawLatexNDC(0.13,0.92,"#font[62]{CMS} #font[52]{Simulation}")
data = finalHistograms["mvv_nominal"]
data.Scale(sf/data.Integral())
data.SetMarkerColor(ROOT.kBlack)
data.SetMarkerStyle(7)
data.Draw("same")
c.SetLogy()


l = ROOT.TLegend(0.17,0.2,0.6,0.33)
l.AddEntry(data,"simulation","lp")
l.AddEntry(finalHistograms["histo_nominal"],"template","l")
l.AddEntry(histogram_pt_up,"#propto m_{jj}","l")
l.AddEntry(histogram_opt_up,"#propto 1/m_{jj}","l")
l.Draw("same")

tmplabel = "nonRes"
if 'Jets' in sampleTypes: tmplabel="Jets"
if options.output.find('HPLP')!=-1: tmplabel+='_HPLP'
if options.output.find('HPHP')!=-1: tmplabel+='_HPHP'
if 'W' in sampleTypes: tmplabel="W"+tmplabel
if 'Z' in sampleTypes: tmplabel="Z"+tmplabel
c.SaveAs("debug_mVV_kernels_"+tmplabel+".pdf")
print "for debugging save","debug_mVV_kernels_"+tmplabel+".pdf"

########################################################


f.Close()



