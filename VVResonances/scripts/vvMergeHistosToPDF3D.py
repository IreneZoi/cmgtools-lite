#!/usr/bin/env python

import ROOT
from array import array
import os, sys, re, optparse,pickle,shutil,json, math

def getListOfBinsLowEdge(hist,dim):
    axis =0
    N = 0
    if dim =="x":
        axis= hist.GetXaxis()
        N = hist.GetNbinsX()
    if dim =="y":
        axis = hist.GetYaxis()
        N = hist.GetNbinsY()
    if dim =="z":
        axis = hist.GetZaxis()
        N = hist.GetNbinsZ()
    if axis==0:
        return {}
    
    mmin = axis.GetXmin()
    mmax = axis.GetXmax()
    r =[]
    for i in range(1,N+2):
        #v = mmin + i * (mmax-mmin)/float(N)
        r.append(axis.GetBinLowEdge(i))
    return r
    
def getBinning(histo):
    b = []
    for i in range(1,histo.GetNbinsX()+2):
        b.append(histo.GetBinLowEdge(i))
    return array('d',b)

def makeHisto(name,fx,nhistox,fy,nhistoy,fz,nhistoz,fout):
    histox=fx.Get(nhistox)
    histoy=fy.Get(nhistoy)
    histoz=fz.Get(nhistoz)
    try:
        if histox == None: 
            print "histox == None"
            raise TypeError
        if histoy == None: 
            print "histoy == None"
            raise TypeError
        if histoz == None:
            print "histoz == None"
            raise TypeError
    except TypeError: return -1
    hxx = histox.ProjectionX("projX")
    hxy = histox.ProjectionY("projY")
    binningx= getBinning(hxx)
    binningz= getBinning(hxy)
    binningz2= getBinning(histoz)
    binningz3= getBinning(histoy.ProjectionY("projY2"))
    h=ROOT.TH3F(name,name,len(binningx)-1,binningx,len(binningx)-1,binningx,len(binningz)-1,binningz)
    for k in range(1,histoz.GetNbinsX()+2):
     for j in range(1,histoy.GetNbinsX()+2):
        for i in range(1,histox.GetNbinsX()+2):
            c = histoz.GetBinContent(k)*histox.GetBinContent(i,k)*histoy.GetBinContent(j,k)
            #if c == 0:
            #    print histoz.GetBinCenter(k)
            h.SetBinContent(i,j,k,c)
    fout.cd()
    h.Write()

def add_sigmoid_shape(h):

 xbins = array("f",getListOfBinsLowEdge(h,"x"))
 zbins = array("f",getListOfBinsLowEdge(h,"z"))

 hname = h.GetName()
 hnew_up = ROOT.TH3F(hname+'_OPT3Up',hname+'_OPT3Up',len(xbins)-1,xbins,len(xbins)-1,xbins,len(zbins)-1,zbins)
 hnew_down = ROOT.TH3F(hname+'_OPT3Down',hname+'_OPT3Down',len(xbins)-1,xbins,len(xbins)-1,xbins,len(zbins)-1,zbins)

 for bx in range(1,h.GetNbinsX()+1):
  for by in range(1,h.GetNbinsY()+1):
   for bz in range(1,h.GetNbinsZ()+1):
   
    x = h.GetXaxis().GetBinCenter(bx)
    y = h.GetYaxis().GetBinCenter(by)
    z = h.GetZaxis().GetBinCenter(bz)
    n_in = h.GetBinContent(bx,by,bz)
    
    factor = 1./(1+math.exp( ((x-235.)/10.)+((y-235.)/10.)+(-(z-2000.)/100.)) )
   
    hnew_up.SetBinContent(bx,by,bz,n_in*factor) 
    hnew_down.SetBinContent(bx,by,bz,n_in/(factor))    
    
 hnew_up.Scale(1./hnew_up.Integral())
 hnew_down.Scale(1./hnew_down.Integral()) 
 
 return hnew_up,hnew_down
 
parser = optparse.OptionParser()
parser.add_option("-s","--systX",dest="systX",default='',help="Comma   separated and semicolon separated systs for p0 ")
parser.add_option("-S","--systY",dest="systY",default='',help="Comma   separated and semicolon separated systs for p1 ")
parser.add_option("-C","--systCommon",dest="systCommon",default='',help="Comma   separated and semicolon separated systs for p2")
parser.add_option("-i","--inputX",dest="inputX",default='erfexp',help="Comma   separated and semicolon separated systs for p2")
parser.add_option("-I","--inputY",dest="inputY",default='erfexp',help="Comma   separated and semicolon separated systs for p2")
parser.add_option("-z","--inputZ",dest="inputZ",default='erfexp',help="Comma   separated and semicolon separated systs for p2")
parser.add_option("-o","--output",dest="output",help="Output ROOT File",default='')


(options,args) = parser.parse_args()
#define output dictionary

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

inputx=ROOT.TFile(options.inputX)
inputy=ROOT.TFile(options.inputY)
inputz=ROOT.TFile(options.inputZ)
output=ROOT.TFile(options.output,"RECREATE")

print options.inputX
print options.inputY
print options.inputZ
print options.output

print "   - nominal"
makeHisto("histo",inputx,"histo_nominal",inputy,"histo_nominal",inputz,"histo_nominal",output)
print "   - pt (x,y) up/down"
makeHisto("histo_PTXYUp",inputx,"histo_nominal_PTUp",inputy,"histo_nominal_PTUp",inputz,"histo_nominal",output)
makeHisto("histo_PTXYDown",inputx,"histo_nominal_PTDown",inputy,"histo_nominal_PTDown",inputz,"histo_nominal",output)
print "   - pt z up/down"
makeHisto("histo_PTZUp",inputx,"histo_nominal",inputy,"histo_nominal",inputz,"histo_nominal_PTUp",output)
makeHisto("histo_PTZDown",inputx,"histo_nominal",inputy,"histo_nominal",inputz,"histo_nominal_PTDown",output)
print "   - pt x up/down"
makeHisto("histo_PTXUp",inputx,"histo_nominal_PTUp",inputy,"histo_nominal",inputz,"histo_nominal",output)
makeHisto("histo_PTXDown",inputx,"histo_nominal_PTDown",inputy,"histo_nominal",inputz,"histo_nominal",output)
print "   - pt y up/down"
makeHisto("histo_PTYUp",inputx,"histo_nominal",inputy,"histo_nominal_PTUp",inputz,"histo_nominal",output)
makeHisto("histo_PTYDown",inputx,"histo_nominal",inputy,"histo_nominal_PTUp",inputz,"histo_nominal",output)
print "   - Opt (x,y) up/down"
makeHisto("histo_OPTXYUp",inputx,"histo_nominal_OPTUp",inputy,"histo_nominal_OPTUp",inputz,"histo_nominal",output)
makeHisto("histo_OPTXYDown",inputx,"histo_nominal_OPTDown",inputy,"histo_nominal_OPTDown",inputz,"histo_nominal",output)
print "   - Opt z up/down"
makeHisto("histo_OPTZUp",inputx,"histo_nominal",inputy,"histo_nominal",inputz,"histo_nominal_OPTUp",output)
makeHisto("histo_OPTZDown",inputx,"histo_nominal",inputy,"histo_nominal",inputz,"histo_nominal_OPTDown",output)
print "   - pt x,y,z up/down"
makeHisto("histo_PTUp",inputx,"histo_nominal_PTUp",inputy,"histo_nominal_PTUp",inputz,"histo_nominal_PTUp",output)
makeHisto("histo_PTDown",inputx,"histo_nominal_PTDown",inputy,"histo_nominal_PTDown",inputz,"histo_nominal_PTDown",output)
print "   - Opt x,y,z up/down"
makeHisto("histo_OPTUp",inputx,"histo_nominal_OPTUp",inputy,"histo_nominal_OPTUp",inputz,"histo_nominal_OPTUp",output)
makeHisto("histo_OPTDown",inputx,"histo_nominal_OPTDown",inputy,"histo_nominal_OPTDown",inputz,"histo_nominal_OPTDown",output)
print "   - Opt x up/down"
makeHisto("histo_OPTXUp",inputx,"histo_nominal_OPTUp",inputy,"histo_nominal",inputz,"histo_nominal",output)
makeHisto("histo_OPTXDown",inputx,"histo_nominal_OPTDown",inputy,"histo_nominal",inputz,"histo_nominal",output)
print "   - Opt y up/down"
makeHisto("histo_OPTYUp",inputx,"histo_nominal",inputy,"histo_nominal_OPTUp",inputz,"histo_nominal",output)
makeHisto("histo_OPTYDown",inputx,"histo_nominal",inputy,"histo_nominal_OPTUp",inputz,"histo_nominal",output)
print "   - Opt2 x,y,z up/down"
makeHisto("histo_OPT2Up",inputx,"histo_nominal_OPT2Up",inputy,"histo_nominal_OPT2Up",inputz,"histo_nominal_OPT2Up",output)
makeHisto("histo_OPT2Down",inputx,"histo_nominal_OPT2Down",inputy,"histo_nominal_OPT2Down",inputz,"histo_nominal_OPT2Down",output)

print "   - herwig z"
makeHisto("histo_altshapeZUp",inputx,"histo_nominal",inputy,"histo_nominal",inputz,"histo_altshapeUp",output)
makeHisto("histo_altshapeZDown",inputx,"histo_nominal",inputy,"histo_nominal",inputz,"histo_altshapeDown",output)
print "   - madgraph z"
makeHisto("histo_altshape2ZUp",inputx,"histo_nominal",inputy,"histo_nominal",inputz,"histo_altshape2Up",output)
makeHisto("histo_altshape2ZDown",inputx,"histo_nominal",inputy,"histo_nominal",inputz,"histo_altshape2Down",output)
print "   - powheg z"
makeHisto("histo_altshape3ZUp",inputx,"histo_nominal",inputy,"histo_nominal",inputz,"histo_NLO",output)
makeHisto("histo_altshape3ZDown",inputx,"histo_nominal",inputy,"histo_nominal",inputz,"histo_NLODown",output)
print "   - herwig pdf"
makeHisto("histo_altshapeUp",inputx,"histo_altshapeUp",inputy,"histo_altshapeUp",inputz,"histo_altshapeUp",output)
makeHisto("histo_altshapeDown",inputx,"histo_altshapeDown",inputy,"histo_altshapeDown",inputz,"histo_altshapeDown",output)
print "   - madgraph pdf"
makeHisto("histo_altshape2Up",inputx,"histo_altshape2Up",inputy,"histo_altshape2Up",inputz,"histo_altshape2Up",output)
makeHisto("histo_altshape2Down",inputx,"histo_altshape2Down",inputy,"histo_altshape2Down",inputz,"histo_altshape2Down",output)
print "   - powheg pdf"
makeHisto("histo_altshape3Up",inputx,"histo_NLO",inputy,"histo_NLO",inputz,"histo_NLO",output)
makeHisto("histo_altshape3Down",inputx,"histo_NLODown",inputy,"histo_NLODown",inputz,"histo_NLODown",output)

#print "   - Scale (x,y,z) up/down"
#makeHisto("histo_ScaleUp",inputx,"histo_nominal_ScaleXUp",inputy,"histo_nominal_ScaleXUp",inputz,"histo_nominal_scale_up",output)
#makeHisto("histo_ScaleDown",inputx,"histo_nominal_ScaleXDown",inputy,"histo_nominal_ScaleXDown",inputz,"histo_nominal_scale_down",output)
#print "   - Res (x,y,z) up/down"
#makeHisto("histo_ResUp",inputx,"histo_nominal_ResXUp",inputy,"histo_nominal_ResXUp",inputz,"histo_nominal_res_up",output)
#makeHisto("histo_ResDown",inputx,"histo_nominal_ResXDown",inputy,"histo_nominal_ResXDown",inputz,"histo_nominal_res_down",output)

print "   - nominal (herwig)"
makeHisto("histo_herwig",inputx,"histo_altshapeUp",inputy,"histo_altshapeUp",inputz,"histo_altshapeUp",output)
print "   - pt (x,y) up/down"
makeHisto("histo_herwig_PTXYUp",inputx,"histo_altshapeUp_PTUp",inputy,"histo_altshapeUp_PTUp",inputz,"histo_altshapeUp",output)
makeHisto("histo_herwig_PTXYDown",inputx,"histo_altshapeUp_PTDown",inputy,"histo_altshapeUp_PTDown",inputz,"histo_altshapeUp",output)
print "   - pt z up/down"
makeHisto("histo_herwig_PTZUp",inputx,"histo_altshapeUp",inputy,"histo_altshapeUp",inputz,"histo_altshapeUp_PTUp",output)
makeHisto("histo_herwig_PTZDown",inputx,"histo_altshapeUp",inputy,"histo_altshapeUp",inputz,"histo_altshapeUp_PTDown",output)
print "   - Opt (x,y) up/down"
makeHisto("histo_herwig_OPTXYUp",inputx,"histo_altshapeUp_OPTUp",inputy,"histo_altshapeUp_OPTUp",inputz,"histo_altshapeUp",output)
makeHisto("histo_herwig_OPTXYDown",inputx,"histo_altshapeUp_OPTDown",inputy,"histo_altshapeUp_OPTDown",inputz,"histo_altshapeUp",output)
print "   - Opt z up/down"
makeHisto("histo_herwig_OPTZUp",inputx,"histo_altshapeUp",inputy,"histo_altshapeUp",inputz,"histo_altshapeUp_OPTUp",output)
makeHisto("histo_herwig_OPTZDown",inputx,"histo_altshapeUp",inputy,"histo_altshapeUp",inputz,"histo_altshapeUp_OPTDown",output)
print "   - pt x,y,z up/down"
makeHisto("histo_herwig_PTUp",inputx,"histo_altshapeUp_PTUp",inputy,"histo_altshapeUp_PTUp",inputz,"histo_altshapeUp_PTUp",output)
makeHisto("histo_herwig_PTDown",inputx,"histo_altshapeUp_PTDown",inputy,"histo_altshapeUp_PTDown",inputz,"histo_altshapeUp_PTDown",output)
print "   - Opt x,y,z up/down"
makeHisto("histo_herwig_OPTUp",inputx,"histo_altshapeUp_OPTUp",inputy,"histo_altshapeUp_OPTUp",inputz,"histo_altshapeUp_OPTUp",output)
makeHisto("histo_herwig_OPTDown",inputx,"histo_altshapeUp_OPTDown",inputy,"histo_altshapeUp_OPTDown",inputz,"histo_altshapeUp_OPTDown",output)
print "   - pythia pdf 4 herwig"
makeHisto("histo_herwig_altshapeUp",inputx,"histo_nominal",inputy,"histo_nominal",inputz,"histo_nominal",output)
makeHisto("histo_herwig_altshapeDown",inputx,"histo_altshapeUp_altshapeDown",inputy,"histo_altshapeUp_altshapeDown",inputz,"histo_altshapeUp_altshapeDown",output)
print "   - madgraph pdf 4 herwig"
makeHisto("histo_herwig_altshape2Up",inputx,"histo_altshape2Up",inputy,"histo_altshape2Up",inputz,"histo_altshape2Up",output)
makeHisto("histo_herwig_altshape2Down",inputx,"histo_altshapeUp_altshape2Down",inputy,"histo_altshapeUp_altshape2Down",inputz,"histo_altshapeUp_altshape2Down",output)
print "   - powheg pdf 4 herwig"
makeHisto("histo_herwig_altshape3Up",inputx,"histo_NLO",inputy,"histo_NLO",inputz,"histo_NLO",output)
makeHisto("histo_herwig_altshape3Down",inputx,"histo_altshapeUp_NLODown",inputy,"histo_altshapeUp_NLODown",inputz,"histo_altshapeUp_NLODown",output)

print "   - nominal (madgraph)"
makeHisto("histo_madgraph",inputx,"histo_altshape2Up",inputy,"histo_altshape2Up",inputz,"histo_altshape2Up",output)
print "   - pt (x,y) up/down"
makeHisto("histo_madgraph_PTXYUp",inputx,"histo_altshape2_PTUp",inputy,"histo_altshape2_PTUp",inputz,"histo_altshape2Up",output)
makeHisto("histo_madgraph_PTXYDown",inputx,"histo_altshape2_PTDown",inputy,"histo_altshape2_PTDown",inputz,"histo_altshape2Up",output)
print "   - pt z up/down"
makeHisto("histo_madgraph_PTZUp",inputx,"histo_altshape2Up",inputy,"histo_altshape2Up",inputz,"histo_altshape2_PTUp",output)
makeHisto("histo_madgraph_PTZDown",inputx,"histo_altshape2Up",inputy,"histo_altshape2Up",inputz,"histo_altshape2_PTDown",output)
print "   - Opt (x,y) up/down"
makeHisto("histo_madgraph_OPTXYUp",inputx,"histo_altshape2_OPTUp",inputy,"histo_altshape2_OPTUp",inputz,"histo_altshape2Up",output)
makeHisto("histo_madgraph_OPTXYDown",inputx,"histo_altshape2_OPTDown",inputy,"histo_altshape2_OPTDown",inputz,"histo_altshape2Up",output)
print "   - Opt z up/down"
makeHisto("histo_madgraph_OPTZUp",inputx,"histo_altshape2Up",inputy,"histo_altshape2Up",inputz,"histo_altshape2_OPTUp",output)
makeHisto("histo_madgraph_OPTZDown",inputx,"histo_altshape2Up",inputy,"histo_altshape2Up",inputz,"histo_altshape2_OPTDown",output)
print "   - pt x,y,z up/down"
makeHisto("histo_madgraph_PTUp",inputx,"histo_altshape2_PTUp",inputy,"histo_altshape2_PTUp",inputz,"histo_altshape2_PTUp",output)
makeHisto("histo_madgraph_PTDown",inputx,"histo_altshape2_PTDown",inputy,"histo_altshape2_PTDown",inputz,"histo_altshape2_PTDown",output)
print "   - Opt x,y,z up/down"
makeHisto("histo_madgraph_OPTUp",inputx,"histo_altshape2_OPTUp",inputy,"histo_altshape2_OPTUp",inputz,"histo_altshape2_OPTUp",output)
makeHisto("histo_madgraph_OPTDown",inputx,"histo_altshape2_OPTDown",inputy,"histo_altshape2_OPTDown",inputz,"histo_altshape2_OPTDown",output)
print "   - pythia pdf 4 madgraph"
makeHisto("histo_madgraph_altshapeUp",inputx,"histo_nominal",inputy,"histo_nominal",inputz,"histo_nominal",output)
makeHisto("histo_madgraph_altshapeDown",inputx,"histo_altshape2_altshapeDown",inputy,"histo_altshape2_altshapeDown",inputz,"histo_altshape2_altshapeDown",output)
print "   - herwig pdf 4 madgraph"
makeHisto("histo_madgraph_altshape2Up",inputx,"histo_altshapeUp",inputy,"histo_altshapeUp",inputz,"histo_altshapeUp",output)
makeHisto("histo_madgraph_altshape2Down",inputx,"histo_altshape2_altshape2Down",inputy,"histo_altshape2_altshape2Down",inputz,"histo_altshape2_altshape2Down",output)
print "   - powheg pdf 4 madgraph"
makeHisto("histo_madgraph_altshape3Up",inputx,"histo_NLO",inputy,"histo_NLO",inputz,"histo_NLO",output)
makeHisto("histo_madgraph_altshape3Down",inputx,"histo_altshape2_NLODown",inputy,"histo_altshape2_NLODown",inputz,"histo_altshape2_NLODown",output)

print "Write file "+options.output

inputx.Close()
inputy.Close()
inputz.Close()
output.Close()

print "Update out file with sigmoid shapes:"
print

output=ROOT.TFile(options.output,"UPDATE")
print "- sigmoid for nominal"
histo_up,histo_down = add_sigmoid_shape(output.histo)
histo_up.Write()
histo_down.Write()
try: 
  print "- sigmoid for herwig"
  histo_herwig_up,histo_herwig_down = add_sigmoid_shape(output.histo_herwig)
  histo_herwig_up.Write()
  histo_herwig_down.Write()
except:
    print('histo_herwig not in output file!')
try:
 print "- sigmoid for madgraph"
 histo_madgraph_up,histo_madgraph_down = add_sigmoid_shape(output.histo_madgraph)
 histo_madgraph_up.Write()
 histo_madgraph_down.Write()
except:
     print('histo_madgraph not in output file!')

output.Close()
