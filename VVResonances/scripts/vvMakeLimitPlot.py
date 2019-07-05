#!/usr/bin/env python

#vvMakeLimitPlot.py limits.root -x 1200 -X 5200 -b 0 -s VprimeWV --hvt 1 ---HVTworkspace workspace_JJ_VprimeWV_13TeV.root
import ROOT
import optparse, time, sys, math
from CMGTools.VVResonances.plotting.CMS_lumi import *
from CMGTools.VVResonances.plotting.tdrstyle import *
from array import array
parser = optparse.OptionParser()
parser.add_option("-o","--output",dest="output",default='limitPlot',help="Limit plot")
parser.add_option("-s","--signal",dest="sig",type=str,help="Signal sample",default='BulkGWW')
parser.add_option("--sigscale",dest="sigscale",type=float,help="maximum y",default=0.001)
parser.add_option("-x","--minX",dest="minX",type=float,help="minimum x",default=1000.0)
parser.add_option("-X","--maxX",dest="maxX",type=float,help="maximum x",default=5000.0)
parser.add_option("-y","--minY",dest="minY",type=float,help="minimum y",default=0.0001)
parser.add_option("-Y","--maxY",dest="maxY",type=float,help="maximum y",default=10)
parser.add_option("-b","--blind",dest="blind",type=int,help="Not do observed ",default=1)
parser.add_option("-l","--log",dest="log",type=int,help="Log plot",default=1)

parser.add_option("-t","--titleX",dest="titleX",default='M_{X} [GeV]',help="title of x axis")
parser.add_option("-T","--titleY",dest="titleY",default='#sigma x BR(X #rightarrow WW) [pb]  ',help="title of y axis")

parser.add_option("-p","--period",dest="period",default='2017',help="period")
parser.add_option("-f","--final",dest="final",type=int, default=1,help="Preliminary or not")
parser.add_option("--hvt","--hvt",dest="hvt",type=int, default=0,help="do HVT (1) or do BulkG (2)")
parser.add_option("--HVTworkspace","--HVTworkspace",dest="HVTworkspace",default="workspace_JJ_VprimeWV_13TeV.root",help="HVT workspace with spline interpolation")

#    parser.add_option("-x","--minMVV",dest="minMVV",type=float,help="minimum MVV",default=1000.0)
#    parser.add_option("-X","--maxMVV",dest="maxMVV",type=float,help="maximum MVV",default=13000.0)

(options,args) = parser.parse_args()
#define output dictionary

setTDRStyle()

masses = array('d',[i*100. for i in range(8,60)])
scaleLimits = {}
for m in masses:
 scaleLimits[str(int(m))] = options.sigscale

if options.hvt>0:
 fin = ROOT.TFile.Open(options.HVTworkspace,"READ")
 w = fin.Get("w")
  

# if options.hvt == 1:
#  filenameTHWp = "$CMSSW_BASE/src/CMGTools/VVResonances/scripts/theoryXsec/WprimeWZ.root"
#  filenameTHZp = "$CMSSW_BASE/src/CMGTools/VVResonances/scripts/theoryXsec/ZprimeWW.root"
# else:
#  filenameTHWp = "$CMSSW_BASE/src/CMGTools/VVResonances/scripts/theoryXsec/BulkG.root"
#  filenameTHZp = "$CMSSW_BASE/src/CMGTools/VVResonances/scripts/theoryXsec/BulkG.root"
     
# thFileWp       = ROOT.TFile.Open(filenameTHWp,'READ')   
# thFileZp       = ROOT.TFile.Open(filenameTHZp,'READ')  
# print "Opening file " ,thFileWp.GetName()
# gtheoryWp      = thFileWp.Get("gtheory")
# gtheoryWpUP    = thFileWp.Get("gtheoryUP")
# gtheoryWpDOWN  = thFileWp.Get("gtheoryDOWN")
# gtheoryWpSHADE = thFileWp.Get("grshade")
# print "Opening file " ,thFileZp.GetName()
# gtheoryZp      = thFileZp.Get("gtheory")
# gtheoryZpUP    = thFileZp.Get("gtheoryUP")
# gtheoryZpDOWN  = thFileZp.Get("gtheoryDOWN")
# gtheoryZpSHADE = thFileZp.Get("grshade")

 xsecTot = array('d',[])
 xsecTotUp = array('d',[])
 xsecTotDown = array('d',[])
 shade_x = array('d',[])
 shade_y = array('d',[])

 for m in masses:

  argset = ROOT.RooArgSet()
  MH=w.var("MH")
  argset.add(MH)
  MH.setVal(m)
  if options.hvt == 1:
   func1 = w.function('ZprimeWW_JJ_HPHP_13TeV_2016_sigma')
   func2 = w.function('WprimeWZ_JJ_HPHP_13TeV_2016_sigma')  
  else: 
   func1 = w.function('BulkGravWW_JJ_HPLP_13TeV_2016_sigma')
   func2 = w.function('BulkGravZZ_JJ_HPLP_13TeV_2016_sigma')  
  scaleLimits[str(int(m))] = func1.getVal(argset)+func2.getVal(argset) 

 spline_x_wp = []
 spline_y_wp = []
 spline_y_wpUP = []
 spline_y_wpDOWN = []  
# for i in range(gtheoryWp.GetN()):
 
#  x = ROOT.Double(0.)
#  y = ROOT.Double(0.)
#  gtheoryWp.GetPoint(i,x,y)
#  spline_y_wp.append(y)
#  spline_x_wp.append(x*1000.)
  
#  x = ROOT.Double(0.)
#  y = ROOT.Double(0.)
#  gtheoryWpUP.GetPoint(i,x,y)  
#  spline_y_wpUP.append(y)

#  x = ROOT.Double(0.)
#  y = ROOT.Double(0.)
#  gtheoryWpDOWN.GetPoint(i,x,y)  
#  spline_y_wpDOWN.append(y)
    
 spline_y_zp = [] 
 spline_x_zp = []
 spline_y_zpUP = []
 spline_y_zpDOWN = [] 
# for i in range(gtheoryZp.GetN()):
 
#  x = ROOT.Double(0.)
#  y = ROOT.Double(0.)
#  gtheoryZp.GetPoint(i,x,y)
#  spline_y_zp.append(y)
#  spline_x_zp.append(x*1000.)

#  x = ROOT.Double(0.)
#  y = ROOT.Double(0.)
#  gtheoryZpUP.GetPoint(i,x,y)  
#  spline_y_zpUP.append(y)

#  x = ROOT.Double(0.)
#  y = ROOT.Double(0.)
#  gtheoryZpDOWN.GetPoint(i,x,y)  
#  spline_y_zpDOWN.append(y)
       
# spline_zp=ROOT.RooSpline1D("Zprime_sigma","Zprime_sigma",MH,len(spline_x_zp),array('d',spline_x_zp),array('d',spline_y_zp))  
# spline_wp=ROOT.RooSpline1D("Wprime_sigma","Wprime_sigma",MH,len(spline_x_wp),array('d',spline_x_wp),array('d',spline_y_wp))  
# spline_zpUP=ROOT.RooSpline1D("Zprime_sigmaUP","Zprime_sigmaUP",MH,len(spline_x_zp),array('d',spline_x_zp),array('d',spline_y_zpUP))  
# spline_wpUP=ROOT.RooSpline1D("Wprime_sigmaUP","Wprime_sigmaUP",MH,len(spline_x_wp),array('d',spline_x_wp),array('d',spline_y_wpUP)) 
# spline_zpDOWN=ROOT.RooSpline1D("Zprime_sigmaDOWN","Zprime_sigmaDOWN",MH,len(spline_x_zp),array('d',spline_x_zp),array('d',spline_y_zpDOWN))  
# spline_wpDOWN=ROOT.RooSpline1D("Wprime_sigmaDOWN","Wprime_sigmaDOWN",MH,len(spline_x_wp),array('d',spline_x_wp),array('d',spline_y_wpDOWN)) 
  
 for m in masses:
 
  MH.setVal(m)
  
#  tot = spline_zp.getVal(argset)+spline_wp.getVal(argset)
#  xsecTot.append(tot)
  
#  uncUp_zp = spline_zpUP.getVal(argset)-spline_zp.getVal(argset)
#  uncUp_wp = spline_wpUP.getVal(argset)-spline_wp.getVal(argset)
#  xsecTotUp.append( tot+math.sqrt(uncUp_zp*uncUp_zp+uncUp_wp*uncUp_wp) )
  
#  uncDown_zp = spline_zp.getVal(argset)-spline_zpDOWN.getVal(argset)
#  uncDown_wp = spline_wp.getVal(argset)-spline_wpDOWN.getVal(argset)
#  xsecTotDown.append( tot-math.sqrt(uncDown_zp*uncDown_zp+uncDown_wp*uncDown_wp) )
    
  shade_x.append(m)
#  shade_y.append( tot+math.sqrt(uncUp_zp*uncUp_zp+uncUp_wp*uncUp_wp) )

 for i in range( len(masses)-1, -1, -1 ):
  shade_x.append(masses[i])
#  shade_y.append(xsecTotDown[i])
       
# gtheory = ROOT.TGraphErrors(len(masses),masses,xsecTot)
# gtheory.SetLineColor(ROOT.kRed)
# gtheory.SetLineWidth(3)
# gtheoryUP = ROOT.TGraphErrors(len(masses),masses,xsecTotUp)
# gtheoryUP.SetLineColor(ROOT.kRed-2)
# gtheoryUP.SetLineWidth(3)
# gtheoryDOWN = ROOT.TGraphErrors(len(masses),masses,xsecTotDown)
# gtheoryDOWN.SetLineColor(ROOT.kRed-2)
# gtheoryDOWN.SetLineWidth(3)
# gtheorySHADE = ROOT.TGraphErrors(len(shade_x),shade_x,shade_y)
# gtheorySHADE.SetLineColor(ROOT.kRed-2)
# gtheorySHADE.SetLineWidth(3)
 
f=ROOT.TFile(args[0])
limit=f.Get("limit")
data={}

def rescaleaxis(g,scale=1000.):
    N = g.GetN()
    x = g.GetX()
    for i in range(N):
        x[i] *= scale
    g.GetHistogram().Delete()
    g.SetHistogram(0)
    return
    
for event in limit:
    if float(event.mh)<options.minX or float(event.mh)>options.maxX:
        continue
    
    if not (event.mh in data.keys()):
        data[event.mh]={}

    lim = event.limit*scaleLimits[str(int(event.mh))]
    if event.quantileExpected<0:            
        data[event.mh]['obs']=lim
    if event.quantileExpected>0.02 and event.quantileExpected<0.03:            
        data[event.mh]['-2sigma']=lim
    if event.quantileExpected>0.15 and event.quantileExpected<0.17:            
        data[event.mh]['-1sigma']=lim
    if event.quantileExpected>0.49 and event.quantileExpected<0.51:            
        data[event.mh]['exp']=lim
    if event.quantileExpected>0.83 and event.quantileExpected<0.85:            
        data[event.mh]['+1sigma']=lim
    if event.quantileExpected>0.974 and event.quantileExpected<0.976:            
        data[event.mh]['+2sigma']=lim

mean=ROOT.TGraphAsymmErrors()
mean.SetName("mean")

band68=ROOT.TGraphAsymmErrors()
band68.SetName("band68")
band95=ROOT.TGraphAsymmErrors()
band95.SetName("band95")
bandObs=ROOT.TGraph()
bandObs.SetName("bandObs")

line_plus1=ROOT.TGraph()
line_plus1.SetName("line_plus1")

line_plus2=ROOT.TGraph()
line_plus2.SetName("line_plus2")

line_minus1=ROOT.TGraph()
line_minus1.SetName("line_minus1")

line_minus2=ROOT.TGraph()
line_minus2.SetName("line_minus2")



N=0
for mass,info in data.iteritems():
    print 'Setting mass',mass,info

    if not ('exp' in info.keys() and '+1sigma' in info.keys() and '+2sigma' in info.keys() and '-1sigma' in info.keys() and '-2sigma' in info.keys()):
        print 'Incomplete file'
        continue
    if options.blind==0 and not ('obs' in info.keys()):
        print 'Incomplete file'
        continue    
    
    mean.SetPoint(N,mass,info['exp'])    
    band68.SetPoint(N,mass,info['exp'])
    band95.SetPoint(N,mass,info['exp'])
    line_plus1.SetPoint(N,mass,info['+1sigma'])
    line_plus2.SetPoint(N,mass,info['+2sigma'])
    line_minus1.SetPoint(N,mass,info['-1sigma'])
    line_minus2.SetPoint(N,mass,info['-2sigma'])

    if options.blind==0: bandObs.SetPoint(N,mass,info['obs'])
    band68.SetPointError(N,0.0,0.0,info['exp']-info['-1sigma'],info['+1sigma']-info['exp'])
    band95.SetPointError(N,0.0,0.0,info['exp']-info['-2sigma'],info['+2sigma']-info['exp'])
    N=N+1


mean.Sort()
band68.Sort()
band95.Sort()
if options.blind==0: bandObs.Sort()
line_plus1.Sort()    
line_plus2.Sort()    
line_minus1.Sort()    
line_minus2.Sort()    




band68.SetFillColor(ROOT.kGreen)
band68.SetLineWidth(3)
band68.SetLineColor(ROOT.kWhite)
band68.SetLineStyle(0)
band68.SetMarkerStyle(0)

band95.SetFillColor(ROOT.kYellow)
band95.SetLineColor(ROOT.kWhite)

bandObs.SetLineWidth(3)
bandObs.SetLineColor(ROOT.kBlack)
bandObs.SetMarkerStyle(20)

mean.SetLineWidth(2)
mean.SetLineColor(ROOT.kBlack)
mean.SetLineStyle(2)

line_plus1.SetLineWidth(1)
line_plus1.SetLineColor(ROOT.kGreen+1)

line_plus2.SetLineWidth(1)
line_plus2.SetLineColor(ROOT.kOrange-2)

line_minus1.SetLineWidth(1)
line_minus1.SetLineColor(ROOT.kGreen+1)

line_minus2.SetLineWidth(1)
line_minus2.SetLineColor(ROOT.kOrange-2)

if not options.hvt:
 '''
 gtheory = ROOT.TGraphErrors(1)
 gtheory.SetLineColor(ROOT.kRed)
 gtheory.SetLineWidth(3)
 gtheoryUP = ROOT.TGraphErrors(1)
 gtheoryUP.SetLineColor(ROOT.kRed-2)
 gtheoryUP.SetLineWidth(3)
 gtheoryDOWN = ROOT.TGraphErrors(1)
 gtheoryDOWN.SetLineColor(ROOT.kRed-2)
 gtheoryDOWN.SetLineWidth(3)
 gtheorySHADE = ROOT.TGraphErrors(1)
 gtheorySHADE.SetLineColor(ROOT.kRed-2)
 gtheorySHADE.SetLineWidth(3)
 

 filenameTH = "$CMSSW_BASE/src/CMGTools/VVResonances/scripts/theoryXsec/%s.root"%options.sig
 thFile       = ROOT.TFile.Open(filenameTH,'READ')   
 print "Opening file " ,thFile.GetName()
 gtheory      = thFile.Get("gtheory")
 gtheoryUP    = thFile.Get("gtheoryUP")
 gtheoryDOWN  = thFile.Get("gtheoryDOWN")
 gtheorySHADE = thFile.Get("grshade")
 
 rescaleaxis(gtheory	)
 rescaleaxis(gtheoryUP	)
 rescaleaxis(gtheoryDOWN )
 rescaleaxis(gtheorySHADE)
 '''
'''
gtheory     .SetName("%s_gtheory"    %options.sig)
gtheoryUP   .SetName("%s_gtheoryUP"  %options.sig)
gtheoryDOWN .SetName("%s_gtheoryDOWN"%options.sig)
gtheorySHADE.SetName("%s_grshade"    %options.sig)
gtheorySHADE.SetLineColor(0)
gtheorySHADE.SetFillColor(ROOT.kRed)
gtheorySHADE.SetFillStyle(3013)
gtheoryUP.SetLineColor(ROOT.kRed)
gtheoryDOWN.SetLineColor(ROOT.kRed)
gtheoryUP.SetLineWidth(1)
gtheoryDOWN.SetLineWidth(1)
# thFile.Close()
'''
#plotting information
H_ref = 600; 
W_ref = 800; 
W = W_ref
H = H_ref

T = 0.08*H_ref
B = 0.12*H_ref 
L = 0.12*W_ref
R = 0.04*W_ref
c=ROOT.TCanvas("c","c",50,50,W,H)
c.SetFillColor(0)
c.SetBorderMode(0)
c.SetFrameFillStyle(0)
c.SetFrameBorderMode(0)
c.SetLeftMargin( L/W )
c.SetRightMargin( R/W )
c.SetTopMargin( T/H )
c.SetBottomMargin( B/H )
c.SetTickx(0)
c.SetTicky(0)
c.GetWindowHeight()
c.GetWindowWidth()
c.SetLogy()
# c.SetGrid()
c.SetLogy()
c.cd()


if "Wprime"  in options.sig: 
  ltheory="#sigma_{TH}#timesBR(W'#rightarrowWZ) HVT_{B}"
  ytitle ="#sigma x BR(W' #rightarrow WZ) [pb]  "
  xtitle = "M_{W'} [GeV]"
if "BulkGWW" in options.sig: 
  ltheory="#sigma_{TH}#timesBR(G_{Bulk}#rightarrowWW) #tilde{k}=0.5"
  ytitle ="#sigma x BR(G_{Bulk} #rightarrow WW) [pb]  "
  xtitle = "M_{G_{Bulk}} [GeV]"
if "BulkGZZ" in options.sig: 
  ltheory="#sigma_{TH}#timesBR(G_{Bulk}#rightarrowZZ) #tilde{k}=0.5"  
  ytitle ="#sigma x BR(G_{Bulk} #rightarrow ZZ) [pb]  "
  xtitle = "M_{G_{Bulk}} [GeV]"
if "Zprime"  in options.sig: 
  ltheory="#sigma_{TH}#timesBR(Z'#rightarrowWW) HVT_{B}"
  ytitle ="#sigma x BR(Z' #rightarrow WW) [pb]  "
  xtitle = "M_{Z'} [GeV]"
if "Vprime"  in options.sig: 
  ltheory="#sigma_{TH}#timesBR(V'#rightarrowWV) HVT_{B}"
  ytitle ="#sigma x BR(V' #rightarrow WV) [pb]  "
  xtitle = "M_{V'} [GeV]"
if "BulkGVV"  in options.sig: 
  ltheory="#sigma_{TH}#timesBR(G_{Bulk}#rightarrowVV) #tilde{k}=0.5"  
  ytitle ="#sigma x BR(G_{Bulk} #rightarrow VV) [pb]  "
  xtitle = "M_{G_{Bulk}} [GeV]"
    
frame=c.DrawFrame(options.minX,options.minY,options.maxX,options.maxY)
frame.GetXaxis().SetTitle(xtitle)
frame.GetXaxis().SetTitleOffset(0.9)
frame.GetXaxis().SetTitleSize(0.05)

frame.GetYaxis().SetTitle(ytitle)
frame.GetYaxis().SetTitleSize(0.05)
frame.GetYaxis().SetTitleOffset(1.15)

frame.Draw()
band95.Draw("3same")
band68.Draw("3same")
# band68.Draw("XLsame")
line_plus1.Draw("Lsame")
line_plus2.Draw("Lsame")
line_minus1.Draw("Lsame")
line_minus2.Draw("Lsame")
mean.Draw("Lsame")
#gtheory.Draw("Lsame")
#gtheorySHADE.Draw("Fsame")

c.SetLogy(options.log)
c.Draw()



leg  = ROOT.TLegend(0.498995,0.6602591,0.9446734,0.9011917)
leg2 = ROOT.TLegend(0.498995,0.6602591,0.9446734,0.9011917)
leg.SetTextSize(0.028)
leg.SetLineColor(1)
leg.SetShadowColor(0)
leg.SetLineStyle(1)
leg.SetLineWidth(1)
leg.SetFillColor(ROOT.kWhite)
# leg.SetFillStyle(0)
leg.SetMargin(0.35)
leg2.SetTextSize(0.028)
leg2.SetLineColor(1)
leg2.SetShadowColor(0)
leg2.SetLineStyle(1)
leg2.SetLineWidth(1)
leg2.SetFillColor(0)
leg2.SetFillStyle(0)
leg2.SetMargin(0.35)
leg.SetBorderSize(1)

if options.blind==0: leg.AddEntry(bandObs, "Observed", "Lp")
leg.AddEntry(band68, "Expected #pm 1 std. deviation", "f")
leg.AddEntry(band95 , "Expected #pm 2 std. deviation", "f")
#leg.AddEntry(gtheory, ltheory, "L")

if not options.blind: leg2.AddEntry(bandObs, " ", "")
leg2.AddEntry(mean, " ", "L")
leg2.AddEntry(mean, " ", "L")
#leg2.AddEntry(gtheory, " ", "")      
      




if options.final:
    cmslabel_final(c,options.period,11)
else:
    cmslabel_prelim(c,options.period,11)
leg.Draw()
leg2.Draw()
c.Update()
c.RedrawAxis()

if options.blind==0:
    bandObs.Draw("PLsame")
c.SaveAs(options.output+options.sig+".png")    
c.SaveAs(options.output+options.sig+".pdf")    
c.SaveAs(options.output+options.sig+".C")    

fout=ROOT.TFile(options.output+options.sig+".root","RECREATE")
fout.cd()
c.Write()
band68.Write()
band95.Write()
bandObs.Write()
line_plus1.Write()    
line_plus2.Write()    
line_minus1.Write()    
line_minus2.Write()    

fout.Close()
f.Close()
