#python runFitPlots_vjets_signal_oneyear_splitRes.py -n workspace_JJ_ZprimeZH_HPLP_13TeV_2016.root  -l sigonly -i 2016/JJ_nonRes_HPLP.root -M 2000 -s
import ROOT
ROOT.gROOT.SetBatch(True)
import os, sys, re, optparse,pickle,shutil,json
import time
from array import array
import math
import CMS_lumi
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gROOT.ProcessLine(".x tdrstyle.cc");

#ROOT.gSystem.Load("Util_cxx.so")
#from ROOT import draw_error_band

#python runFitPlots_vjets_signal_bigcombo_splitRes.py -n workspace_combo_"+signalName+".root  -l comboHPHP -i /afs/cern.ch/user/j/jngadiub/public/2016/JJ_nonRes_HPHP.root -M 1200 -s
#python runFitPlots_vjets_signal_bigcombo_splitRes.py -n workspace_combo_"+signalName+".root  -l comboHPLP -i /afs/cern.ch/user/j/jngadiub/public/2016/JJ_nonRes_HPLP.root -M 1200

addTT = False 
doFit = False
parser = optparse.OptionParser()
parser.add_option("-o","--output",dest="output",help="Output folder name",default='')
parser.add_option("-n","--name",dest="name",help="Input workspace",default='workspace.root')
parser.add_option("-i","--input",dest="input",help="Input nonRes histo",default='JJ_HPHP.root')
parser.add_option("-x","--xrange",dest="xrange",help="set range for x bins in projection",default="0,-1")
parser.add_option("-y","--yrange",dest="yrange",help="set range for y bins in projection",default="0,-1")
parser.add_option("-z","--zrange",dest="zrange",help="set range for z bins in projection",default="0,-1")
parser.add_option("-p","--projection",dest="projection",help="choose which projection should be done",default="xyz")
parser.add_option("-d","--data",dest="data",action="store_true",help="make also postfit plots",default=True)
parser.add_option("-l","--label",dest="label",help="add extra label such as pythia or herwig",default="")
parser.add_option("--log",dest="log",help="write fit result to log file",default="fit_results.log")
parser.add_option("--pdfz",dest="pdfz",help="name of pdfs lie PTZUp etc",default="")
parser.add_option("--pdfx",dest="pdfx",help="name of pdfs lie PTXUp etc",default="")
parser.add_option("--pdfy",dest="pdfy",help="name of pdfs lie PTYUp etc",default="")
parser.add_option("-s","--signal",dest="fitSignal",action="store_true",help="do S+B fit",default=False)
parser.add_option("-t","--addTop",dest="addTop",action="store_true",help="Fit top",default=False)
parser.add_option("-M","--mass",dest="signalMass",type=float,help="signal mass",default=1560.)
parser.add_option("--prelim",dest="prelim",type=int,help="add preliminary label",default=0)

(options,args) = parser.parse_args()
ROOT.gStyle.SetOptStat(0)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
#colors = [ROOT.kBlack,ROOT.kPink-1,ROOT.kAzure+1,ROOT.kAzure+1,210,210,ROOT.kMagenta,ROOT.kMagenta,ROOT.kOrange,ROOT.kOrange,ROOT.kViolet,ROOT.kViolet]
colors = [ROOT.kGray+2,ROOT.kRed,ROOT.kBlue,ROOT.kGray+1,210,210,ROOT.kMagenta,ROOT.kMagenta,ROOT.kOrange,ROOT.kOrange,ROOT.kViolet,ROOT.kViolet]

signalName = "ZprimeZH"
period = "2016"
if options.name.find("2017")!=-1: period = "2017"

def get_canvas(cname):

 #change the CMS_lumi variables (see CMS_lumi.py)
 CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
 CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
 CMS_lumi.writeExtraText = 1
 if options.prelim==0:
    CMS_lumi.extraText = " "
 else:
    CMS_lumi.extraText = "Preliminary"

 H_ref = 600 
 W_ref = 600 
 W = W_ref
 H  = H_ref

 iPeriod = 0

 # references for T, B, L, R
 T = 0.08*H_ref
 B = 0.12*H_ref 
 L = 0.12*W_ref
 R = 0.04*W_ref

 canvas = ROOT.TCanvas(cname,cname,50,50,W,H)
 canvas.SetFillColor(0)
 canvas.SetBorderMode(0)
 canvas.SetFrameFillStyle(0)
 canvas.SetFrameBorderMode(0)
 canvas.SetLeftMargin( L/W )
 canvas.SetRightMargin( R/W )
 canvas.SetTopMargin( T/H )
 canvas.SetBottomMargin( B/H )
 canvas.SetTickx()
 canvas.SetTicky()
 
 return canvas

def get_pad(name):

 #change the CMS_lumi variables (see CMS_lumi.py)
 CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
 CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
 CMS_lumi.lumi_13TeV = "77.3 fb^{-1}"
 CMS_lumi.writeExtraText = 1
 CMS_lumi.lumi_sqrtS = "13 TeV (2016+2017)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
 if options.prelim==0:
    CMS_lumi.extraText = " "
 if options.prelim==1:
    CMS_lumi.extraText = "Preliminary"
 if options.prelim==2:
     CMS_lumi.extraText = "Supplementary"

 CMS_lumi.lumiText = "(13 TeV)"
 iPos = 0
 if( iPos==0 ): CMS_lumi.relPosX = 0.014
 CMS_lumi.relPosX=0.05
 H_ref = 600 
 W_ref = 600 
 W = W_ref
 H  = H_ref

 iPeriod = 0
 iPeriod = 4
 # references for T, B, L, R
 T = 0.08*H_ref
 B = 0.12*H_ref 
 L = 0.12*W_ref
 R = 0.04*W_ref

 pad = ROOT.TPad(name, name, 0, 0.3, 1, 1.0)
 pad.SetFillColor(0)
 pad.SetBorderMode(0)
 pad.SetFrameFillStyle(0)
 pad.SetFrameBorderMode(0)
 #pad.SetLeftMargin( L/W )
 #pad.SetRightMargin( R/W )
 pad.SetTopMargin( T/H )
 #pad.SetBottomMargin( B/H )
 pad.SetTickx()
 pad.SetTicky()
 
 
 return pad


def getListFromRange(xyzrange):
    r=[]
    a,b = xyzrange.split(",")
    r.append(float(a))
    r.append(float(b))
    return r


def getListOfBins(hist,dim):
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
    bins ={}
    for i in range(1,N+1): bins[i] = axis.GetBinCenter(i) 
    
    return bins


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
    r=[]
    for i in range(1,N+2): r.append(axis.GetBinLowEdge(i)) 
    
    return array("d",r)


def getListOfBinsWidth(hist,dim):
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
    r ={}
    for i in range(0,N+2):
        #v = mmin + i * (mmax-mmin)/float(N)
        r[i] = axis.GetBinWidth(i) 
    return r 

    
def reduceBinsToRange(Bins,r):
    if r[0]==0 and r[1]==-1:
        return Bins
    result ={}
    for key, value in Bins.iteritems():
        if value >= r[0] and value <=r[1]:
            result[key]=value
    return result


def MakePlots(histos,hdata,hsig,axis,nBins,errors=None):
   
    extra1 = ''
    extra2 = ''
    htitle = ''
    xtitle = ''
    ymin = 0
    ymax = 0
    xrange = options.xrange
    yrange = options.yrange
    zrange = options.zrange
    if options.xrange == '0,-1': xrange = '55,215'
    if options.yrange == '0,-1': yrange = '55,215'
    if options.zrange == '0,-1': zrange = '1126,5500'
    if axis=='z':
     htitle = "Z-Proj. x : "+options.xrange+" y : "+options.yrange
     xtitle = "Dijet invariant mass [GeV]"
     ymin = 0.2
     ymax = hdata.GetMaximum()*10
     extra1 = xrange.split(',')[0]+' < m_{jet1} < '+ xrange.split(',')[1]+' GeV'
     extra2 = yrange.split(',')[0]+' < m_{jet2} < '+ yrange.split(',')[1]+' GeV'
    elif axis=='x':
     htitle = "X-Proj. y : "+options.yrange+" z : "+options.zrange
     xtitle = " m_{jet1} [GeV]"
     ymin = 0.02
     ymax = hdata.GetMaximum()*1.3
     extra1 = yrange.split(',')[0]+' < m_{jet2} < '+ yrange.split(',')[1]+' GeV'
     extra2 = zrange.split(',')[0]+' < m_{jj} < '+ zrange.split(',')[1]+' GeV'
    elif axis=='y':
     htitle = "Y-Proj. x : "+options.xrange+" z : "+options.zrange
     xtitle = " m_{jet2} [GeV]"
     ymin = 0.02
     ymax = hdata.GetMaximum()*1.3
     extra1 = xrange.split(',')[0]+' < m_{jet1} < '+ xrange.split(',')[1]+' GeV'
     extra2 = zrange.split(',')[0]+' < m_{jj} < '+ zrange.split(',')[1]+' GeV'
                   
    #leg = ROOT.TLegend(0.450436242,0.5531968,0.7231544,0.8553946)
    leg = ROOT.TLegend(0.50809045,0.5063636,0.7622613,0.8520979)
    leg.SetTextSize(0.05)
    c = ROOT.TCanvas('c')
    pad1 = get_pad("pad1") #ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    if axis == 'z': pad1.SetLogy()
    pad1.SetBottomMargin(0.01)    
    pad1.SetTopMargin(0.1) 
    pad1.Draw()
    pad1.cd()	 
 
    histos[0].SetMinimum(ymin)
    histos[0].SetMaximum(ymax) 
    histos[0].SetTitle(htitle)
    histos[0].SetLineColor(colors[0])
    histos[0].SetLineWidth(2)
    histos[0].GetXaxis().SetTitle(xtitle)
    histos[0].GetYaxis().SetTitleOffset(1.3)
    histos[0].GetYaxis().SetTitle("Events")
    histos[0].GetYaxis().SetTitleOffset(1.3)
    histos[0].GetYaxis().SetTitle("Events/ 2 GeV")
    if axis == 'z': histos[0].GetYaxis().SetTitle("Events/ 100 GeV")
    histos[0].GetYaxis().SetTitleSize(0.06)
    histos[0].GetYaxis().SetLabelSize(0.06)
    histos[0].Draw('HIST')

 
    histos[1].SetLineColor(colors[1])
    histos[1].SetLineWidth(2)
    
    histos[2].SetLineColor(colors[2])
    histos[2].SetLineWidth(2)
    
    
    if options.addTop:
      histos[3].SetLineColor(colors[3])
      histos[3].SetLineWidth(2)
      leg.AddEntry(histos[3],"t","l")
	   
    for i in range(4,len(histos)):
        histos[i].SetLineColor(colors[i])
        histos[i].Draw("histsame")
        name = histos[i].GetName().split("_")
        leg.AddEntry(histos[i],name[2],"l")

    hdata.SetMarkerStyle(20)
    hdata.SetMarkerColor(ROOT.kBlack)
    hdata.SetLineColor(ROOT.kBlack)
    hdata.SetMarkerSize(0.7)
    
    if errors!=None:
        errors[0].SetFillColor(colors[0])
        errors[0].SetFillStyle(3001)
        errors[0].SetLineColor(colors[0])
        errors[0].SetLineWidth(0)
        errors[0].SetMarkerSize(0)
    
    scaling = 1000.
    eff = 0.1
    if purity =="HPHP":
        scaling = 5
        eff = 0.02
    
    if hsig:
      if hsig.Integral()!=0.:   
        #hsig.Scale(1/hsig.Integral())
        #hsig.Scale(histos[2].Integral()*0.4)
        #hsig.Scale(scaling/hsig.Integral()*0.013146*77300.*eff)
        print hsig.Integral()
        hsig.Scale(scaling/0.0865311263651)#/0.270828488383)
        #hsig.Scale(scaling/1.12624917194)
        #hsig.Scale(scaling/7.8149570876e-09)
      hsig.SetFillColor(ROOT.kGreen-6)
      hsig.SetLineColor(ROOT.kBlack)
      hsig.SetLineStyle(1)
      hsig.Draw("HISTsame")
     #leg.AddEntry(hsig,"Signal pdf","F")
    
    #errors[0].Draw("E5same")
    if errors!=None:
        if axis=="z":
            errors[0].Draw("2same")
        else:
            errors[0].Draw("3same")
    histos[0].Draw("samehist")
    if options.addTop: histos[3].Draw("histsame") 
    histos[1].SetLineStyle(7)
    histos[2].SetLineStyle(6)
    histos[1].Draw("histsame") 
    histos[2].Draw("histsame")
    hdata.Draw("samePE0")         
    leg.SetLineColor(0)
    
    
    leg.AddEntry(hdata,"Data","ep")
    leg.AddEntry(histos[0],"Signal+background fit","l")
    if errors!=None:
        leg.AddEntry(errors[0],"#pm 1#sigma unc.","f")
    leg.AddEntry(histos[1],"W+jets, t#bar{t}","l")
    leg.AddEntry(histos[2],"Z+jets","l")
    
    text = "G_{bulk} (%.1f TeV) #rightarrow WW (#times %i)"%(options.signalMass/1000.,scaling)
    if (options.signalMass%1000.)==0:
       text = "G_{bulk} (%i TeV) #rightarrow WW (#times %i)"%(options.signalMass/1000.,scaling) 
    
    if signalName.find("ZprimeWW")!=-1:
        text = "Z' (%.1f TeV) #rightarrow WW (#times %i)"%(options.signalMass/1000.,scaling)
        if (options.signalMass%1000.)==0:
            text = "Z' (%i TeV) #rightarrow WW (#times %i)"%(options.signalMass/1000.,scaling) 
            
    if signalName.find("ZprimeZH")!=-1:
        text = "Z' (%.1f TeV) #rightarrow ZH (#times %i)"%(options.signalMass/1000.,scaling)
        if (options.signalMass%1000.)==0:
            text = "Z' (%i TeV) #rightarrow ZH (#times %i)"%(options.signalMass/1000.,scaling) 
     
    if options.fitSignal==True: 
        leg.AddEntry(hsig,text,"F")
    leg.Draw("same")
    
    #errors[0].Draw("E2same")
    
    chi2 = getChi2proj(histos[0],hdata)
    print hdata.GetEntries(),hdata.Integral()
    print "Projection %s: Chi2/ndf = %.2f/%i"%(axis,chi2[0],chi2[1]),"= %.2f"%(chi2[0]/chi2[1])," prob = ",ROOT.TMath.Prob(chi2[0],chi2[1])

    pt = ROOT.TPaveText(0.18,0.06,0.54,0.17,"NDC")
    pt.SetTextFont(42)
    pt.SetTextSize(0.05)
    pt.SetTextAlign(12)
    pt.SetFillColor(0)
    pt.SetBorderSize(0)
    pt.SetFillStyle(0)
    pt.AddText("Chi2/ndf = %.2f/%i = %.2f"%(chi2[0],chi2[1],chi2[0]/chi2[1]))
    pt.AddText("Prob = %.3f"%ROOT.TMath.Prob(chi2[0],chi2[1]))
    #pt.Draw()

    pt2 = ROOT.TPaveText(0.55,0.29,0.99,0.4,"NDC")
    pt2.SetTextFont(42)
    pt2.SetTextSize(0.05)
    pt2.SetTextAlign(12)
    pt2.SetFillColor(0)
    pt2.SetBorderSize(0)
    pt2.SetFillStyle(0)
    pt2.AddText(extra1)
    pt2.AddText(extra2)
    pt2.Draw()

    pt3 = ROOT.TPaveText(0.65,0.39,0.99,0.52,"NDC")
    pt3.SetTextFont(42)
    pt3.SetTextSize(0.05)
    pt3.SetTextAlign(12)
    pt3.SetFillColor(0)
    pt3.SetBorderSize(0)
    pt3.SetFillStyle(0)
    pt3.AddText("%s category"%purity)
    #pt3.AddText(extra1)
    #pt3.AddText(extra2)
    pt3.Draw()

    CMS_lumi.CMS_lumi(pad1, 4, 10)
        
    pad1.Modified()
    pad1.Update()
    
    c.Update()
    c.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0.01)
    pad2.SetBottomMargin(0.4)
    #pad2.SetGridy()
    pad2.Draw()
    pad2.cd()
    
    #for ratio
    #graphs = addRatioPlot(hdata,histos[0],nBins,errors[0])
    #graphs[1].Draw("AP")
    #graphs[0].Draw("E3same")
    #graphs[1].Draw("Psame")
    
    #for pulls
    if errors ==None: errors=[0,0];
    if options.name.find('sigonly')!=-1: graphs = addPullPlot(hdata,hsig,nBins,errors[0])
    else:
        graphs = addPullPlot(hdata,histos[0],nBins,errors[0])
    # graphs = addRatioPlot(hdata,histos[0],nBins,errors[0])
    graphs[0].Draw("HIST")

    pad2.Modified()
    pad2.Update()
    
    c.cd()
    c.Update()
    c.Modified()
    c.Update()
    c.cd()
    c.SetSelected(c)
    #errors[0].Draw("E2same")
    #CMS_lumi.CMS_lumi(c, 0, 11)
    #c.cd()
    #c.Update()
    #c.RedrawAxis()
    #frame = c.GetFrame()
    #frame.Draw()
    if options.prelim==0:
        c.SaveAs(options.output+"PostFit"+options.label+"_"+htitle.replace(' ','_').replace('.','_').replace(':','_').replace(',','_')+".png")
        c.SaveAs(options.output+"PostFit"+options.label+"_"+htitle.replace(' ','_').replace('.','_').replace(':','_').replace(',','_')+".pdf")
        c.SaveAs(options.output+"PostFit"+options.label+"_"+htitle.replace(' ','_').replace('.','_').replace(':','_').replace(',','_')+".C")
        c.SaveAs(options.output+"PostFit"+options.label+"_"+htitle.replace(' ','_').replace('.','_').replace(':','_').replace(',','_')+".root")
    else:
        c.SaveAs(options.output+"PostFit"+options.label+"_"+htitle.replace(' ','_').replace('.','_').replace(':','_').replace(',','_')+"_prelim.png")
        c.SaveAs(options.output+"PostFit"+options.label+"_"+htitle.replace(' ','_').replace('.','_').replace(':','_').replace(',','_')+"_prelim.pdf")
        c.SaveAs(options.output+"PostFit"+options.label+"_"+htitle.replace(' ','_').replace('.','_').replace(':','_').replace(',','_')+"_prelim.C")
        c.SaveAs(options.output+"PostFit"+options.label+"_"+htitle.replace(' ','_').replace('.','_').replace(':','_').replace(',','_')+"_prelim.root")
    

def doZprojection(pdfs,data,norm_nonres,norm_Wres,pdf_sig,norm_sig,norm_Zres,norm_TThad): 
    data = ROOT.RooDataHist("data","data",args_ws,data)
    # do some z projections
    h=[]
    lv=[]
    lv1_sig=[]
    h1_sig = 0    
    dh = ROOT.TH1F("dh","dh",len(zBinslowedge)-1,zBinslowedge)
    neventsPerBin_1 = {}
    for zk,zv in zBins_redux.iteritems():
        neventsPerBin_1[zk]=0 
   
    for p in pdfs:
        h.append( ROOT.TH1F("h_"+p.GetName(),"h_"+p.GetName(),len(zBinslowedge)-1,zBinslowedge))
        lv.append({})
    for i in range(0,len(pdfs)):
        for zk,zv in zBins_redux.iteritems(): lv[i][zv]=0  
	
    if pdf_sig:
        h1_sig = ROOT.TH1F("h1_"+pdf_sig.GetName(),"h1_"+pdf_sig.GetName(),len(zBinslowedge)-1,zBinslowedge)
        lv1_sig.append({})
	for zk,zv in zBins_redux.iteritems(): lv1_sig[0][zv]=0
    	
    for xk, xv in xBins_redux.iteritems():
         MJ1.setVal(xv)
         for yk, yv in yBins_redux.iteritems():
             MJ2.setVal(yv)
             for zk,zv in zBins_redux.iteritems():
                 MJJ.setVal(zv)
                 #dh.Fill(zv,data.weight(argset))
		 neventsPerBin_1[zk] += data.weight(argset)
                 i=0
                 binV = zBinsWidth[zk]*xBinsWidth[xk]*yBinsWidth[yk]
                 for p in pdfs:
                    lv[i][zv] += p.getVal(argset)*binV
                    i+=1
		 if pdf_sig:
                     for s in range(0,100):
                        MJJ.setVal(zv-zBinsWidth[zk]//2.+s*zBinsWidth[zk]/100.)
                        lv1_sig[0][zv] += pdf_sig.getVal(argset)*binV/100.

		 		    
    for i in range(0,len(pdfs)):
        for zk,zv in zBins_redux.iteritems():
           if "nonRes" in str(pdfs[i].GetName()) : h[i].Fill(zv,lv[i][zv]*norm_nonres[0])
	   elif "Wjet" in str(pdfs[i].GetName()) : h[i].Fill(zv,lv[i][zv]*norm_Wres[0])
	   elif "Zjet" in str(pdfs[i].GetName()) : h[i].Fill(zv,lv[i][zv]*norm_Zres[0])
	   elif "TThad" in str(pdfs[i].GetName()): h[i].Fill(zv,lv[i][zv]*norm_TThad[0])
	   else: h[i].Fill(zv,lv[i][zv]*(norm_Wres[0]+norm_Zres[0]+norm_TThad[0]+norm_nonres[0]))
    
    if pdf_sig:
        for zk,zv in zBins_redux.iteritems(): h1_sig.Fill(zv,lv1_sig[0][zv]*norm_sig[0])    
     

    htot_sig = ROOT.TH1F("htot_sig","htot_sig",len(zBinslowedge)-1,zBinslowedge)
    htot_sig.Add(h1_sig)

    htot_nonres = ROOT.TH1F("htot_nonres","htot_nonres",len(zBinslowedge)-1,zBinslowedge)
    htot_nonres.Add(h[0])

    htot_Wres = ROOT.TH1F("htot_Wres","htot_Wres",len(zBinslowedge)-1,zBinslowedge)
    htot_Wres.Add(h[1])
    
    htot_Zres = ROOT.TH1F("htot_Zres","htot_Zres",len(zBinslowedge)-1,zBinslowedge)
    htot_Zres.Add(h[2])
    
    htot_TThad = ROOT.TH1F("htot_TThad","htot_TThad",len(zBinslowedge)-1,zBinslowedge)
    htot_TThad.Add(h[3])
            	    
    htot = ROOT.TH1F("htot","htot",len(zBinslowedge)-1,zBinslowedge)
    htot.Add(htot_nonres)
    htot.Add(htot_Wres)
    htot.Add(htot_Zres)
    if options.addTop: htot.Add(htot_TThad)
    htot.Add(htot_sig)

    hfinals = []
    hfinals.append(htot)
    hfinals.append(htot_Wres)
    hfinals.append(htot_Zres)
    if options.addTop: hfinals.append(htot_TThad)
    for i in range(10,len(h)): hfinals.append(h[i])    
    for b,v in neventsPerBin_1.iteritems(): dh.SetBinContent(b,neventsPerBin_1[b])
    dh.SetBinErrorOption(ROOT.TH1.kPoisson)
    if doFit:
        errors = draw_error_band(htot,norm_nonres[0]+norm_Wres[0]+norm_Zres[0]+norm_TThad[0],math.sqrt(norm_nonres[1]*norm_nonres[1]+norm_Wres[1]*norm_Wres[1]+norm_Zres[1]*norm_Zres[1]+norm_TThad[1]*norm_TThad[1]),pdfs[-1],zBinslowedge,'z')
        MakePlots(hfinals,dh,htot_sig,'z',zBinslowedge,errors)
    else: MakePlots(hfinals,dh,htot_sig,'z',zBinslowedge)
        
def doXprojection(pdfs,data1,norm1_nonres,norm1_Wres,pdf1_sig,norm1_sig,norm1_Zres,norm1_TThad):
    data1 = ROOT.RooDataHist("data1","data1",args_ws,data1)
    # do some x projections
    h=[]
    lv=[]
    h1_sig=0
    lv1_sig=[]
    h2_sig = 0    
    lv2_sig=[]
    dh = ROOT.TH1F("dh","dh",len(xBinslowedge)-1,xBinslowedge)
    neventsPerBin_1 = {}
    for xk,xv in xBins_redux.iteritems():
        neventsPerBin_1[xk]=0 
    
    for p in pdfs:
        h.append( ROOT.TH1F("h_"+p.GetName(),"h_"+p.GetName(),len(xBinslowedge)-1,xBinslowedge))
        lv.append({})
    for i in range(0,len(pdfs)):
        for xk,xv in xBins_redux.iteritems(): lv[i][xv]=0   
	
    if pdf1_sig:
        h1_sig = ROOT.TH1F("h1_"+pdf1_sig.GetName(),"h1_"+pdf1_sig.GetName(),len(xBinslowedge)-1,xBinslowedge)
        lv1_sig.append({})
        for xk,xv in xBins_redux.iteritems(): lv1_sig[0][xv]=0 

    for xk, xv in xBins_redux.iteritems():
         MJ1.setVal(xv)
         for yk, yv in yBins_redux.iteritems():
             MJ2.setVal(yv)
             for zk,zv in zBins_redux.iteritems():
                 MJJ.setVal(zv)
                 #dh.Fill(zv,data.weight(argset))
		 neventsPerBin_1[xk] += data1.weight(argset)
                 i=0
                 binV = zBinsWidth[zk]*xBinsWidth[xk]*yBinsWidth[yk]
                 for p in pdfs:
                    lv[i][xv] += p.getVal(argset)*binV
                    i+=1
		 if pdf1_sig: lv1_sig[0][xv] += pdf1_sig.getVal(argset)*binV 
            
		    
    for i in range(0,len(pdfs)):
        for xk,xv in xBins_redux.iteritems():
           if "nonRes" in str(pdfs[i].GetName()) : h[i].Fill(xv,lv[i][xv]*norm1_nonres[0])
	   elif "Wjet" in str(pdfs[i].GetName()) : h[i].Fill(xv,lv[i][xv]*norm1_Wres[0])
	   elif "Zjet" in str(pdfs[i].GetName()) : h[i].Fill(xv,lv[i][xv]*norm1_Zres[0])
	   elif "TThad" in str(pdfs[i].GetName()): h[i].Fill(xv,lv[i][xv]*norm1_TThad[0])
	   else: h[i].Fill(xv,lv[i][xv]*(norm1_Wres[0]+norm1_Zres[0]+norm1_TThad[0]+norm1_nonres[0]))
	       
    if pdf1_sig:
        for xk,xv in xBins_redux.iteritems(): h1_sig.Fill(xv,lv1_sig[0][xv]*norm1_sig[0])    
    
    
	   
    htot_sig = ROOT.TH1F("htot_sig","htot_sig",len(xBinslowedge)-1,xBinslowedge)
    htot_sig.Add(h1_sig)
   
    htot_nonres = ROOT.TH1F("htot_nonres","htot_nonres",len(xBinslowedge)-1,xBinslowedge)
    htot_nonres.Add(h[0])
 
    htot_Wres = ROOT.TH1F("htot_Wres","htot_Wres",len(xBinslowedge)-1,xBinslowedge)
    htot_Wres.Add(h[1])
   
    htot_Zres = ROOT.TH1F("htot_Zres","htot_Zres",len(xBinslowedge)-1,xBinslowedge)
    htot_Zres.Add(h[2])
    
    htot_TThad = ROOT.TH1F("htot_TThad","htot_TThad",len(xBinslowedge)-1,xBinslowedge)
    htot_TThad.Add(h[3])
            	    
    htot = ROOT.TH1F("htot","htot",len(xBinslowedge)-1,xBinslowedge)
    htot.Add(htot_nonres)
    htot.Add(htot_Wres)
    htot.Add(htot_Zres)
    if options.addTop: htot.Add(htot_TThad)
    htot.Add(htot_sig)

    hfinals = []
    hfinals.append(htot)
    hfinals.append(htot_Wres)
    hfinals.append(htot_Zres)
    if options.addTop: hfinals.append(htot_TThad)
    for i in range(10,len(h)): hfinals.append(h[i])    
    for b,v in neventsPerBin_1.iteritems(): dh.SetBinContent(b,neventsPerBin_1[b])
    dh.SetBinErrorOption(ROOT.TH1.kPoisson)
    if doFit:
        errors = draw_error_band(htot,norm1_nonres[0]+norm1_Wres[0]+norm1_Zres[0]+norm1_TThad[0],math.sqrt(norm1_nonres[1]*norm1_nonres[1]+norm1_Wres[1]*norm1_Wres[1]+norm1_Zres[1]*norm1_Zres[1]+norm1_TThad[1]*norm1_TThad[1]),pdfs[-1],xBinslowedge,'x')
        MakePlots(hfinals,dh,htot_sig,'x',xBinslowedge,errors)
    else: MakePlots(hfinals,dh,htot_sig,'x',xBinslowedge)

def doYprojection(pdfs,data1,norm1_nonres,norm1_Wres,pdf1_sig,norm1_sig,norm1_Zres,norm1_TThad):
    data1 = ROOT.RooDataHist("data1","data1",args_ws,data1)
    # do some y projections
    h=[]
    lv=[]
    lv1_sig=[]
    h1_sig = 0    
    lv2_sig=[]
    h2_sig = 0   
    dh = ROOT.TH1F("dh","dh",len(yBinslowedge)-1,yBinslowedge)
    neventsPerBin_1 = {}
    for yk,yv in yBins_redux.iteritems():
        neventsPerBin_1[yk]=0 
  
    for p in pdfs:
        h.append( ROOT.TH1F("h_"+p.GetName(),"h_"+p.GetName(),len(yBinslowedge)-1,yBinslowedge))
        lv.append({})
    for i in range(0,len(pdfs)):
        for yk,yv in yBins_redux.iteritems(): lv[i][yv]=0  
	
    if pdf1_sig:
        h1_sig = ROOT.TH1F("h1_"+pdf1_sig.GetName(),"h1_"+pdf1_sig.GetName(),len(yBinslowedge)-1,yBinslowedge)
        lv1_sig.append({})
        for yk,yv in yBins_redux.iteritems(): lv1_sig[0][yv]=0;
    for xk, xv in xBins_redux.iteritems():
         MJ1.setVal(xv)
         for yk, yv in yBins_redux.iteritems():
             MJ2.setVal(yv)
             for zk,zv in zBins_redux.iteritems():
                 MJJ.setVal(zv)
                 #dh.Fill(zv,data.weight(argset))
		 neventsPerBin_1[yk] += data1.weight(argset)
                 i=0
                 binV = zBinsWidth[zk]*xBinsWidth[xk]*yBinsWidth[yk]
                 for p in pdfs:
                    lv[i][yv] += p.getVal(argset)*binV
                    i+=1
		 if pdf1_sig: lv1_sig[0][yv] += pdf1_sig.getVal(argset)*binV 
    for i in range(0,len(pdfs)):
        for yk,yv in yBins_redux.iteritems():
           if "nonRes" in str(pdfs[i].GetName()) : h[i].Fill(yv,lv[i][yv]*norm1_nonres[0])
	   elif "Wjet" in str(pdfs[i].GetName()) : h[i].Fill(yv,lv[i][yv]*norm1_Wres[0])
	  
	   elif "Zjet" in str(pdfs[i].GetName()) : h[i].Fill(yv,lv[i][yv]*norm1_Zres[0])
	   
	   elif "TThad" in str(pdfs[i].GetName()): h[i].Fill(yv,lv[i][yv]*norm1_TThad[0])
	   
	   else: h[i].Fill(yv,lv[i][yv]*(norm1_Wres[0]+norm1_Zres[0]+norm1_TThad[0]+norm1_nonres[0]))
	       
    if pdf1_sig:
        for yk,yv in yBins_redux.iteritems(): h1_sig.Fill(yv,lv1_sig[0][yv]*norm1_sig[0])    
   
	    
    htot_sig = ROOT.TH1F("htot_sig","htot_sig",len(yBinslowedge)-1,yBinslowedge)
    htot_sig.Add(h1_sig)

    htot_nonres = ROOT.TH1F("htot_nonres","htot_nonres",len(yBinslowedge)-1,yBinslowedge)
    htot_nonres.Add(h[0])
    
    htot_Wres = ROOT.TH1F("htot_Wres","htot_Wres",len(yBinslowedge)-1,yBinslowedge)
    htot_Wres.Add(h[1])
   
    htot_Zres = ROOT.TH1F("htot_Wres","htot_Wres",len(yBinslowedge)-1,yBinslowedge)
    htot_Zres.Add(h[2])
    
    htot_TThad = ROOT.TH1F("htot_TThad","htot_TThad",len(yBinslowedge)-1,yBinslowedge)
    htot_TThad.Add(h[3])
             	    
    htot = ROOT.TH1F("htot","htot",len(yBinslowedge)-1,yBinslowedge)
    htot.Add(htot_nonres)
    htot.Add(htot_Wres)
    htot.Add(htot_Zres)
    if options.addTop:htot.Add(htot_TThad)
    htot.Add(htot_sig)

    hfinals = []
    hfinals.append(htot)
    hfinals.append(htot_Wres)
    hfinals.append(htot_Zres)
    if options.addTop:hfinals.append(htot_TThad)
    for i in range(10,len(h)): hfinals.append(h[i])    
    for b,v in neventsPerBin_1.iteritems(): dh.SetBinContent(b,neventsPerBin_1[b])
    dh.SetBinErrorOption(ROOT.TH1.kPoisson)
    if doFit:
        errors = draw_error_band(htot,norm1_nonres[0]+norm1_Wres[0]+norm1_Zres[0]+norm1_TThad[0],math.sqrt(norm1_nonres[1]*norm1_nonres[1]+norm1_Wres[1]*norm1_Wres[1]+norm1_Zres[1]*norm1_Zres[1]+norm1_TThad[1]*norm1_TThad[1]),pdfs[-1],yBinslowedge,'y')
        MakePlots(hfinals,dh,htot_sig,'y',yBinslowedge,errors)
    else: MakePlots(hfinals,dh,htot_sig,'y',yBinslowedge)
 

def addPullPlot(hdata,hpostfit,nBins,error_band):
    #print "make pull plots: (data-fit)/sigma_data"
    N = hdata.GetNbinsX()
    gpost = ROOT.TGraphErrors(0)
    gt = ROOT.TH1F("gt","gt",len(nBins)-1,nBins)
    for i in range(1,N+1):
        m = hdata.GetXaxis().GetBinCenter(i)
        #ypostfit = (hdata.GetBinContent(i) - hpostfit.GetBinContent(i))/hdata.GetBinErrorUp(i)
        if hpostfit.GetBinContent(i) <= hdata.GetBinContent(i):
            if error_band !=0: ypostfit = (hdata.GetBinContent(i) - hpostfit.GetBinContent(i))/ ROOT.TMath.Abs(hdata.GetBinErrorUp(i))
            else: ypostfit = (hdata.GetBinContent(i) - hpostfit.GetBinContent(i))/ hdata.GetBinErrorUp(i)
        else:
            if error_band!=0: ypostfit = (hdata.GetBinContent(i) - hpostfit.GetBinContent(i))/ ROOT.TMath.Sqrt(ROOT.TMath.Abs( pow(hdata.GetBinErrorUp(i),2) - pow(error_band.GetErrorYlow(i-1),2) ))
            else: ypostfit = (hdata.GetBinContent(i) - hpostfit.GetBinContent(i))/ hdata.GetBinErrorUp(i)
        gpost.SetPoint(i-1,m,ypostfit)
        gt.SetBinContent(i,ypostfit)
	#print "bin",i,"x",m,"data",hdata.GetBinContent(i),"post fit",hpostfit.GetBinContent(i),"err data",hdata.GetBinErrorUp(i),"err fit",error_band.GetBinError(i),"pull postfit",ypostfit
	#print "bin",i,"x",m,"data",hdata.GetBinContent(i),"post fit",hpostfit.GetBinContent(i),"err data",hdata.GetBinErrorUp(i),"err fit",error_band.GetErrorYhigh(i-1),"pull postfit",ypostfit
		
    gpost.SetLineColor(colors[1])
    gpost.SetMarkerColor(colors[1])
    gpost.SetFillColor(ROOT.kGray+3)
    gpost.SetMarkerSize(1)
    gpost.SetMarkerStyle(20)
    gt.SetFillColor(ROOT.kGray+3)
    gt.SetLineColor(ROOT.kGray+3)
    
    #gt = ROOT.TH1F("gt","gt",hdata.GetNbinsX(),hdata.GetXaxis().GetXmin(),hdata.GetXaxis().GetXmax())
    #gt = ROOT.TH1F("gt","gt",len(nBins)-1,nBins)
    gt.SetTitle("")
    #gt.SetMinimum(0.5);
    #gt.SetMaximum(1.5);
    gt.SetMinimum(-2.499);
    gt.SetMaximum(2.499);
    gt.SetDirectory(0);
    gt.SetStats(0);
    gt.SetLineStyle(0);
    gt.SetMarkerStyle(20);
    gt.GetXaxis().SetTitle(hpostfit.GetXaxis().GetTitle());
    gt.GetXaxis().SetLabelFont(42);
    gt.GetXaxis().SetLabelOffset(0.02);
    gt.GetXaxis().SetLabelSize(0.17);
    gt.GetXaxis().SetTitleSize(0.15);
    gt.GetXaxis().SetTitleOffset(1.2);
    gt.GetXaxis().SetTitleFont(42);
    gt.GetYaxis().SetTitle("#frac{Data-fit}{#sigma}");
    gt.GetYaxis().CenterTitle(True);
    gt.GetYaxis().SetNdivisions(205);
    gt.GetYaxis().SetLabelFont(42);
    gt.GetYaxis().SetLabelOffset(0.007);
    gt.GetYaxis().SetLabelSize(0.15);
    gt.GetYaxis().SetTitleSize(0.15);
    gt.GetYaxis().SetTitleOffset(0.4);
    gt.GetYaxis().SetTitleFont(42);
    gt.GetXaxis().SetNdivisions(505)
    #gpre.SetHistogram(gt);
    #gpost.SetHistogram(gt);       
    return [gt] 
    
def addRatioPlot(hdata,hpostfit,nBins,error_band):
    #print "make pull plots: (data-fit)/sigma_data"
      N = hdata.GetNbinsX()
      gpost = ROOT.TGraphErrors(0)
      gt = ROOT.TH1F("gt","gt",len(nBins)-1,nBins)
      for i in range(1,N+1):
          m = hdata.GetXaxis().GetBinCenter(i)
          ypostfit = (hdata.GetBinContent(i)/hpostfit.GetBinContent(i))
          gpost.SetPoint(i-1,m,ypostfit)
          gt.SetBinContent(i,ypostfit)
          print "bin",i,"x",m,"data",hdata.GetBinContent(i),"post fit",hpostfit.GetBinContent(i),"err data",hdata.GetBinErrorUp(i),"err fit",error_band.GetBinError(i),"pull postfit",ypostfit
		
      gpost.SetLineColor(colors[1])
      gpost.SetMarkerColor(colors[1])
      gpost.SetFillColor(ROOT.kBlue)
      gpost.SetMarkerSize(1)
      gpost.SetMarkerStyle(20)
      gt.SetFillColor(ROOT.kBlue)
      gt.SetLineColor(ROOT.kBlue)
    
      #gt = ROOT.TH1F("gt","gt",hdata.GetNbinsX(),hdata.GetXaxis().GetXmin(),hdata.GetXaxis().GetXmax())
      #gt = ROOT.TH1F("gt","gt",len(nBins)-1,nBins)
      gt.SetTitle("")
      #gt.SetMinimum(0.5);
      #gt.SetMaximum(1.5);
      gt.SetMinimum(0.001);
      gt.SetMaximum(1.999);
      gt.SetDirectory(0);
      gt.SetStats(0);
      gt.SetLineStyle(0);
      gt.SetMarkerStyle(20);
      gt.GetXaxis().SetTitle(hpostfit.GetXaxis().GetTitle());
      gt.GetXaxis().SetLabelFont(42);
      gt.GetXaxis().SetLabelOffset(0.02);
      gt.GetXaxis().SetLabelSize(0.17);
      gt.GetXaxis().SetTitleSize(0.15);
      gt.GetXaxis().SetTitleOffset(1);
      gt.GetXaxis().SetTitleFont(42);
      gt.GetYaxis().SetTitle("#frac{data}{fit}");
      gt.GetYaxis().CenterTitle(True);
      gt.GetYaxis().SetNdivisions(205);
      gt.GetYaxis().SetLabelFont(42);
      gt.GetYaxis().SetLabelOffset(0.007);
      gt.GetYaxis().SetLabelSize(0.15);
      gt.GetYaxis().SetTitleSize(0.15);
      gt.GetYaxis().SetTitleOffset(0.4);
      gt.GetYaxis().SetTitleFont(42);
      gt.GetXaxis().SetNdivisions(505)
      #gpre.SetHistogram(gt);
      #gpost.SetHistogram(gt);       
      return [gt]


def draw_error_band(histo_central,norm1,err_norm1,rpdf1,x_min,proj):
    
    rand = ROOT.TRandom3(1234);
    number_errorband = 5
    syst = [0 for i in range(number_errorband)]
      
    value = [0 for x in range(len(x_min))]  
    number_point = len(value)
    
    par_pdf1 = rpdf1.getParameters(argset)  
    
    
    iter = par_pdf1.createIterator()
    var = iter.Next()
    print var.GetName()
    
   # par_pdf1.Print()
   # par_pdf2.Print()
      
    for j in range(number_errorband):
    
     #print j
     syst[j] = ROOT.TGraph(number_point+1);
     
     #paramters value are randomized using rfres and this can be done also if they are not decorrelate
     par_tmp = ROOT.RooArgList(fitresult_bkg_only.randomizePars())
     iter = par_pdf1.createIterator()
     var = iter.Next()

     while var:
      index = par_tmp.index(var.GetName())
      if index != -1:
       #print "pdf1",var.GetName(), var.getVal()
       var.setVal(par_tmp.at(index).getVal())     
       #print " ---> new value: ",var.getVal()
      var = iter.Next()
      
           
     norm1_tmp = rand.Gaus(norm1,err_norm1); #new poisson random number of events
     value = [0 for i in range(number_point)]
     for xk, xv in xBins_redux.iteritems():
       MJ1.setVal(xv)
       for yk, yv in yBins_redux.iteritems():
        MJ2.setVal(yv)
	for zk,zv in zBins_redux.iteritems():
         MJJ.setVal(zv)
	 binV = zBinsWidth[zk]*xBinsWidth[xk]*yBinsWidth[yk]
	 if proj == 'z': value[zk-1] += (norm1_tmp*rpdf1.getVal(argset)*binV)
	 elif proj == 'y': value[yk-1] += (norm1_tmp*rpdf1.getVal(argset)*binV)
	 elif proj == 'x': value[xk-1] += (norm1_tmp*rpdf1.getVal(argset)*binV)

     for ix,x in enumerate(x_min): 
      syst[j].SetPoint(ix, x, value[ix])

    #Try to build and find max and minimum for each point --> not the curve but the value to do a real envelope -> take one 2sigma interval        
    errorband = ROOT.TGraphAsymmErrors()#ROOT.TH1F("errorband","errorband",len(x_min)-1,x_min)

    val = [0 for i in range(number_errorband)]
    for ix,x in enumerate(x_min):
    
     for j in range(number_errorband):
      val[j]=(syst[j]).GetY()[ix]
     val.sort()
     errorband.SetPoint(ix,x_min[ix]+histo_central.GetBinWidth(ix+1)/2.,histo_central.GetBinContent(ix+1))
     #errorband.SetBinContent(ix+1,histo_central.GetBinContent(ix+1))
     #print "set bin content error band "+str(histo_central.GetBinContent(ix+1))+" for bin "+str(ix+1)
     errup = (val[int(0.84*number_errorband)]-histo_central.GetBinContent(ix+1)) #ROOT.TMath.Abs
     errdn = ( histo_central.GetBinContent(ix+1)-val[int(0.16*number_errorband)])
     print "error up "+str(errup)+" error down "+str(errdn)
     #if errup > errdn: errorband.SetBinError(ix+1,errup)
     #else: errorband.SetBinError(ix+1,errdn)
     #print ix,ix+1,histo_central.GetBinContent(ix+1),errorband.GetBinContent(ix+1)
     errorband.SetPointError(ix,histo_central.GetBinWidth(ix+1)/2.,histo_central.GetBinWidth(ix+1)/2.,ROOT.TMath.Abs(errdn),ROOT.TMath.Abs(errup))
    errorband.SetFillColor(ROOT.kBlack)
    errorband.SetFillStyle(3008)
    errorband.SetLineColor(ROOT.kGreen)
    errorband.SetMarkerSize(0)
       
    return [errorband]
	
	
def builtFittedPdf(pdfs,coefficients):
    result = RooAddPdf(pdfs,coefficients)
    return result


def getChi2fullModel(pdf,data,norm):
    pr=[]
    dr=[]
    for xk, xv in xBins.iteritems():
         MJ1.setVal(xv)
         for yk, yv in yBins.iteritems():
             MJ2.setVal(yv)
             for zk,zv in zBins.iteritems():
                 MJJ.setVal(zv)
                 dr.append(data.weight(argset))
                 binV = zBinsWidth[zk]*xBinsWidth[xk]*yBinsWidth[yk]
                 pr.append( pdf.getVal(argset)*binV*norm)
    ndof = 0
    chi2 = 0
    for i in range(0,len(pr)):
        if dr[i] < 10e-10:
            continue
        ndof+=1
        #chi2+= pow((dr[i] - pr[i]),2)/pr[i]
        chi2+= 2*( pr[i] - dr[i] + dr[i]* ROOT.TMath.Log(dr[i]/pr[i]))

    return [chi2,ndof]

def getChi2proj(histo_pdf,histo_data):
    pr=[]
    dr=[]
    for b in range(1,histo_pdf.GetNbinsX()+1):
     dr.append(histo_data.GetBinContent(b))
     pr.append(histo_pdf.GetBinContent(b))
    
    ndof = 0
    chi2 = 0
    for i in range(0,len(pr)):
        if dr[i] < 10e-10:
            continue
        ndof+=1
        #chi2+= pow((dr[i] - pr[i]),2)/pr[i]
	#print i,dr[i],pr[i],(dr[i] - pr[i]),pow((dr[i] - pr[i]),2)/pr[i],(dr[i] - pr[i])/histo_data.GetBinError(i+1)
        chi2+= 2*( pr[i] - dr[i] + dr[i]* ROOT.TMath.Log(dr[i]/pr[i]))

    return [chi2,ndof]   

if __name__=="__main__":
     finMC = ROOT.TFile(options.input,"READ");
     hinMC = finMC.Get("nonRes");
     purity = options.input.replace('.root','').split('_')[-1]   
                        
     #################################################
     xBins= getListOfBins(hinMC,"x")
     xBinslowedge = getListOfBinsLowEdge(hinMC,'x')
     xBinsWidth   = getListOfBinsWidth(hinMC,"x")
     print "x bins:"
     print xBins
     print "x bins low edge:"
     print xBinslowedge
     print "x bins width:"
     print xBinsWidth
     
     #################################################
     print
     yBins= getListOfBins(hinMC,"y")
     yBinslowedge = getListOfBinsLowEdge(hinMC,'y')     
     yBinsWidth   = getListOfBinsWidth(hinMC,"y")
     print "y bins:"
     print yBins
     print "y bins low edge:"
     print yBinslowedge
     print "y bins width:"
     print yBinsWidth
     
     #################################################
     print 
     zBins= getListOfBins(hinMC,"z")
     zBinslowedge = getListOfBinsLowEdge(hinMC,'z')
     zBinsWidth   = getListOfBinsWidth(hinMC,"z")
     print "z bins:"
     print zBins
     print "z bins low edge:"
     print zBinslowedge
     print "z bins width:"
     print zBinsWidth
            
     #################################################                
     print 
     print "open file " +options.name
     f = ROOT.TFile(options.name,"READ")
     workspace = f.Get("w")
     f.Close()
     #workspace.Print()

     model = workspace.pdf("model_b") 
     model_b = workspace.pdf("model_b")
     if options.fitSignal: model = workspace.pdf("model_s")
     data_all = workspace.data("data_obs")
     data_all.Print()
     data1 = workspace.data("data_obs").reduce("CMS_channel==CMS_channel::JJ_"+purity+"_13TeV_"+period)
     data1.Print()
     print
     print "Observed number of events in",purity,"category:",data1.sumEntries(),"("+period+")"
     print model
     args  = model.getComponents()
     pdf1Name = "pdf_binJJ_"+purity+"_13TeV_"+period+"_bonly"
     print pdf1Name
     if options.fitSignal:
      pdf1Name = "pdf_binJJ_"+purity+"_13TeV_"+period
     print args
     print "Expected number of QCD events:",(args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_nonRes"].getVal(),"("+period+")"
     print "Expected number of W+jets events:",(args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_Wjets"].getVal(),"("+period+")"
     print "Expected number of Z+jets events:",(args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_Zjets"].getVal(),"("+period+")"
     
     qcd_exp        = (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_nonRes"].getVal()
     wjets_exp      = (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_Wjets"].getVal()
     zjets_exp      = (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_Zjets"].getVal()
         
     
     
     if options.addTop:
       print "Expected number of tt events:",(args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_TThad"].getVal(),"("+period+")"
     if options.fitSignal:
      workspace.var("MH").setVal(options.signalMass)
      workspace.var("MH").setConstant(1)
      #workspace.var("r").setRange(0,1000)
      #workspace.var("r").setVal(9)
      print "Expected signal yields:",(args[pdf1Name].getComponents())["n_exp_final_binJJ_"+purity+"_13TeV_"+period+"_proc_"+signalName+""].getVal(),"("+period+")"
     print 

     #################################################
     print
     if doFit:
        fitresult = model.fitTo(data_all,ROOT.RooFit.SumW2Error(not(options.data)),ROOT.RooFit.Minos(0),ROOT.RooFit.Verbose(0),ROOT.RooFit.Save(1),ROOT.RooFit.NumCPU(8))  
        
        fitresult_bkg_only = model_b.fitTo(data_all,ROOT.RooFit.SumW2Error(not(options.data)),ROOT.RooFit.Minos(0),ROOT.RooFit.Verbose(0),ROOT.RooFit.Save(1),ROOT.RooFit.NumCPU(8))  
        
       
        
  
            
     print 
     print period+" Prefit nonRes pdf:"
     pdf1_nonres_shape_prefit = args["nonResNominal_JJ_"+purity+"_13TeV_"+period]
     pdf1_nonres_shape_prefit.Print()
     print
    
     
     print period+" Postfit nonRes pdf:"
     pdf1_nonres_shape_postfit  = args["shapeBkg_nonRes_JJ_"+purity+"_13TeV_2016"]
     pdf1_nonres_shape_postfit.Print()
     pdf1_nonres_shape_postfit.funcList().Print()
     pdf1_nonres_shape_postfit.coefList().Print()
     print
    
     
     print period+" Postfit W+jets res pdf:"
     pdf1_Wres_shape_postfit  = args["shapeBkg_Wjets_JJ_"+purity+"_13TeV_2016"]
     pdf1_Wres_shape_postfit.Print()
    
     
     print period+" Postfit Z+jets res pdf:"
     pdf1_Zres_shape_postfit  = args["shapeBkg_Zjets_JJ_"+purity+"_13TeV_2016"]
     pdf1_Zres_shape_postfit.Print()
     
     
     if options.addTop:
        print period+" Postfit tt res pdf:"
        pdf1_TThad_shape_postfit  = args["shapeBkg_TThad_JJ_"+purity+"_13TeV_2016"]
        pdf1_TThad_shape_postfit.Print()
        
     

     pdf1_signal_postfit = 0
     if options.fitSignal:
      print period+" Signal pdf:"     
      pdf1_signal_postfit  = args["shapeSig_"+signalName+"_JJ_"+purity+"_13TeV_2016"]
      pdf1_signal_postfit.Print()
      print
      
      
     print "Full 2016 post-fit pdf:"     
     pdf1_shape_postfit  = args[pdf1Name+"_nuis"]
     pdf1_shape_postfit.Print()
     print
    
		    
     allpdfsz = [] #let's have always pre-fit and post-fit as firt elements here, and add the optional shapes if you want with options.pdf
     allpdfsz.append(pdf1_nonres_shape_postfit)
     allpdfsz.append(pdf1_Wres_shape_postfit)
     allpdfsz.append(pdf1_Zres_shape_postfit)
     
     if options.addTop:
      print "add top"
      allpdfsz.append(pdf1_TThad_shape_postfit)
     else:
       allpdfsz.append(pdf1_Zres_shape_postfit)
     allpdfsz.append(pdf1_shape_postfit)
     for p in options.pdfz.split(","):
         if p == '': continue
         print "add pdf:",p
         args[p].Print()
         allpdfsz.append(args[p])

     allpdfsx = [] #let's have always pre-fit and post-fit as firt elements here, and add the optional shapes if you want with options.pdf
     allpdfsx.append(pdf1_nonres_shape_postfit)
     allpdfsx.append(pdf1_Wres_shape_postfit)
     allpdfsx.append(pdf1_Zres_shape_postfit)
     if options.addTop:
      allpdfsx.append(pdf1_TThad_shape_postfit)
     else:
       allpdfsx.append(pdf1_Zres_shape_postfit) #dummy, not used
     allpdfsx.append(pdf1_shape_postfit)
     for p in options.pdfx.split(","):
         if p == '': continue
	 print "add pdf:",p
	 args[p].Print()
         allpdfsx.append(args[p])

     allpdfsy = [] #let's have always pre-fit and post-fit as firt elements here, and add the optional shapes if you want with options.pdf
     allpdfsy.append(pdf1_nonres_shape_postfit)
     allpdfsy.append(pdf1_Wres_shape_postfit)
     allpdfsy.append(pdf1_Zres_shape_postfit)
     if options.addTop:
      allpdfsy.append(pdf1_TThad_shape_postfit)
     else:
      allpdfsy.append(pdf1_Zres_shape_postfit)
     allpdfsy.append(pdf1_shape_postfit)
     for p in options.pdfy.split(","):
         if p == '': continue
         print "add pdf:",p
         args[p].Print()
         allpdfsy.append(args[p])
      
      
     norm1_nonres = [0,0]
     norm1_nonres[0] = (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_nonRes"].getVal()
     norm1_Wres = [0,0]
     norm1_Wres[0] = (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_Wjets"].getVal()
     norm1_Zres = [0,0]
     norm1_Zres[0] = (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_Zjets"].getVal()
     norm1_TThad = [0,0]
     norm1_sig = [0,0]
     norm1_sig[0] = (args[pdf1Name].getComponents())["n_exp_final_binJJ_"+purity+"_13TeV_"+period+"_proc_"+signalName+""].getVal()
     if doFit: 
        print
        (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_nonRes"].dump()     
        norm1_nonres[1] = (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_nonRes"].getPropagatedError(fitresult)
        
                    
        print
        (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_Wjets"].dump()
        norm1_Wres[1] = (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_Wjets"].getPropagatedError(fitresult)
        
        print
        (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_Zjets"].dump()
        norm1_Zres[1] = (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_Zjets"].getPropagatedError(fitresult)
        
        
        if options.addTop:
            print
            (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_TThad"].dump()
            norm1_TThad[0] = (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_TThad"].getVal()
            norm1_TThad[1] = (args[pdf1Name].getComponents())["n_exp_binJJ_"+purity+"_13TeV_"+period+"_proc_TThad"].getPropagatedError(fitresult)
        
        else:
            norm1_TThad    = [0,0]
            norm1_TThad[0] = 1.
            norm1_TThad[1] = 1.
            norm2_TThad    = [0,0]
            norm2_TThad[0] = 1.
            norm2_TThad[1] = 1.
        
        norm1_sig = 0
        if options.fitSignal:
            print
            #(args[pdfName].getComponents())["n_exp_final_binJJ_"+purity+"_13TeV_proc_"+signalName+""].dump()
            norm1_sig[1] = (args[pdf1Name].getComponents())["n_exp_final_binJJ_"+purity+"_13TeV_"+period+"_proc_"+signalName+""].getPropagatedError(fitresult)
            
                    
        print
        print "QCD normalization after fit: ",norm1_nonres[0],"+/-",norm1_nonres[1],"("+period+")"
        print "W+jets normalization after fit: ",norm1_Wres[0],"+/-",norm1_Wres[1],"("+period+")"
        print "Z+jets normalization after fit: ",norm1_Zres[0],"+/-",norm1_Zres[1],"("+period+")"
        if options.addTop: print "tt normalization after fit: ",norm1_TThad[0],"+/-",norm1_TThad[1],"("+period+")"
        if options.fitSignal: print "Signal yields after fit: ",norm1_sig[0],"+/-",norm1_sig[1],"("+period+")"
        print
     
                   
    
     #################################################
     # get variables from workspace 
     MJ1= workspace.var("MJ1");
     MJ2= workspace.var("MJ2");
     MJJ= workspace.var("MJJ");
     #del workspace
    
     argset = ROOT.RooArgSet();
     argset.add(MJJ);
     argset.add(MJ2);
     argset.add(MJ1);
     args_ws = argset
     x = getListFromRange(options.xrange)
     y = getListFromRange(options.yrange)
     z = getListFromRange(options.zrange)     
     
     xBins_redux = reduceBinsToRange(xBins,x)
     yBins_redux = reduceBinsToRange(yBins,y)
     zBins_redux = reduceBinsToRange(zBins,z)
     print "x bins reduced:"
     print xBins_redux
     print "y bins reduced:"
     print yBins_redux
     print "z bins reduced:"
     print zBins_redux
     print 

     #################################################     
     #make projections onto MJJ axis
     if options.projection =="z": doZprojection(allpdfsz,data1,norm1_nonres,norm1_Wres,pdf1_signal_postfit,norm1_sig,norm1_Zres,norm1_TThad)
         
     #make projections onto MJ1 axis
     if options.projection =="x": doXprojection(allpdfsx,data1,norm1_nonres,norm1_Wres,pdf1_signal_postfit,norm1_sig,norm1_Zres,norm1_TThad)
                  
     #make projections onto MJ2 axis
     if options.projection =="y": doYprojection(allpdfsy,data1,norm1_nonres,norm1_Wres,pdf1_signal_postfit,norm1_sig,norm1_Zres,norm1_TThad)
         
     if options.projection =="xyz":
      doZprojection(allpdfsz,data1,norm1_nonres,norm1_Wres,pdf1_signal_postfit,norm1_sig,norm1_Zres,norm1_TThad)
      doXprojection(allpdfsx,data1,norm1_nonres,norm1_Wres,pdf1_signal_postfit,norm1_sig,norm1_Zres,norm1_TThad)
      doYprojection(allpdfsy,data1,norm1_nonres,norm1_Wres,pdf1_signal_postfit,norm1_sig,norm1_Zres,norm1_TThad)
     
  
