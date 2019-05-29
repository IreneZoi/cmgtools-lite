import os,sys,optparse,time
import ROOT as rt
rt.gROOT.ProcessLine(".x tdrstyle.cc")
import CMS_lumi
#setTDRStyle()
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptTitle(0)
rt.gROOT.SetBatch(True)
from time import sleep

# Run from command line with
#python Projections3DHisto.py --mc JJ_nonRes_VH_HPLP.root,nonRes -k JJ_2016_nonRes_3D_VH_HPLP.root,histo -o control-plots-HPLP-pythia
#old
#python Projections3DHisto.py --mc 2017/JJ_nonRes_HPLP.root,nonRes -k 2016/JJ_nonRes_3D_HPLP_2017_copy.root,histo -o control-plots-HPLP-pythia
#python Projections3DHisto.py --mc 2017/JJ_nonRes_HPHP.root,nonRes -k 2016/JJ_nonRes_3D_HPHP_2017_copy.root,histo -o control-plots-HPHP-pythia
#python Projections3DHisto.py --mc 2016/JJ_nonRes_LPLP_altshapeUp.root,nonRes -k 2016/JJ_nonRes_3D_LPLP_fixed.root,histo_altshapeUp -o control-plots-LPLPfixed-herwig
#python Projections3DHisto.py --mc 2016/JJ_nonRes_LPLP_altshapeUp.root,nonRes -k 2016/JJ_nonRes_3D_LPLP.root,histo_altshapeUp -o control-plots-LPLP-herwig
#python Projections3DHisto.py --mc 2016/JJ_nonRes_LPLP_altshapeUp.root,nonRes -k JJ_nonRes_3D_LPLP.root,histo -o control-plots-LPLPnew-herwig
#python Projections3DHisto.py --mc JJ_nonRes_LPLP_nominal.root,nonRes -k JJ_nonRes_3D_LPLP.root,histo -o control-plots-LPLP-pythia

def get_canvas(cname):

 #change the CMS_lumi variables (see CMS_lumi.py)
 CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
 CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
 CMS_lumi.writeExtraText = 1
 CMS_lumi.extraText = "Simulation"
 #need to be fixed
 CMS_lumi.lumi_sqrtS = "13 TeV (2016)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

 iPos = 11
 if( iPos==0 ): CMS_lumi.relPosX = 0.12

 H_ref = 600 
 W_ref = 800 
 W = W_ref
 H  = H_ref

 iPeriod = 0

 # references for T, B, L, R
 T = 0.08*H_ref
 B = 0.12*H_ref 
 L = 0.12*W_ref
 R = 0.04*W_ref

 canvas = rt.TCanvas(cname,cname,50,50,W,H)
 canvas.SetFillColor(0)
 canvas.SetBorderMode(0)
 canvas.SetFrameFillStyle(0)
 canvas.SetFrameBorderMode(0)
 canvas.SetLeftMargin( L/W )
 canvas.SetRightMargin( R/W )
 canvas.SetTopMargin( T/H )
 canvas.SetBottomMargin( B/H )
 canvas.SetTickx(0)
 canvas.SetTicky(0)
 
 return canvas

parser = optparse.OptionParser()
parser.add_option("--mc","--mc",dest="mc",help="File with mc events and histo name (separated by comma)",default='JJ_nonRes_HPHP_nominal.root,nonRes')
parser.add_option("-k","--kernel",dest="kernel",help="File with kernel and histo name (separated by comma)",default='JJ_nonRes_3D_HPHP.root,histo')
parser.add_option("-o","--outdir",dest="outdir",help="Output directory for plots",default='control-plots')
parser.add_option("-l","--label",dest="label",help="MC type label (Pythia8, Herwig, Madgraph, Powheg)",default='Pythia8')
(options,args) = parser.parse_args()

#void Projections3DHisto(std::string dataFile, std::string hdataName, std::string fitFile, std::string hfitName, std::string outDirName){

os.system('rm -rf %s'%options.outdir)
os.system('mkdir %s'%options.outdir)

kfile,kname = options.kernel.split(',')
fin = rt.TFile.Open(kfile,"READ")
hin = fin.Get(kname)
hin.Scale(1./hin.Integral())

MCfile,MCname = options.mc.split(',')
finMC = rt.TFile.Open(MCfile,"READ")
hinMC = finMC.Get(MCname)
hinMC.Scale(1./hinMC.Integral())

binsx = hin.GetNbinsX()
xmin = hin.GetXaxis().GetXmin()
xmax = hin.GetXaxis().GetXmax()
print "xmin",xmin,"xmax",xmax,"binsx",binsx

binsy = hin.GetNbinsY()
ymin = hin.GetYaxis().GetXmin()
ymax = hin.GetYaxis().GetXmax()
print "ymin",ymin,"ymax",ymax,"binsy",binsy

binsz = hin.GetNbinsZ()
zmin = hin.GetZaxis().GetXmin()
zmax = hin.GetZaxis().GetXmax()
print "zmin",zmin,"zmax",zmax,"binsz",binsz

hx = []
hy = []
hz = []
hxMC = []
hyMC = []
hzMC = []

pullsx = []
pullsy = []
pullsz = []

zbinMin = [1,1,hin.GetZaxis().FindBin(1530)+1,hin.GetZaxis().FindBin(2546)+1]
zbinMax = [binsz,hin.GetZaxis().FindBin(1530),hin.GetZaxis().FindBin(2546),binsz]
colors = [1,99,9,8,94]

scale = [1.,0.8,3.0,30.]

for i in range(4):

 print "Plotting mJ projections",zbinMin[i],zbinMax[i]
 pname = "px_%i"%i
 hx.append( hin.ProjectionX(pname,1,binsy,zbinMin[i],zbinMax[i]) )
 pname = "py_%i"%i
 hy.append( hin.ProjectionY(pname,1,binsx,zbinMin[i],zbinMax[i]) )
 pname = "px_MC%i"%i
 hxMC.append( hinMC.ProjectionX(pname,1,binsy,zbinMin[i],zbinMax[i]) )
 pname = "pyMC_%i"%i
 hyMC.append( hinMC.ProjectionY(pname,1,binsx,zbinMin[i],zbinMax[i]) )

 pname = "pullsx_%i"%i
 pullsx.append( rt.TH1F(pname,pname,40,-10,10) )
 pname = "pullsy_%i"%i
 pullsy.append( rt.TH1F(pname,pname,40,-10,10) )
        
for i in range(4):
 for b in range(1,binsx+1):
  if hxMC[i].GetBinContent(b) != 0: pullsx[i].Fill( (hxMC[i].GetBinContent(b)-hx[i].GetBinContent(b))/hxMC[i].GetBinError(b) )
  if hyMC[i].GetBinContent(b) != 0: pullsy[i].Fill( (hyMC[i].GetBinContent(b)-hy[i].GetBinContent(b))/hyMC[i].GetBinError(b) )

for i in range(4):

 hx[i].Scale(scale[i])
 hy[i].Scale(scale[i])
 hxMC[i].Scale(scale[i])
 hyMC[i].Scale(scale[i])
    
 hx[i].SetLineColor(colors[i])
 hx[i].SetMarkerColor(colors[i])
 hy[i].SetLineColor(colors[i])
 hy[i].SetMarkerColor(colors[i])
 hxMC[i].SetLineColor(colors[i])
 hxMC[i].SetMarkerColor(colors[i])
 hxMC[i].SetMarkerStyle(20)
 hxMC[i].SetMarkerSize(0.5)
 hyMC[i].SetLineColor(colors[i])
 hyMC[i].SetMarkerColor(colors[i]) 
 hyMC[i].SetMarkerStyle(20)
 hyMC[i].SetMarkerSize(0.5)
 
 pullsx[i].SetLineColor(colors[i])
 pullsx[i].SetLineWidth(2)
 pullsx[i].SetMarkerSize(0)
 pullsy[i].SetLineColor(colors[i])
 pullsy[i].SetLineWidth(2)
 pullsy[i].SetMarkerSize(0)
  
hx[0].SetMinimum(0)
hx[0].SetMaximum(0.03)
hy[0].SetMinimum(0)
hy[0].SetMaximum(0.03)

#leg = rt.TLegend(0.6,0.6,0.85,0.8)
leg = rt.TLegend(0.51,0.60,0.76,0.85)
leg.SetBorderSize(0)
leg.SetTextSize(0.035)
leg.AddEntry(hxMC[0],"Simulation (%s)"%options.label,"LP")
leg.AddEntry(hx[0],"Template","L")
for i in range(1,4):
 leg.AddEntry(hx[i],"%.1f < m_{jj} < %.1f TeV"%( hin.GetZaxis().GetBinLowEdge(zbinMin[i])/1000.,hin.GetZaxis().GetBinUpEdge(zbinMax[i])/1000.) )
 
cx = get_canvas("cx")
cx.cd()
for i in range(4):
 hx[i].Draw("HISTsame")
 hxMC[i].Draw("PEsame")
hx[0].GetXaxis().SetTitle("m_{jet1} (proj. x) [GeV]")
leg.Draw()

CMS_lumi.CMS_lumi(cx, 0, 11)
cx.cd()
cx.Update()
cx.RedrawAxis()
frame = cx.GetFrame()
frame.Draw()
cx.SaveAs(options.outdir+"/cx.png","pdf")

cy = get_canvas("cy")
cy.cd()
hy[0].GetXaxis().SetTitle("m_{jet2} (proj. y) [GeV]")
hy[0].GetXaxis().SetTitleSize(hx[0].GetXaxis().GetTitleSize())
hy[0].GetXaxis().SetTitleOffset(hx[0].GetXaxis().GetTitleOffset())
for i in range(4): 
 hy[i].Draw("HISTsame")
 hyMC[i].Draw("PEsame")
leg.Draw()

CMS_lumi.CMS_lumi(cy, 0, 11)
cy.cd()
cy.Update()
cy.RedrawAxis()
frame = cy.GetFrame()
frame.Draw()
cy.SaveAs(options.outdir+"/cy.png","pdf")

'''
labelsXY = ['All m_{jj} bins']
for i in range(1,4): labelsXY.append( "%.1f < m_{jj} < %.1f TeV"%( hin.GetZaxis().GetBinLowEdge(zbinMin[i])/1000.,hin.GetZaxis().GetBinUpEdge(zbinMax[i])/1000.) )

for i in range(4):

 pt = rt.TPaveText(0.1436782,0.7690678,0.4224138,0.8644068,"brNDC")
 pt.SetTextFont(42)
 pt.SetTextSize(0.042)
 pt.SetTextAlign(22)
 pt.SetFillColor(0)
 pt.SetBorderSize(1)
 pt.SetFillStyle(0)
 pt.SetLineWidth(2)
 pt.AddText(labelsXY[i])
  
 cname = "cpullsx_%i"%i
 cpullsx = get_canvas(cname)
 cpullsx.cd() 
 pullsx[i].Draw("PE")

 f = rt.TF1("func","gaus(0)",-10,10)
 f.SetParameter(0,10)
 f.SetParError(0,5)
 f.SetParameter(1,0)
 f.SetParError(1,0.5)
 f.SetParameter(2,4)
 f.SetParError(2,2)
 pullsx[i].Fit("func")
 
 pt.Draw()

 cpullsx.SaveAs(options.outdir+"/"+cpullsx.GetName()+".png","pdf")

for i in range(4):

 pt = rt.TPaveText(0.1436782,0.7690678,0.4224138,0.8644068,"brNDC")
 pt.SetTextFont(42)
 pt.SetTextSize(0.042)
 pt.SetTextAlign(22)
 pt.SetFillColor(0)
 pt.SetBorderSize(1)
 pt.SetFillStyle(0)
 pt.SetLineWidth(2)
 pt.AddText(labelsXY[i])
  
 cname = "cpullsy_%i"%i
 cpullsy = get_canvas(cname)
 cpullsy.cd() 
 pullsy[i].Draw("PE")

 f = rt.TF1("func","gaus(0)",-10,10)
 f.SetParameter(0,10)
 f.SetParError(0,5)
 f.SetParameter(1,0)
 f.SetParError(1,0.5)
 f.SetParameter(2,4)
 f.SetParError(2,2)
 pullsy[i].Fit("func")
 
 pt.Draw()

 cpullsy.SaveAs(options.outdir+"/"+cpullsy.GetName()+".png","pdf")
'''

#xbinMin[5] = {1,hin.GetXaxis().FindBin(55),hin.GetXaxis().FindBin(70),hin.GetXaxis().FindBin(100),hin.GetXaxis().FindBin(150)}
#xbinMax[5] = {binsx,hin.GetXaxis().FindBin(70),hin.GetXaxis().FindBin(100),hin.GetXaxis().FindBin(150),binsx}
xbinMin = [1,hin.GetXaxis().FindBin(55),hin.GetXaxis().FindBin(73)+1,hin.GetXaxis().FindBin(103)+1,hin.GetXaxis().FindBin(167)+1]
xbinMax = [binsx,hin.GetXaxis().FindBin(73),hin.GetXaxis().FindBin(103),hin.GetXaxis().FindBin(167),binsx]
#ybinMin = hin.GetXaxis().FindBin(55)
#ybinMax = hin.GetXaxis().FindBin(73)
#float scalez[5] = {1.,1.,0.1,0.01,0.001}
scalez = [1.,1.,0.1,0.01,0.001]

for i in range(5):

 print "Plotting mJJ projections",xbinMin[i],xbinMax[i]
 pname = "pz_%i"%i
 hz.append( hin.ProjectionZ(pname,xbinMin[i],xbinMax[i],xbinMin[i],xbinMax[i]) )
 pname = "pzMC_%i"%i
 hzMC.append( hinMC.ProjectionZ(pname,xbinMin[i],xbinMax[i],xbinMin[i],xbinMax[i]) )

 pname = "pullsz_%i"%i
 pullsz.append( rt.TH1F(pname,pname,40,-10,10) )
   
for i in range(5):
 for b in range(1,binsx+1):
  if hzMC[i].GetBinContent(b) != 0: pullsz[i].Fill( (hzMC[i].GetBinContent(b)-hz[i].GetBinContent(b))/hzMC[i].GetBinError(b) )

for i in range(5):
 
 hz[i].Scale(scalez[i])
 hzMC[i].Scale(scalez[i])
   
 hz[i].SetLineColor(colors[i])
 hz[i].SetMarkerColor(colors[i])
 hzMC[i].SetLineColor(colors[i])
 hzMC[i].SetMarkerColor(colors[i])
 hzMC[i].SetMarkerStyle(20)
 hzMC[i].SetMarkerSize(0.5)

 pullsz[i].SetLineColor(colors[i])
 pullsz[i].SetLineWidth(2)
 pullsz[i].SetMarkerSize(0)


#leg2 = rt.TLegend(0.6,0.6,0.85,0.85)
leg2 = rt.TLegend(0.51,0.65,0.76,0.90)
leg2.SetBorderSize(0)
leg2.SetTextSize(0.035)
leg2.AddEntry(hzMC[0],"Simulation (%s)"%options.label,"LP")
leg2.AddEntry(hz[0],"Template","L")
for i in range(1,5):
 leg2.AddEntry(hz[i],"%i < m_{jet} < %i GeV"%(  hin.GetXaxis().GetBinLowEdge(xbinMin[i]),hin.GetXaxis().GetBinUpEdge(xbinMax[i])) )
 #else: leg2.AddEntry(hz[i],"%i < m_{jet} < %i GeV"%(  hin.GetXaxis().GetBinLowEdge(xbinMin[i]),hin.GetXaxis().GetBinUpEdge(xbinMax[i])) )

cz = get_canvas("cz")
cz.SetLogy()
cz.cd()
hz[0].SetMinimum(1E-11)
hz[0].SetMaximum(50.0)
for i in range(5):
 hz[i].Draw("HISTsame")
 hzMC[i].Draw("PEsame")
hz[0].GetXaxis().SetTitle("m_{jj} (proj. z) [GeV]")
leg2.Draw()

CMS_lumi.CMS_lumi(cz, 0, 11)
cz.cd()
cz.Update()
cz.RedrawAxis()
frame = cz.GetFrame()
frame.Draw()
cz.SaveAs(options.outdir+"/cz.png","pdf")

labelsZ = ["All m_{jet} bins"]
for i in range(1,5):
 labelsZ.append("%i < m_{jet} < %i GeV"%(hin.GetXaxis().GetBinLowEdge(xbinMin[i]),hin.GetXaxis().GetBinUpEdge(xbinMax[i])))
''' 
for i in range(5):

 pt = rt.TPaveText(0.1436782,0.7690678,0.4224138,0.8644068,"brNDC")
 pt.SetTextFont(42)
 pt.SetTextSize(0.042)
 pt.SetTextAlign(22)
 pt.SetFillColor(0)
 pt.SetBorderSize(1)
 pt.SetFillStyle(0)
 pt.SetLineWidth(2)
 pt.AddText(labelsZ[i])
 

 
 cname = "cpullsz_%i"%i
 cpullsz = get_canvas(cname)
 cpullsz.cd() 
 pullsz[i].Draw("PE")

 f = rt.TF1("func","gaus(0)",-10,10)
 f.SetParameter(0,10)
 f.SetParError(0,5)
 f.SetParameter(1,0)
 f.SetParError(1,0.5)
 f.SetParameter(2,4)
 f.SetParError(2,2)
 pullsz[i].Fit("func")
 
 pt.Draw()

 cpullsz.SaveAs(options.outdir+"/"+cpullsz.GetName()+".png","pdf")
''' 

#xbinMin = [hin.GetXaxis().FindBin(173)+1]
#xbinMax = [binsx]
  
hin_PTUp = fin.Get("histo_PTUp")
hin_PTUp.Scale(1./hin_PTUp.Integral())
hin_PTDown = fin.Get("histo_PTDown")
hin_PTDown.Scale(1./hin_PTDown.Integral())
hin_OPTUp = fin.Get("histo_OPTUp")
hin_OPTUp.Scale(1./hin_OPTUp.Integral())
hin_OPTDown = fin.Get("histo_OPTDown")
hin_OPTDown.Scale(1./hin_OPTDown.Integral())
print "irene commented out the stuff concerning the alternative shape because hists are empty!!  "
#irene commented out alt shape because hists are empty!!
#hin_altshapeUp = fin.Get("histo_altshapeUp")
#hin_altshapeUp.Scale(1./hin_altshapeUp.Integral())
#hin_altshapeDown = fin.Get("histo_altshapeDown")
#hin_altshapeDown.Scale(1./hin_altshapeDown.Integral())
#hin_altshape2Up = fin.Get("histo_altshape2Up")
#hin_altshape2Up.Scale(1./hin_altshape2Up.Integral())
#hin_altshape2Down = fin.Get("histo_altshape2Down")
#hin_altshape2Down.Scale(1./hin_altshape2Down.Integral())
#end of irene intervention
#hin_altshape3Up = fin.Get("histo_altshape3Up")
#hin_altshape3Up.Scale(1./hin_altshape3Up.Integral())
#hin_altshape3Down = fin.Get("histo_altshape3Down")
#hin_altshape3Down.Scale(1./hin_altshape3Down.Integral())
hin_OPT3Up = fin.histo_OPT3Up
hin_OPT3Up.Scale(1./hin_OPT3Up.Integral())
hin_OPT3Down = fin.histo_OPT3Down
hin_OPT3Down.Scale(1./hin_OPT3Down.Integral())

hz_PTUp = hin_PTUp.ProjectionZ("pz_PTUp",xbinMin[0],xbinMax[0],xbinMin[0],xbinMax[0])
hz_PTUp.SetLineColor(rt.kMagenta)
hz_PTUp.Scale(1./hz_PTUp.Integral())
hz_PTDown = hin_PTDown.ProjectionZ("pz_PTDown",xbinMin[0],xbinMax[0],xbinMin[0],xbinMax[0])
hz_PTDown.SetLineColor(rt.kMagenta)
hz_PTDown.Scale(1./hz_PTDown.Integral())
hz_OPTUp = hin_OPTUp.ProjectionZ("pz_OPTUp",xbinMin[0],xbinMax[0],xbinMin[0],xbinMax[0])
hz_OPTUp.SetLineColor(210)
hz_OPTUp.Scale(1./hz_OPTUp.Integral())
hz_OPTDown = hin_OPTDown.ProjectionZ("pz_OPTDown",xbinMin[0],xbinMax[0],xbinMin[0],xbinMax[0])
hz_OPTDown.SetLineColor(210)
hz_OPTDown.Scale(1./hz_OPTDown.Integral())
#irene commented out alt shape because hists are empty!! 
#hz_altshapeUp = hin_altshapeUp.ProjectionZ("pz_altshapeUp",xbinMin[0],xbinMax[0],xbinMin[0],xbinMax[0])
#hz_altshapeUp.SetLineColor(rt.kBlue)
#hz_altshapeUp.Scale(1./hz_altshapeUp.Integral())
#hz_altshapeDown = hin_altshapeDown.ProjectionZ("pz_altshapeDown",xbinMin[0],xbinMax[0],xbinMin[0],xbinMax[0])
#hz_altshapeDown.SetLineColor(rt.kBlue)
#hz_altshapeDown.Scale(1./hz_altshapeDown.Integral())
#hz_altshape2Up = hin_altshape2Up.ProjectionZ("pz_altshape2Up",xbinMin[0],xbinMax[0],xbinMin[0],xbinMax[0])
#hz_altshape2Up.SetLineColor(rt.kRed)
#hz_altshape2Up.Scale(1./hz_altshape2Up.Integral())
#hz_altshape2Down = hin_altshape2Down.ProjectionZ("pz_altshape2Down",xbinMin[0],xbinMax[0],xbinMin[0],xbinMax[0])
#hz_altshape2Down.SetLineColor(rt.kRed)
#hz_altshape2Down.Scale(1./hz_altshape2Down.Integral())
#end of irene intervention
#hz_altshape3Up = hin_altshape3Up.ProjectionZ("pz_altshape3Up",xbinMin[0],xbinMax[0],xbinMin[0],xbinMax[0])
#hz_altshape3Up.SetLineColor(rt.kOrange+1)
#hz_altshape3Down = hin_altshape3Down.ProjectionZ("pz_altshape3Down",xbinMin[0],xbinMax[0],xbinMin[0],xbinMax[0])
#hz_altshape3Down.SetLineColor(rt.kOrange+1)
hz_OPT3Up = hin_OPT3Up.ProjectionZ("pz_OPT3Up",xbinMin[0],xbinMax[0],xbinMin[0],xbinMax[0])
hz_OPT3Up.SetLineColor(rt.kViolet-6)
hz_OPT3Up.Scale(1./hz_OPT3Up.Integral())
hz_OPT3Down = hin_OPT3Down.ProjectionZ("pz_OPT3Down",xbinMin[0],xbinMax[0],xbinMin[0],xbinMax[0])
hz_OPT3Down.SetLineColor(rt.kViolet-6)
hz_OPT3Down.Scale(1./hz_OPT3Down.Integral())

hzMC[0].Scale(1./hzMC[0].Integral())
hz[0].Scale(1./hz[0].Integral())
#leg3 = rt.TLegend(0.6,0.55,0.95,0.8)
leg3 = rt.TLegend(0.53,0.55,0.78,0.89)
leg3.SetBorderSize(0)
leg3.SetTextSize(0.035)
leg3.AddEntry(hzMC[0],"Simulation (%s)"%(options.label),"LP")
leg3.AddEntry(hz[0],"Template","L")
leg3.AddEntry(hz_PTUp,"#propto m_{jj} up/down","L")
leg3.AddEntry(hz_OPTUp,"#propto 1/m_{jj} up/down","L")
#irene commented out alt shape because hists are empty!!         
#leg3.AddEntry(hz_altshapeUp,"HERWIG up/down","L")
#leg3.AddEntry(hz_altshape2Up,"MADGRAPH+PYTHIA up/down","L")
#end of irene intervention
#leg3.AddEntry(hz_altshape3Up,"POWHEG up/down","L")
leg3.AddEntry(hz_OPT3Up,"m_{jj} turn-on up/down","L")

czSyst = get_canvas("czSyst")
czSyst.cd()
czSyst.SetLogy()

#hz[4].SetLineColor(rt.kBlack)
#hz[4].Scale(1./0.001)
hz[0].SetMinimum(1E-06)
hz[0].SetMaximum(10.0)
hz[0].Draw("HIST")
hz_PTUp.Draw("HISTsame")
hz_PTDown.Draw("HISTsame") 
hz_OPTUp.Draw("HISTsame")
hz_OPTDown.Draw("HISTsame")
#irene commented out alt shape because hists are empty!!         
#hz_altshapeUp.Draw("HISTsame")
#hz_altshapeDown.Draw("HISTsame")
#hz_altshape2Up.Draw("HISTsame")
#hz_altshape2Down.Draw("HISTsame")
#end of irene intervention  
#hz_altshape3Up.Draw("HISTsame")
#hz_altshape3Down.Draw("HISTsame")
#hzMC[4].SetLineColor(rt.kBlack)
#hzMC[4].SetMarkerColor(rt.kBlack)
#hzMC[4].Scale(1./0.001)
hz_OPT3Up.Draw("HISTsame")
hz_OPT3Down.Draw("HISTsame")
hzMC[0].Draw("same")
leg3.Draw()

CMS_lumi.CMS_lumi(czSyst, 0, 11)
czSyst.cd()
czSyst.Update()
czSyst.RedrawAxis()
frame = czSyst.GetFrame()
frame.Draw()
czSyst.SaveAs(options.outdir+"/czSyst.png","pdf")
# sleep(10000)
hx_PTUp = hin_PTUp.ProjectionX("px_PTUp",1,binsy,zbinMin[0],zbinMax[0])
hx_PTUp.SetLineColor(rt.kMagenta)
hx_PTUp.Scale(1./hx_PTUp.Integral())
hx_PTDown = hin_PTDown.ProjectionX("px_PTDown",1,binsy,zbinMin[0],zbinMax[0])
hx_PTDown.SetLineColor(rt.kMagenta)
hx_PTDown.Scale(1./hx_PTDown.Integral())
hx_OPTUp = hin_OPTUp.ProjectionX("px_OPTUp",1,binsy,zbinMin[0],zbinMax[0])
hx_OPTUp.SetLineColor(210)
hx_OPTUp.Scale(1./hx_OPTUp.Integral())
hx_OPTDown = hin_OPTDown.ProjectionX("px_OPTDown",1,binsy,zbinMin[0],zbinMax[0])
hx_OPTDown.SetLineColor(210)
hx_OPTDown.Scale(1./hx_OPTDown.Integral())
#irene commented out alt shape because hists are empty!!  
#hx_altshapeUp = hin_altshapeUp.ProjectionX("px_altshapeUp",1,binsy,zbinMin[0],zbinMax[0])
#hx_altshapeUp.SetLineColor(rt.kBlue)
#hx_altshapeUp.Scale(1./hx_altshapeUp.Integral())
#hx_altshapeDown = hin_altshapeDown.ProjectionX("px_altshapeDown",1,binsy,zbinMin[0],zbinMax[0])
#hx_altshapeDown.SetLineColor(rt.kBlue)
#hx_altshapeDown.Scale(1./hx_altshapeDown.Integral())
#hx_altshape2Up = hin_altshape2Up.ProjectionX("px_altshape2Up",1,binsy,zbinMin[0],zbinMax[0])
#hx_altshape2Up.SetLineColor(rt.kRed)
#hx_altshape2Up.Scale(1./hx_altshape2Up.Integral())
#hx_altshape2Down = hin_altshape2Down.ProjectionX("px_altshape2Down",1,binsy,zbinMin[0],zbinMax[0])
#hx_altshape2Down.SetLineColor(rt.kRed)
#hx_altshape2Down.Scale(1./hx_altshape2Down.Integral())
#end of irene intervention  
#hx_altshape3Up = hin_altshape3Up.ProjectionX("px_altshape3Up",1,binsy,zbinMin[0],zbinMax[0])
#hx_altshape3Up.SetLineColor(rt.kOrange+1)
#hx_altshape3Down = hin_altshape3Down.ProjectionX("px_altshape3Down",1,binsy,zbinMin[0],zbinMax[0])
#hx_altshape3Down.SetLineColor(rt.kOrange+1)
hx_OPT3Up = hin_OPT3Up.ProjectionX("px_OPT3Up",1,binsy,zbinMin[0],zbinMax[0])
hx_OPT3Up.SetLineColor(rt.kViolet-6)
hx_OPT3Up.Scale(1./hx_OPT3Up.Integral())
hx_OPT3Down = hin_OPT3Down.ProjectionX("px_OPT3Down",1,binsy,zbinMin[0],zbinMax[0])
hx_OPT3Down.SetLineColor(rt.kViolet-6)
hx_OPT3Down.Scale(1./hx_OPT3Down.Integral())


hxMC[0].Scale(1./hxMC[0].Integral())
hx[0].Scale(1./hx[0].Integral())
#leg3 = rt.TLegend(0.6,0.55,0.95,0.8)
leg3 = rt.TLegend(0.53,0.50,0.78,0.84)
leg3.SetBorderSize(0)
leg3.SetTextSize(0.035)
leg3.AddEntry(hxMC[0],"Simulation (%s)"%(options.label),"LP")
leg3.AddEntry(hx[0],"Template","L")
leg3.AddEntry(hx_PTUp,"#propto m_{jj} up/down","L")
leg3.AddEntry(hx_OPTUp,"#propto 1/m_{jj} up/down","L")
#irene commented out alt shape because hists are empty!!  
#leg3.AddEntry(hx_altshapeUp,"HERWIG up/down","L")
#leg3.AddEntry(hx_altshape2Up,"MADGRAPH+PYTHIA up/down","L")
#end of irene intervention
#leg3.AddEntry(hx_altshape3Up,"POWHEG up/down","L")
leg3.AddEntry(hx_OPT3Up,"m_{jj} turn-on up/down","L")

cxSyst = get_canvas("cxSyst")
cxSyst.cd()

hx[0].SetMinimum(0)
hx[0].SetMaximum(0.04)
hx[0].Draw("HIST")
hx_PTUp.Draw("HISTsame")
hx_PTDown.Draw("HISTsame") 
hx_OPTUp.Draw("HISTsame")
hx_OPTDown.Draw("HISTsame")
#irene commented out alt shape because hists are empty!!  
#hx_altshapeUp.Draw("HISTsame")
#hx_altshapeDown.Draw("HISTsame")
#hx_altshape2Up.Draw("HISTsame")
#hx_altshape2Down.Draw("HISTsame")
#end of irene intervention
#hx_altshape3Up.Draw("HISTsame")
#hx_altshape3Down.Draw("HISTsame")
hx_OPT3Up.Draw("HISTsame")
hx_OPT3Down.Draw("HISTsame")
hxMC[0].Draw("same")
leg3.Draw()

CMS_lumi.CMS_lumi(cxSyst, 0, 11)
cxSyst.cd()
cxSyst.Update()
cxSyst.RedrawAxis()
frame = cxSyst.GetFrame()
frame.Draw()
cxSyst.SaveAs(options.outdir+"/cxSyst.png","pdf")


hy_PTUp = hin_PTUp.ProjectionY("py_PTUp",1,binsy,zbinMin[0],zbinMax[0])
hy_PTUp.SetLineColor(rt.kMagenta)
hy_PTUp.Scale(1./hy_PTUp.Integral())
hy_PTDown = hin_PTDown.ProjectionY("py_PTDown",1,binsy,zbinMin[0],zbinMax[0])
hy_PTDown.SetLineColor(rt.kMagenta)
hy_PTDown.Scale(1./hy_PTDown.Integral())
hy_OPTUp = hin_OPTUp.ProjectionY("py_OPTUp",1,binsy,zbinMin[0],zbinMax[0])
hy_OPTUp.SetLineColor(210)
hy_OPTUp.Scale(1./hy_OPTUp.Integral())
hy_OPTDown = hin_OPTDown.ProjectionY("py_OPTDown",1,binsy,zbinMin[0],zbinMax[0])
hy_OPTDown.SetLineColor(210)
hy_OPTDown.Scale(1./hy_OPTDown.Integral())
#irene commented out alt shape because hists are empty!!  
#hy_altshapeUp = hin_altshapeUp.ProjectionY("py_altshapeUp",1,binsy,zbinMin[0],zbinMax[0])
#hy_altshapeUp.SetLineColor(rt.kBlue)
#hy_altshapeUp.Scale(1./hy_altshapeUp.Integral())
#hy_altshapeDown = hin_altshapeDown.ProjectionY("py_altshapeDown",1,binsy,zbinMin[0],zbinMax[0])
#hy_altshapeDown.SetLineColor(rt.kBlue)
#hy_altshapeDown.Scale(1./hy_altshapeDown.Integral())
#hy_altshape2Up = hin_altshape2Up.ProjectionY("py_altshape2Up",1,binsy,zbinMin[0],zbinMax[0])
#hy_altshape2Up.SetLineColor(rt.kRed)
#hy_altshape2Up.Scale(1./hy_altshape2Up.Integral())
#hy_altshape2Down = hin_altshape2Down.ProjectionY("py_altshape2Down",1,binsy,zbinMin[0],zbinMax[0])
#hy_altshape2Down.SetLineColor(rt.kRed)
#hy_altshape2Down.Scale(1./hy_altshape2Down.Integral())
#end of irene intervention 
#hy_altshape3Up = hin_altshape3Up.ProjectionY("py_altshape3Up",1,binsy,zbinMin[0],zbinMax[0])
#hy_altshape3Up.SetLineColor(rt.kOrange+1)
#hy_altshape3Down = hin_altshape3Down.ProjectionY("py_altshape3Down",1,binsy,zbinMin[0],zbinMax[0])
#hy_altshape3Down.SetLineColor(rt.kOrange+1)
hy_OPT3Up = hin_OPT3Up.ProjectionY("py_OPT3Up",xbinMin[0],xbinMax[0],zbinMin[0],zbinMax[0])
hy_OPT3Up.SetLineColor(rt.kViolet-6)
hy_OPT3Up.Scale(1./hy_OPT3Up.Integral())
hy_OPT3Down = hin_OPT3Down.ProjectionY("py_OPT3Down",xbinMin[0],xbinMax[0],zbinMin[0],zbinMax[0])
hy_OPT3Down.SetLineColor(rt.kViolet-6)
hy_OPT3Down.Scale(1./hy_OPT3Down.Integral())

hyMC[0].Scale(1./hyMC[0].Integral())
hy[0].Scale(1./hy[0].Integral())
#leg3 = rt.TLegend(0.6,0.55,0.95,0.8)
leg3 = rt.TLegend(0.53,0.50,0.78,0.84)
leg3.SetBorderSize(0)
leg3.SetTextSize(0.035)
leg3.AddEntry(hyMC[0],"Simulation (%s)"%(options.label),"LP")
leg3.AddEntry(hy[0],"Template","L")
leg3.AddEntry(hy_PTUp,"#propto m_{jj} up/down","L")
leg3.AddEntry(hy_OPTUp,"#propto 1/m_{jj} up/down","L")
#irene commented out alt shape because hists are empty!! 
#leg3.AddEntry(hy_altshapeUp,"HERWIG up/down","L")
#leg3.AddEntry(hy_altshape2Up,"MADGRAPH+PYTHIA up/down","L")
#end of irene intervention 
#leg3.AddEntry(hy_altshape3Up,"POWHEG up/down","L")
leg3.AddEntry(hy_OPT3Up,"m_{jj} turn-on up/down","L")

cySyst = get_canvas("cySyst")
cySyst.cd()

hy[0].SetMinimum(0)
hy[0].SetMaximum(0.04)
hy[0].Draw("HIST")
hy_PTUp.Draw("HISTsame")
hy_PTDown.Draw("HISTsame") 
hy_OPTUp.Draw("HISTsame")
hy_OPTDown.Draw("HISTsame")
#irene commented out alt shape because hists are empty!! 
#hy_altshapeUp.Draw("HISTsame")
#hy_altshapeDown.Draw("HISTsame")
#hy_altshape2Up.Draw("HISTsame")
#hy_altshape2Down.Draw("HISTsame")
#end of irene intervention 
#hy_altshape3Up.Draw("HISTsame")
#hy_altshape3Down.Draw("HISTsame")
hy_OPT3Up.Draw("HISTsame")
hy_OPT3Down.Draw("HISTsame")
hyMC[0].Draw("same")
leg3.Draw()

CMS_lumi.CMS_lumi(cySyst, 0, 11)
cySyst.cd()
cySyst.Update()
cySyst.RedrawAxis()
frame = cySyst.GetFrame()
frame.Draw()
cySyst.SaveAs(options.outdir+"/cySyst.png","pdf")



'''
TCanvas* cxz = new TCanvas("cxz","cxz")
cxz.cd()
TH2F* hxz = (TH2F*)hin.Project3D("zx")
hxz.Draw("COLZ")

cxz.SaveAs(TString(outDirName)+TString("/")+TString("cxz.png"),"pdf")

TCanvas* cyz = new TCanvas("cyz","cyz")
cyz.cd()
TH2F* hyz = (TH2F*)hin.Project3D("zy")
hyz.Draw("COLZ")

cyz.SaveAs(TString(outDirName)+TString("/")+TString("cyz.png"),"pdf")

}
'''
