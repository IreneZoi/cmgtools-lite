import ROOT
import os,sys
from  CMGTools.VVResonances.plotting.CMS_lumi import *
#import CMS_lumi
#ROOT.gROOT.ProcessLine(".x tdrstyle.cc")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLegendBorderSize(0)

# usage:
#    python makeProjections.py purity year
# where purity can be "VV_HPHP", "VH_HPLP" and so on
# and year can be 2016, 2017 or 2018
path="../../plots/"

def getCanvas(w=800,h=600):
 H_ref = 600 
 W_ref = 600 
 W = W_ref
 H  = H_ref

 iPeriod = 0

 # references for T, B, L, R
 T = 0.08*H_ref
 B = 0.12*H_ref 
 L = 0.13*W_ref
 R = 0.04*W_ref
 cname = "c"
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


def GetName(histname):
    if histname.find("OPT3")!=-1:
        return "m_{jj} turn-on up/down"
    if histname.find("OPT")!=-1:
        return "#propto 1/m_{jj} up/down"
    if histname.find("PT")!=-1:
        return "#propto m_{jj} up/down"
    if histname.find("altshape2")!=-1:
        return "HERWIG up/down"
    if histname.find("altshape")!=-1:
        return "MADGRAPH+PYTHIA up/down"
    else:
        return "Template"

def plotStack(backgrounds,xrange,yrange,zrange,projection,purity):
    c = getCanvas()
    #leg = ROOT.TLegend(0.12,0.91,0.3,0.7)
    leg = ROOT.TLegend(0.6,0.89,0.87,0.7)
    if projection == "z":
        leg = ROOT.TLegend(0.6,0.89,0.89,0.7)
    leg.SetTextFont(42)
    
    p_bkg =[]
    labels  =[]
    stack = ROOT.THStack("stackplot","stackplot")
    projections=[]
    for bkg in backgrounds:
        print " make projection of histogram "+str(bkg.GetName())
        if projection == "x":
            print "zrange "+str(zrange[0][1]) + " to "+str(zrange[0][2])
            p = makeProjections(bkg,yrange[0],zrange[0],projection,bkg.GetName().replace("_all",""))
        if projection == "y":
            p = makeProjections(bkg,xrange[0],zrange[0],projection,bkg.GetName().replace("_all",""))
        if projection == "z":
            p = makeProjections(bkg,xrange[0],yrange[0],projection,bkg.GetName().replace("_all",""))
        p_bkg.append(p) 
        if bkg.GetName().find("W")!=-1:
            p.SetFillColor(ROOT.kRed)
            labels.append("W+jets")
        if bkg.GetName().find("Z")!=-1:
            p.SetFillColor(ROOT.kGreen)
            labels.append("Z+jets")
        if bkg.GetName().find("tt")!=-1:
            p.SetFillColor(ROOT.kBlue)
            labels.append("t#bar{t}+jets")
        if bkg.GetName().find("nonRes")!=-1:
            p.SetFillColor(ROOT.kGray)
            labels.append("qcd multijet")
        projections.append(p)
        stack.Add(p)
        
    for i in range(len(backgrounds)):
        leg.AddEntry(projections[len(backgrounds)-1-i],labels[len(backgrounds)-1-i],"f")
    
    addText = ROOT.TLatex()
    addText.SetTextFont(42)
    label_proj_2 = str(yrange[0][1])+" < m_{jet2} < "+str(yrange[0][2])
    label_proj_1 = str(xrange[0][1])+" < m_{jet1} < "+str(xrange[0][2])
    label_proj_3 = str(zrange[0][1])+" < m_{jj} < "    +str(zrange[0][2])+" GeV"

    
    stack.Draw("hist")
    stack.GetXaxis().SetTitle("Softdrop m_{jet1} [GeV]")
    stack.GetYaxis().SetTitle("Events / 2 GeV")
    stack.GetYaxis().SetTitleOffset(1.35)
    stack.GetXaxis().SetTitleOffset(1.05)
    stack.GetXaxis().SetNdivisions(5,5,0)
    stack.GetXaxis().SetLabelSize(0.05)
    stack.GetYaxis().SetLabelSize(0.05)
    stack.GetXaxis().SetTitleSize(0.05)
    stack.GetYaxis().SetTitleSize(0.05)
    
    leg.Draw("same")
    if projection=="z":
        c.SetLogy()
        addText.DrawLatexNDC(0.55,0.65,"#scale[0.8]{"+label_proj_1+"}")
        addText.DrawLatexNDC(0.55,0.6,"#scale[0.8]{"+label_proj_2+"}")
        stack.GetXaxis().SetTitle("m_{jj} [GeV]")
        stack.GetYaxis().SetTitle("Events / 100 GeV")
        #stack.GetYaxis().SetRangeUser(0,10e5)
    else:
        if projection =="y":
            stack.GetXaxis().SetTitle("Softdrop m_{jet2} [GeV]")
            addText.DrawLatexNDC(0.55,0.65,"#scale[0.8]{"+label_proj_1+"}")
            addText.DrawLatexNDC(0.55,0.6,"#scale[0.8]{"+label_proj_3+"}")
        #addText.DrawLatexNDC(0.12,0.65,"#scale[0.8]{"+label_proj_2+"}")
        #addText.DrawLatexNDC(0.12,0.6,"#scale[0.8]{"+label_proj_3+"}")
        else:
            addText.DrawLatexNDC(0.55,0.65,"#scale[0.8]{"+label_proj_2+"}")
            addText.DrawLatexNDC(0.55,0.6,"#scale[0.8]{"+label_proj_3+"}")
    c.SaveAs(path+"CP_background_p"+projection+"_xrange"+str(xrange[0][1])+"-"+str(xrange[0][2])+"_yrange"+str(yrange[0][1])+"-"+str(yrange[0][2])+"zrange"+str(zrange[0][1])+"-"+str(zrange[0][2])+"_"+purity+".pdf")



def getLegend(x1=0.410010112,y1=0.523362,x2=0.55202143,y2=0.879833):
  legend = ROOT.TLegend(x1,y1,x2,y2)
  legend.SetTextSize(0.04)
  legend.SetLineColor(0)
  legend.SetShadowColor(0)
  legend.SetLineStyle(1)
  legend.SetLineWidth(1)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.SetTextFont(42)
  legend.SetMargin(0.35)
  return legend

def plotHistos(xrange,yrange,zrange,hists,outname,label,dolog=False,data=None,scales=None,p=""):
    colors=[ROOT.kBlack,ROOT.kBlue,ROOT.kBlue,ROOT.kRed,ROOT.kRed,ROOT.kGreen+2,ROOT.kGreen+2,ROOT.kOrange,ROOT.kOrange,ROOT.kGray+2,ROOT.kGray+2,ROOT.kGreen+2,ROOT.kPink]
    if scales!=None:
        colors=[ROOT.kBlack,ROOT.kBlue,ROOT.kRed,ROOT.kGreen+2,ROOT.kOrange,ROOT.kGray+2,ROOT.kGray+2,ROOT.kGreen+2,ROOT.kPink]
    
    if scales!=None:
        for s in range(0,len(scales)):
            hists[s].Scale(scales[s])
            data[s].Scale(scales[s])
    else:
        for h in hists:
            h.Scale(1/h.Integral())
    
    c= getCanvas()
    #leg = ROOT.TLegend(0.12,0.91,0.4,0.6)
    #if dolog == True:
    leg = getLegend()
    
    
    if data!=None:
        i=0
        for d in data:
            d.SetLineColor(ROOT.kBlack)
            d.SetMarkerColor(ROOT.kBlack)
            d.SetMarkerStyle(1)
            if i==0:
                leg.AddEntry(d,"Simulation (Pythia8)","ep")
            i+=1
    
    
    i=0
    for h in hists:
        h.SetMaximum(0.04)
        if outname.find("HPLP")!=-1:
            h.SetMaximum(0.04)
            h.SetMinimum(1e-6)
        if dolog == True:
            h.SetMaximum(10)
            
        if scales!=None:
            h.SetMaximum(0.3)
        h.SetLineColor(colors[i])
        #h.SetLineWidth(2)
        if i==0:
            leg.AddEntry(h,label[0],"l")
        else:
            if label[i]!=label[i-1]:
                leg.AddEntry(h,label[i],"l")
        h.Draw("histsame")
        i+=1
    leg.Draw("same")
    if dolog:
        c.SetLogy()
       
       
    if data!=None:
        i=0
        for d in data:
            d.SetLineColor(colors[i])
            d.SetMarkerColor(colors[i])
            if scales==None:
                d.Scale(1/d.Integral())
                d.SetLineColor(ROOT.kBlack)
                d.SetMarkerColor(ROOT.kBlack)

            #d.SetMarkerStyle(2)
            print "draw data points "
            d.Draw("PEsame")
            #if i==0:
            #    leg.AddEntry(d,"Simulation (PYTHIA)","ep")
            i+=1
            
   
    addText = ROOT.TLatex()
    addText.SetTextFont(42)
    if p =="NL":
        label_proj_2 = "HPLP"
        label_proj_1 = ""
        label_proj_3 = ""
        addText.DrawLatexNDC(0.45,0.25,"#scale[0.8]{"+label_proj_2+"}")
    else:
        label_proj_2 = str(yrange[1])+" GeV"+" < m_{jet2} < "+str(yrange[2])+" GeV"
        label_proj_1 = str(xrange[1])+" GeV"+" < m_{jet1} < "+str(xrange[2])+" GeV"
        label_proj_3 = str(zrange[1])+" GeV"+" < m_{jj} < "    +str(zrange[2])+" GeV"
        addText.DrawLatexNDC(0.55,0.55,"#scale[0.8]{"+label_proj_1+"}")
        addText.DrawLatexNDC(0.55,0.6,"#scale[0.8]{"+label_proj_2+"}")
        addText.DrawLatexNDC(0.55,0.65,"#scale[0.8]{"+label_proj_3+"}")        
    
    
    CMS_lumi.lumi_sqrtS = "13 TeV (2017)"
    CMS_lumi.CMS_lumi(c, 0, 10)
    c.SaveAs(outname.replace(".pdf","_x"+str(xrange[1])+"-"+str(xrange[2])+"_y"+str(yrange[1])+"-"+str(yrange[2])+"_z"+str(zrange[1])+"-"+str(zrange[2]) +".pdf"))
    
  
  
  
  
  
  
def findBins(histo,axis,r):
    if r == 0:
        return 0
    if r == -1:
        return -1
    if axis=="x":
        p = histo.ProjectionX("tmp")
        N = histo.GetNbinsX()
    if axis=="y":
        p = histo.ProjectionY("tmp")
        N = histo.GetNbinsY()
    if axis=="z":
        p = histo.ProjectionZ("tmp")
        N = histo.GetNbinsZ()
        
        
    for i in range(0,N+1):
        if ( r >= p.GetBinLowEdge(N-i)):
           #print "r "+str(r)
           #print "low edge" +str( p.GetBinLowEdge(N-i))
           #print "bin "+str( N-i)
           #print "N "+str(N)
           #print "i "+str(i)
           return N-i
    return 0


def makeProjections(histo,range1,range2,axis,label):
    if axis=="x":
        #print findBins(histo,"z",range2[2])
        #print range2[1]
        #print findBins(histo,"z",range2[1])
        p = histo.ProjectionX(label,findBins(histo,"y",range1[1]),findBins(histo,"y",range1[2]),findBins(histo,"z",range2[1]),findBins(histo,"z",range2[2]))
        p.GetXaxis().SetTitle("Softdrop m_{jet1} [GeV]")
        p.GetYaxis().SetTitleOffset(1.5)
        p.GetXaxis().SetTitleOffset(1.1)
        p.GetYaxis().SetTitle("A.U.")
        p.GetXaxis().SetNdivisions(5,5,0)
        p.GetXaxis().SetTitleSize(0.05)
        p.GetXaxis().SetLabelSize(0.05)
        p.GetYaxis().SetLabelSize(0.05)
        p.GetYaxis().SetTitleSize(0.05)
        return p
    if axis=="y":
        p = histo.ProjectionY(label,findBins(histo,"x",range1[1]),findBins(histo,"x",range1[2]),findBins(histo,"z",range2[1]),findBins(histo,"z",range2[2]))
        p.GetXaxis().SetTitle("Softdrop m_{jet2} [GeV]")
        p.GetYaxis().SetTitleOffset(1.5)
        p.GetXaxis().SetTitleOffset(1.1)
        p.GetYaxis().SetTitle("A.U.")
        p.GetXaxis().SetNdivisions(5,5,0)
        p.GetXaxis().SetTitleSize(0.05)
        p.GetYaxis().SetTitleSize(0.05)
        p.GetXaxis().SetLabelSize(0.05)
        p.GetYaxis().SetLabelSize(0.05)
        
        return p
    if axis== "z":
        #print "z projection "
        p = histo.ProjectionZ(label,findBins(histo,"x",range1[1]),findBins(histo,"x",range1[2]),findBins(histo,"y",range2[1]),findBins(histo,"y",range2[2]))
        p.GetXaxis().SetTitle(" m_{jj} [GeV]")
        p.GetYaxis().SetTitleOffset(1.5)
        p.GetYaxis().SetTitle("A.U.")
        p.GetXaxis().SetNdivisions(5,5,0)
        p.GetXaxis().SetTitleOffset(1.1)
        p.GetXaxis().SetTitleSize(0.05)
        p.GetYaxis().SetTitleSize(0.05)
        p.GetXaxis().SetLabelSize(0.05)
        p.GetYaxis().SetLabelSize(0.05)
        return p
    
    


#z-Axis : MVV spectrum   => projection X gives the projection of m_jet along the MVV axis
#y-Axis : mjet spectrum  => projection Y gives the projection of MVV along the m_jet dimension 
#x-Axis : mjet spectrum  => projection Y gives the projection of MVV along the m_jet dimension 



if __name__=="__main__":
    '''
# irene trying to make things working without QCD!
    f_kernel = ROOT.TFile(sys.argv[1],"READ")
    f_sample = ROOT.TFile(sys.argv[2],"READ")
    print f_kernel
    print f_sample  
    
    purity = "HPHP"
    if sys.argv[1].find("HPLP")!=-1:
        purity = "HPLP"
    if sys.argv[1].find("LPLP")!=-1:
        purity = "LPLP"
  
    kernel = f_kernel.Get("histo")
    sample = f_sample.Get("nonRes")
    #sample = f_sample.Get("data")
    print sample.Integral()
    print kernel.Integral()
    kernel.Scale(sample.Integral()/kernel.Integral())
    
   #  sample.SetFillColor(0)
    
   #  xrange = [["px0",55,215],["px1",65,85],["px2",85,105],["px3",105,215]]
   #  yrange = [["py0",55,215],["py1",65,85],["py2",85,105],["py3",105,215]]
   #  zrange = [["pz0",1126,5500],["pz1",1126,1300],["pz2",1300,2000],["pz3",2000,5000]]
    
   #  projections_x =[]
   #  data_x =[]
   #  labels_x =[]
   #  for z in zrange:
   #      p = makeProjections(kernel,["py",55,215],z,"x",z[0])
   #      projections_x.append(p)
   #      labels_x.append(str(z[1])+" < m_{jj} < "+str(z[2]))
    
   #      d = makeProjections(sample,["py",55,215],z,"x",z[0]+"_d")
   #      data_x.append(d)
    
   #  scales = [10,10,8,50.]
   #  plotHistos(xrange[0],yrange[0],zrange[0],projections_x,"x_projection"+"_"+purity+".pdf",labels_x,False,data_x,scales)
    
    
    
   #  projections_y =[]
   #  data_y =[]
   #  labels_y =[]
   #  for z in zrange:
   #      p = makeProjections(kernel,["px",55,215],z,"y",z[0])
   #      projections_y.append(p)
   #      labels_y.append(str(z[1])+" < m_{jj} < "+str(z[2]))
    
   #      d = makeProjections(sample,["px",55,215],z,"y",z[0]+"_d")
   #      data_y.append(d)
    
   #  scales = [10,10,8,50.]
   #  plotHistos(xrange[0],yrange[0],zrange[0],projections_y,"y_projection"+"_"+purity+".pdf",labels_y,False,data_y,scales)
    
    
   #  projections_z =[]
   #  data_z =[]
   #  labels_z =[]
   #  for j in range(0,len(xrange)):
   #      p = makeProjections(kernel,xrange[j],yrange[j],"z",zrange[j][0])
   #      projections_z.append(p)
   #      labels_z.append(str(xrange[j][1])+" < #rho_{jet1/2} < "+str(xrange[j][2]))
    
   #      d = makeProjections(sample,xrange[j],yrange[j],"z",zrange[j][0]+"_d")
   #      data_z.append(d)
   #  print projections_z
   #  scales = [1,1,1,0.5]
   #  plotHistos(xrange[0],yrange[0],zrange[0],projections_z,"z_projection"+"_"+purity+".pdf",labels_z,True,data_z,scales)
    
        
        
   # ############### make controlplots of QCD alternative shapes ####################
    #kernel_OPTXYup = f_kernel.Get("histo_OPTXYUp")
    kernel = f_kernel.Get("histo")
    #kernel_OPTXYdown = f_kernel.Get("histo_OPTXYDown")
    
    #kernel_PTup = f_kernel.Get("histo_PTXYUp")
    #kernel_PTdown = f_kernel.Get("histo_PTXYDown")
    
    kernel_PTup   = f_kernel.Get("histo_PTUp")
    kernel_PTdown = f_kernel.Get("histo_PTDown")
    
    kernel_OPTup   = f_kernel.Get("histo_OPTUp")
    kernel_OPTdown = f_kernel.Get("histo_OPTDown")
    
    kernel_altshapeUp   = f_kernel.Get("histo_altshapeUp")
    kernel_altshape2Up   = f_kernel.Get("histo_altshape2Up")
    
    kernel_altshapeDown   = f_kernel.Get("histo_altshapeDown")
    kernel_altshape2Down   = f_kernel.Get("histo_altshape2Down")
    
    kernel_OPT3Up   = f_kernel.Get("histo_OPT3Up")
    kernel_OPT3Down   = f_kernel.Get("histo_OPT3Down")
    
    
    xrange = [["px0",55,215]]
    yrange = [["py0",55,215]]
    zrange = [["pz0",1126,5000]]
    
    p_XY=[]
    XY =[kernel,kernel_OPTup,kernel_OPTdown,kernel_PTup,kernel_PTdown,kernel_altshapeUp,kernel_altshapeDown, kernel_altshape2Up,kernel_altshape2Down,kernel_OPT3Up,kernel_OPT3Down]
    labels_XY =[] 
    data_x = []
    for hist in XY:
        print " make projection of histogram "+str(hist.GetName())
        p = makeProjections(hist,yrange[0],zrange[0],"x",hist.GetName())    
        p_XY.append(p)
        d = makeProjections(sample,xrange[0],zrange[0],"x",zrange[0][0]+"_dx")
        data_x.append(d)
        labels_XY.append(GetName(hist.GetName().replace("histo_","")))

    plotHistos(xrange[0],yrange[0],zrange[0],p_XY,"qcd_altenate_shapes_X"+"_"+purity+".pdf",labels_XY,False,data_x,None,"NL") 
    p_XY=[]
    labels_XY =[] 
    data_y = []
    for hist in XY:
        print " make projection of histogram "+str(hist.GetName())
        p = makeProjections(hist,xrange[0],zrange[0],"y",hist.GetName())    
        p_XY.append(p)
        d = makeProjections(sample,xrange[0],zrange[0],"y",zrange[0][0]+"_dy")
        data_y.append(d)
        labels_XY.append(GetName(hist.GetName().replace("histo_","")))

    plotHistos(xrange[0],yrange[0],zrange[0],p_XY,"qcd_altenate_shapes_Y"+"_"+purity+".pdf",labels_XY,False,data_y,None,"NL") 
    
    
    
    p_Z=[]
    labels_Z =[] 
    data_z = []
    Z =[kernel,kernel_OPTup,kernel_OPTdown,kernel_PTup,kernel_PTdown,kernel_altshapeUp,kernel_altshapeDown, kernel_altshape2Up,kernel_altshape2Down,kernel_OPT3Up,kernel_OPT3Down]
    for hist in Z:
        hist.Scale(1/hist.Integral())
        print " make projection of histogram "+str(hist.GetName())
        p = makeProjections(hist,xrange[0],yrange[0],"z",hist.GetName())    
        p_Z.append(p)  
        d = makeProjections(sample,xrange[0],yrange[0],"z",zrange[0][0]+"_dz")
        data_z.append(d)
        labels_Z.append(GetName(hist.GetName().replace("histo_","")))

    plotHistos(xrange[0],yrange[0],zrange[0],p_Z,"qcd_altenate_shapes_Z"+"_"+purity+".pdf",labels_Z,True,data_z,None,"NL")   
    
    #print kernel_varBinUp.Integral()
    '''
    ############## make stackplots of all backgrounds ####################
    purity = sys.argv[1] 
    year = sys.argv[2] 
    f_Wjets  = ROOT.TFile("../../interactive/JJ_"+str(year)+"_WJets_"+str(purity)+".root","READ")
    f_Zjets  = ROOT.TFile("../../interactive/JJ_"+str(year)+"_ZJets_"+str(purity)+".root","READ")
    f_ttjets = ROOT.TFile("../../interactive/JJ_"+str(year)+"_TTJets_"+str(purity)+".root","READ")
    
    #qcd = f_sample.Get("nonRes")
    Wjets = f_Wjets.Get("WJets")
    Zjets = f_Zjets.Get("ZJets")
    ttbar = f_ttjets.Get("TTJets")
    #  qcd  .SetName("nonRes")
    Wjets.SetName("WJets_all")
    Zjets.SetName("ZJets_all")
    ttbar.SetName("ttbar_all")
    
    lumi = 59690.
    Wjets.Scale(lumi)
    Zjets.Scale(lumi)

    if year == 2018:
     lumi = 59690. #to be checked! https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis                                                                                                                                                                                   
    if year == 2017:
     lumi = 41367.
    elif year == 2016:
     lumi = 35900.

    #   qcd  .Scale(lumi)
    ttbar.Scale(lumi)
    
    #    print qcd
    print Wjets
    print Zjets
    print ttbar
    
    xrange = [["px0",55,215]]
    yrange = [["py0",55,215]]
    zrange = [["pz0",1126,5000]]

    backgrounds = [ttbar,Zjets,Wjets] #,Zjets,Wjets]#,qcd]

    plotStack(backgrounds,xrange,yrange,zrange,"x",year+"_"+purity)
    plotStack(backgrounds,xrange,yrange,zrange,"y",year+"_"+purity)
    plotStack(backgrounds,xrange,yrange,zrange,"z",year+"_"+purity)
    
    
    #irene print "qcd : "+str(qcd.Integral())+" Wjets "+str(Wjets.Integral()) + " Zjets "+str(Zjets.Integral()) + " ttbar "+str(ttbar.Integral())

    #irene print " Wjets "+str(Wjets.Integral()/qcd.Integral())+ " Zjets "+str(Zjets.Integral()/qcd.Integral())+" ttbar "+str(ttbar.Integral()/qcd.Integral())
    #irene print (Wjets.Integral()+Zjets.Integral()+ttbar.Integral())/qcd.Integral()

    print "Z/ W+ttbar "+str(+Zjets.Integral()/(Wjets.Integral()+ttbar.Integral()))

