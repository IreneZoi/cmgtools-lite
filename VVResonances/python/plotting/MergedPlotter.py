import ROOT
import sys
from array import array
import pickle
from CMGTools.VVResonances.python.plotting.PlotterBase import PlotterBase
class MergedPlotter(PlotterBase):

    def __init__(self,plotters):
        super(MergedPlotter,self).__init__()
        self.plotters=plotters
        self.corrFactors=plotters[0].corrFactors

    def applySmoothing(self):
        for plotter in self.plotters:
            plotter.applySmoothing()


    def scan(self,var,cut):
        for plotter in self.plotters:
            plotter.tree.SetScanField(0)
            plotter.tree.Scan(var,cut,"colsize=20")


    def drawTH1(self,var,cuts,lumi,bins,min,max,titlex = "",units = "",drawStyle = "HIST"):
        h=None
        for plotter in self.plotters:
            if h is None:
                h=plotter.drawTH1(var,cuts,lumi,bins,min,max,titlex,units,drawStyle)
            else:
                h.Add(plotter.drawTH1(var,cuts,lumi,bins,min,max,titlex,units,drawStyle))
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        h.GetXaxis().SetTitle(titlex+ " ["+units+"]")
        return h

    def drawTH2(self,var,cuts,lumi,binsx,minx,maxx,binsy,miny,maxy,titlex = "",unitsx = "",titley = "",unitsy = "",drawStyle = "COLZ"):
        h=None
        for plotter in self.plotters:
            if h is None:
                h=plotter.drawTH2(var,cuts,lumi,binsx,minx,maxx,binsy,miny,maxy,titlex,unitsx,titley,unitsy,drawStyle)
            else:
                h.Add(plotter.drawTH2(var,cuts,lumi,binsx,minx,maxx,binsy,miny,maxy,titlex,unitsx,titley,unitsy,drawStyle))

#        h.SetLineStyle(self.linestyle)
#        h.SetLineColor(self.linecolor)
#        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        return h


    def drawProfile(self,var,cuts,lumi,binsx,minx,maxx,miny,maxy,titlex = "",unitsx = "",titley = "",unitsy = "",drawStyle = "COLZ"):
        h=None
        for plotter in self.plotters:
            if h is None:
                h=plotter.drawProfile(var,cuts,lumi,binsx,minx,maxx,miny,maxy,titlex,unitsx,titley,unitsy,drawStyle)
            else:
                h.Add(plotter.drawProfile(var,cuts,lumi,binsx,minx,maxx,miny,maxy,titlex,unitsx,titley,unitsy,drawStyle))

#        h.SetLineStyle(self.linestyle)
#        h.SetLineColor(self.linecolor)
#        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        return h


    def drawTH3(self,var,cuts,lumi,binsx,minx,maxx,binsy,miny,maxy,binsz,minz,maxz,titlex = "",unitsx = "",titley = "",unitsy = "",drawStyle = "COLZ"):
        h=None
        for plotter in self.plotters:
            if h is None:
                h=plotter.drawTH3(var,cuts,lumi,binsx,minx,maxx,binsy,miny,maxy,binsz,minz,maxz,titlex,unitsx,titley,unitsy,drawStyle)
            else:
                h.Add(plotter.drawTH3(var,cuts,lumi,binsx,minx,maxx,binsy,miny,maxy,binsz,minz,maxz,titlex,unitsx,titley,unitsy,drawStyle))

        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        return h


    def drawTH2Binned(self,var,cuts,lumi,binningx,binningy,titlex = "",unitsx = "",titley = "",unitsy = "",drawStyle = "COLZ"):
        h=None
        for plotter in self.plotters:
            if h is None:
                h=plotter.drawTH2Binned(var,cuts,lumi,binningx,binningy,titlex,unitsx,titley,unitsy,drawStyle)
            else:
                h.Add(plotter.drawTH2Binned(var,cuts,lumi,binningx,binningy,titlex,unitsx,titley,unitsy,drawStyle))

#        h.SetLineStyle(self.linestyle)
#        h.SetLineColor(self.linecolor)
#        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        return h


    def drawTH3Binned(self,var,cuts,lumi,binningx,binningy,binningz,titlex = "",unitsx = "",titley = "",unitsy = "",titlez="",unitsz="",drawStyle = "COLZ"):
        h=None
        for plotter in self.plotters:
            if h is None:
                h=plotter.drawTH3Binned(var,cuts,lumi,binningx,binningy,binningz,titlex,unitsx,titley,unitsy,titlez,unitsz,drawStyle)
            else:
                h.Add(plotter.drawTH3Binned(var,cuts,lumi,binningx,binningy,binningz,titlex,unitsx,titley,unitsy,titlez,unitsz,drawStyle))

#        h.SetLineStyle(self.linestyle)
#        h.SetLineColor(self.linecolor)
#        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        h.GetZaxis().SetTitle(titlez+ " ["+unitsz+"]")
        return h

    
    def drawTH1Binned(self,var,cuts,lumi,binningx,titlex = "",unitsx = "",drawStyle = "COLZ"):
        h=None
        for plotter in self.plotters:
            if h is None:
                h=plotter.drawTH1Binned(var,cuts,lumi,binningx,titlex,unitsx,drawStyle)
            else:
                h.Add(plotter.drawTH1Binned(var,cuts,lumi,binningx,titlex,unitsx,drawStyle))
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        return h

    def makeDataSet(self,var,cut,maxN):
        data=self.plotters[0].makeDataSet(var,cut,maxN)
        for i in range(1,len(self.plotters)):
            data.append(self.plotters[i].makeDataSet(var,cut,maxN))
        return data
    
