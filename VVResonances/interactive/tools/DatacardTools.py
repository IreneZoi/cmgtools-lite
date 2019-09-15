#!/usr/bin/env python
import sys, os

class DatacardTools():

 def __init__(self,scales,scalesHiggs,vtag_pt_dependence,lumi_unc,vtag_unc,sfQCD,pseudodata,outlabel):
  
  self.scales=scales
  self.vtag_pt_dependence=vtag_pt_dependence
  self.lumi_unc = lumi_unc
  self.vtag_unc = vtag_unc
  self.sfQCD = sfQCD
  self.pseudodata = pseudodata
  self.outlabel = outlabel
  self.scalesHiggs=scalesHiggs
 
 

 def AddSignal(self,card,dataset,category,sig,resultsDir,ncontrib):

      if sig=='VprimeWV':
       sig1 = 'ZprimeWW'
       card.addMVVSignalParametricShape("%s_MVV"%sig1,"MJJ",resultsDir+"/JJ_%s_%s_MVV.json"%(sig1,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig1,"MJ1" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig1,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig1,"MJ2" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig1,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addParametricYieldHVTBR("%s"%sig1,ncontrib-1,resultsDir+"/JJ_%s_%s_"%(sig1,dataset)+category+"_yield.json","../scripts/theoryXsec/HVTB.json","CX0(pb)","BRWW",1000.,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)
       card.product3D("%s"%sig1,"%s_Wqq1"%sig1,"%s_Wqq2"%sig1,"%s_MVV"%sig1)
       sig2 = 'WprimeWZ'
       card.addMVVSignalParametricShape("%s_MVV"%sig2,"MJJ",resultsDir+"/JJ_%s_%s_MVV.json"%(sig2,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig2,"MJ1" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig2,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig2,"MJ2" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig2,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addParametricYieldHVTBR("%s"%sig2,ncontrib,resultsDir+"/JJ_%s_%s_"%(sig2,dataset)+category+"_yield.json","../scripts/theoryXsec/HVTB.json","CX+(pb),CX-(pb)","BRWZ",1000.,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)
       card.product3D("%s"%sig2,"%s_Wqq1"%sig2,"%s_Wqq2"%sig2,"%s_MVV"%sig2)
      if sig=='BulkGVV':
       sig1 = 'BulkGWW'
       card.addMVVSignalParametricShape("%s_MVV"%sig1,"MJJ",resultsDir+"/JJ_%s_%s_MVV.json"%(sig1,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig1,"MJ1" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig1,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig1,"MJ2" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig1,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addParametricYieldHVTBR("%s"%sig1,ncontrib-1,resultsDir+"/JJ_%s_%s_"%(sig1,dataset)+category+"_yield.json","../scripts/theoryXsec/BulkG.json","sigma","BRWW",1000.,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)
       card.product3D("%s"%sig1,"%s_Wqq1"%sig1,"%s_Wqq2"%sig1,"%s_MVV"%sig1)
       sig2 = 'BulkGZZ'
       card.addMVVSignalParametricShape("%s_MVV"%sig2,"MJJ",resultsDir+"/JJ_%s_%s_MVV.json"%(sig2,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig2,"MJ1" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig2,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig2,"MJ2" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig2,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addParametricYieldHVTBR("%s"%sig2,ncontrib,resultsDir+"/JJ_%s_%s_"%(sig2,dataset)+category+"_yield.json","../scripts/theoryXsec/BulkG.json","sigma","BRZZ",1000.,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)
       card.product3D("%s"%sig2,"%s_Wqq1"%sig2,"%s_Wqq2"%sig2,"%s_MVV"%sig2)
      if sig.find("H")==-1 and sig.find("WZ")==-1:
       card.addMVVSignalParametricShape("%s_MVV"%sig,"MJJ",resultsDir+"/JJ_%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig,"MJ1" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig,"MJ2" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.product3D("%s"%sig,"%s_Wqq1"%sig,"%s_Wqq2"%sig,"%s_MVV"%sig)
       card.addParametricYieldWithUncertainty("%s"%sig,ncontrib,resultsDir+"/JJ_%s_%s_"%(sig,dataset)+category+"_yield.json",1,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)             
       
      if sig.find("H")!=-1:
       card.addMVVSignalParametricShape("%s_MVV_c1"%sig,"MJJ",resultsDir+"/JJ_j1%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeHiggs("%s_Wqq1_c1"%sig,"MJ1" ,resultsDir+"/JJ_Hjet_%s_%s_MJrandom_%s.json"%(sig,dataset,category),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scalesHiggs[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2_c1"%sig,"MJ2" ,resultsDir+"/JJ_Vjet_%s_%s_MJrandom_%s.json"%(sig,dataset,category),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.conditionalProduct2("%s_c1"%sig,"%s_Wqq1_c1"%sig,"%s_Wqq2_c1"%sig,"%s_MVV_c1"%sig,"{MJ1,MJ2}")
       
       card.addMVVSignalParametricShape("%s_MVV_c2"%sig,"MJJ",resultsDir+"/JJ_j2%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1_c2"%sig,"MJ1" ,resultsDir+"/JJ_Vjet_%s_%s_MJrandom_%s.json"%(sig,dataset,category),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addMJJSignalParametricShapeHiggs("%s_Wqq2_c2"%sig,"MJ2" ,resultsDir+"/JJ_Hjet_%s_%s_MJrandom_%s.json"%(sig,dataset,category),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scalesHiggs[dataset])
       card.conditionalProduct2("%s_c2"%sig,"%s_Wqq1_c2"%sig,"%s_Wqq2_c2"%sig,"%s_MVV_c2"%sig,"{MJ1,MJ2}")
       card.sumSimple("%s"%sig,"%s_c1"%sig,"%s_c2"%sig,"0.5")
       card.addParametricYieldWithUncertainty("%s"%sig,ncontrib,resultsDir+"/JJ_%s_%s_%s_yield.json"%(sig,dataset,category),1,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)
       
       if self.pseudodata=="":
          card.addParametricYieldWithUncertainty("%s"%sig,ncontrib,resultsDir+"/JJ_%s_%s_%s_yield.json"%(sig,dataset,category),1,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)
       elif self.outlabel.find("sigonly"):
           card.addParametricYieldWithUncertainty("%s"%sig,ncontrib,resultsDir+"/JJ_%s_%s_%s_yield.json"%(sig,dataset,category),1,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],500.)
       
      if sig.find("WZ")!=-1:
       card.addMVVSignalParametricShape("%s_MVV_c1"%sig,"MJJ",resultsDir+"/JJ_j1%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1_c1"%sig,"MJ1" ,resultsDir+"/JJ_%s_%s_MJrandom_%s.json"%(sig,dataset,category),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2_c1"%sig,"MJ2" ,resultsDir+"/JJ_%s_%s_MJrandom_%s.json"%(sig,dataset,category),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.conditionalProduct2("%s_c1"%sig,"%s_Wqq1_c1"%sig,"%s_Wqq2_c1"%sig,"%s_MVV_c1"%sig,"{MJ1,MJ2}")
       
       card.addMVVSignalParametricShape("%s_MVV_c2"%sig,"MJJ",resultsDir+"/JJ_j2%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1_c2"%sig,"MJ1" ,resultsDir+"/JJ_%s_%s_MJrandom_%s.json"%(sig,dataset,category),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2_c2"%sig,"MJ2" ,resultsDir+"/JJ_%s_%s_MJrandom_%s.json"%(sig,dataset,category),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.conditionalProduct2("%s_c2"%sig,"%s_Wqq1_c2"%sig,"%s_Wqq2_c2"%sig,"%s_MVV_c2"%sig,"{MJ1,MJ2}")
       card.sumSimple("%s"%sig,"%s_c1"%sig,"%s_c2"%sig,"0.5")
       
       if self.pseudodata=="":
          card.addParametricYieldWithUncertainty("%s"%sig,ncontrib,resultsDir+"/JJ_%s_%s_%s_yield.json"%(sig,dataset,category),1,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)
       elif self.outlabel.find("sigonly")!=-1:
           card.addParametricYieldWithUncertainty("%s"%sig,ncontrib,resultsDir+"/JJ_%s_%s_%s_yield.json"%(sig,dataset,category),1,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],500.)

 def AddSingleSignal(self,card,dataset,category,sig,resultsDir,ncontrib):
  #NB since the signals are scaled by xsec BulkG = BulkGWW + BulkGZZ, they need to be rescaled accordingly when producing the limit plot!!! 
  if sig=='ZprimeWW':
   card.addMVVSignalParametricShape("%s_MVV"%sig,"MJJ",resultsDir+"/JJ_%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
   card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig,"MJ1" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig,"MJ2" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   card.addParametricYieldHVTBR("%s"%sig,ncontrib-1,resultsDir+"/JJ_%s_%s_"%(sig,dataset)+category+"_yield.json","../scripts/theoryXsec/HVTB.json","CX0(pb)","BRWW",1000.,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)
   card.product3D("%s"%sig,"%s_Wqq1"%sig,"%s_Wqq2"%sig,"%s_MVV"%sig)
  elif sig=='WprimeWZ':
   card.addMVVSignalParametricShape("%s_MVV"%sig,"MJJ",resultsDir+"/JJ_%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
   card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig,"MJ1" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig,"MJ2" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   card.addParametricYieldHVTBR("%s"%sig,ncontrib,resultsDir+"/JJ_%s_%s_"%(sig,dataset)+category+"_yield.json","../scripts/theoryXsec/HVTB.json","CX+(pb),CX-(pb)","BRWZ",1000.,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)
   card.product3D("%s"%sig,"%s_Wqq1"%sig,"%s_Wqq2"%sig,"%s_MVV"%sig)
  elif sig=='BulkGWW':
   card.addMVVSignalParametricShape("%s_MVV"%sig,"MJJ",resultsDir+"/JJ_%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
   card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig,"MJ1" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig,"MJ2" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   card.addParametricYieldHVTBR("%s"%sig,ncontrib-1,resultsDir+"/JJ_%s_%s_"%(sig,dataset)+category+"_yield.json","../scripts/theoryXsec/BulkG.json","sigma","BRWW",1000.,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)
   card.product3D("%s"%sig,"%s_Wqq1"%sig,"%s_Wqq2"%sig,"%s_MVV"%sig)
  elif sig=='BulkGZZ':
   card.addMVVSignalParametricShape("%s_MVV"%sig,"MJJ",resultsDir+"/JJ_%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
   card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig,"MJ1" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig,"MJ2" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   card.addParametricYieldHVTBR("%s"%sig,ncontrib,resultsDir+"/JJ_%s_%s_"%(sig,dataset)+category+"_yield.json","../scripts/theoryXsec/BulkG.json","sigma","BRZZ",1000.,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)
   card.product3D("%s"%sig,"%s_Wqq1"%sig,"%s_Wqq2"%sig,"%s_MVV"%sig)
  elif sig.find("H")==-1:
   card.addMVVSignalParametricShape("%s_MVV"%sig,"MJJ",resultsDir+"/JJ_%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
   card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig,"MJ1" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig,"MJ2" ,resultsDir+"/JJ_%s_%s_MJrandom_"%(sig,dataset)+category+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   card.addParametricYieldWithUncertainty("%s"%sig,ncontrib,resultsDir+"/JJ_%s_%s_"%(sig,dataset)+category+"_yield.json",1,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)             
   card.product3D("%s"%sig,"%s_Wqq1"%sig,"%s_Wqq2"%sig,"%s_MVV"%sig)
  elif sig.find("H")!=-1:
   # only for test purposes put VH_all only as signal category!!!!!
   card.addMVVSignalParametricShape("%s_MVV_c1"%sig,"MJJ",resultsDir+"/JJ_%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
   card.addMJJSignalParametricShapeHiggs("%s_Wqq1_c1"%sig,"MJ1" ,resultsDir+"/JJ_Hjet_%s_%s_MJrandom_VH_all.json"%(sig,dataset),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   card.addMJJSignalParametricShapeNOEXP("%s_Wqq2_c1"%sig,"MJ2" ,resultsDir+"/JJ_Vjet_%s_%s_MJrandom_VH_all.json"%(sig,dataset),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   #card.product3D("%s_c1"%sig,"%s_Wqq1_c1"%sig,"%s_Wqq2_c1"%sig,"%s_MVV_c1"%sig) 
   card.conditionalProduct2("%s_c1"%sig,"%s_Wqq1_c1"%sig,"%s_Wqq2_c1"%sig,"%s_MVV_c1"%sig,"{MJ1,MJ2}")
    
   card.addMVVSignalParametricShape("%s_MVV_c2"%sig,"MJJ",resultsDir+"/JJ_%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
   card.addMJJSignalParametricShapeNOEXP("%s_Wqq1_c2"%sig,"MJ1" ,resultsDir+"/JJ_Vjet_%s_%s_MJrandom_VH_all.json"%(sig,dataset),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   card.addMJJSignalParametricShapeHiggs("%s_Wqq2_c2"%sig,"MJ2" ,resultsDir+"/JJ_Hjet_%s_%s_MJrandom_VH_all.json"%(sig,dataset),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
   #card.product3D("%s_c2"%sig,"%s_Wqq1_c2"%sig,"%s_Wqq2_c2"%sig,"%s_MVV_c2"%sig) 
   card.conditionalProduct2("%s_c2"%sig,"%s_Wqq1_c2"%sig,"%s_Wqq2_c2"%sig,"%s_MVV_c2"%sig,"{MJ1,MJ2}")
    
   card.sum("%s"%sig,"%s_c1"%sig,"%s_c2"%sig,"0.5")
   
   if self.pseudodata=="":
    card.addParametricYieldWithUncertainty("%s"%sig,ncontrib,resultsDir+"/JJ_%s_%s_VV_%s_yield.json"%(sig,dataset,p),1,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],1.0)
   elif self.outlabel.find("sigonly"):
    card.addParametricYieldWithUncertainty("%s"%sig,ncontrib,resultsDir+"/JJ_%s_%s_VV_%s_yield.json"%(sig,dataset,p),1,'CMS_tau21_PtDependence',self.vtag_pt_dependence[category],100.)
    
 def AddWResBackground(self,card,dataset,category,rootFileMVV,rootFileNorm,resultsDir,ncontrib):
  
       sys.path.append(resultsDir)
       from JJ_WJets_VV_HPLP import Wjets_TTbar_nonRes_l1, Wjets_TTbar_Res_l1, Wjets_TTbar_nonRes_l2, Wjets_TTbar_Res_l2
       from JJ_WJets_VV_HPLP import Zjets_Res_l1, Zjets_Res_l2, Zjets_nonRes_l1, Zjets_nonRes_l2  
       
       # W+jets 

       card.addHistoShapeFromFile("Wjets_mjj_c1",["MJJ"],rootFileMVV,"histo_nominal",['PT:CMS_VV_JJ_Wjets_PTZ_'+category,'OPT:CMS_VV_JJ_Wjets_OPTZ_'+category],False,0)
       card.addMJJSignalShapeNOEXP("Wjets_mjetRes_l1","MJ1","",Wjets_TTbar_Res_l1,{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addGaussianShape("Wjets_mjetNonRes_l2","MJ2",Wjets_TTbar_nonRes_l2)
       card.product3D("Wjets_c1","Wjets_mjetRes_l1","Wjets_mjetNonRes_l2","Wjets_mjj_c1")
      
       # jets + W
       card.addHistoShapeFromFile("Wjets_mjj_c2",["MJJ"],rootFileMVV,"histo_nominal",['PT:CMS_VV_JJ_Wjets_PTZ_'+category,'OPT:CMS_VV_JJ_Wjets_OPTZ_'+category],False,0)
       card.addMJJSignalShapeNOEXP("Wjets_mjetRes_l2","MJ2","",Wjets_TTbar_Res_l2,{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addGaussianShape("Wjets_mjetNonRes_l1","MJ1",Wjets_TTbar_nonRes_l1)
       card.product3D("Wjets_c2","Wjets_mjetRes_l2","Wjets_mjetNonRes_l1","Wjets_mjj_c2")
       card.sumPdf('Wjets',"Wjets_c1","Wjets_c2","CMS_ratio_Wjets_"+category)
     
       print("outlabel "+self.outlabel)
       if self.pseudodata=="":
           card.addFixedYieldFromFile('Wjets',ncontrib,rootFileNorm,"WJets")
       if self.outlabel.find("sigonly")!=-1:
           print("add small yield")
           card.addFixedYieldFromFile('Wjets',ncontrib,rootFileNorm,"WJets",0.000001)

 def AddZResBackground(self,card,dataset,category,rootFileMVV,rootFileNorm,resultsDir,ncontrib):  
       print("add Zres background")
       sys.path.append(resultsDir)
       from JJ_WJets_VV_HPLP import Wjets_TTbar_nonRes_l1, Wjets_TTbar_Res_l1, Wjets_TTbar_nonRes_l2, Wjets_TTbar_Res_l2
       from JJ_WJets_VV_HPLP import Zjets_Res_l1, Zjets_Res_l2, Zjets_nonRes_l1, Zjets_nonRes_l2  
            
       # Z+jets 
       card.addHistoShapeFromFile("Zjets_mjj_c1",["MJJ"],rootFileMVV,"histo_nominal",['PT:CMS_VV_JJ_Zjets_PTZ_'+category,'OPT:CMS_VV_JJ_Zjets_OPTZ_'+category],False,0)
       card.addMJJSignalShapeNOEXP("Zjets_mjetRes_l1","MJ1","",Zjets_Res_l1,{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addGaussianShape("Zjets_mjetNonRes_l2","MJ2",Zjets_nonRes_l2)
       card.product3D("Zjets_c1","Zjets_mjetRes_l1","Zjets_mjetNonRes_l2","Zjets_mjj_c1")
           
       # jets + Z
       card.addHistoShapeFromFile("Zjets_mjj_c2",["MJJ"],rootFileMVV,"histo_nominal",['PT:CMS_VV_JJ_Zjets_PTZ_'+category,'OPT:CMS_VV_JJ_Zjets_OPTZ_'+category],False,0)
       card.addMJJSignalShapeNOEXP("Zjets_mjetRes_l2","MJ2","",Zjets_Res_l2,{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},self.scales[dataset])
       card.addGaussianShape("Zjets_mjetNonRes_l1","MJ1",Zjets_nonRes_l1)
       card.product3D("Zjets_c2","Zjets_mjetRes_l2","Zjets_mjetNonRes_l1","Zjets_mjj_c2")
       card.sumPdf('Zjets',"Zjets_c1","Zjets_c2","CMS_ratio_Zjets_"+category)
      
       if self.pseudodata=="":
             card.addFixedYieldFromFile('Zjets',ncontrib,rootFileNorm,"ZJets") 
       if self.outlabel.find("sigonly")!=-1:
           card.addFixedYieldFromFile('Zjets',ncontrib,rootFileNorm,"ZJets",0.000001)
       print("stop Zres background")
   
 def AddNonResBackground(self,card,dataset,category,rootFile3DPDF,rootFileNorm,ncontrib):
      print("adding HistoShapeFromFile")
      card.addHistoShapeFromFile("nonRes",["MJ1","MJ2","MJJ"],rootFile3DPDF,"histo",['PT:CMS_VV_JJ_nonRes_PT_'+category,'OPT:CMS_VV_JJ_nonRes_OPT_'+category,'OPT3:CMS_VV_JJ_nonRes_OPT3_'+category,'altshape:CMS_VV_JJ_nonRes_altshape_'+category,'altshape2:CMS_VV_JJ_nonRes_altshape2_'+category],False,0) ,    
      print("added HistoShapeFromFile")    
      if self.pseudodata=="":
          card.addFixedYieldFromFile("nonRes",ncontrib,rootFileNorm,"nonRes",self.sfQCD)
      if self.outlabel.find("sigonly")!=-1:
          card.addFixedYieldFromFile("nonRes",ncontrib,rootFileNorm,"nonRes",0.0000000000001)
 
 def AddData(self,card,fileData,histoName,scaleData):
  
      card.importBinnedData(fileData,histoName,["MJ1","MJ2","MJJ"],'data_obs',scaleData)
  
 def AddSigSystematics(self,card,sig,dataset,category,correlate):

      card.addSystematic("CMS_scale_prunedj","param",[0.0,0.02])
      card.addSystematic("CMS_res_prunedj","param",[0.0,0.08])
      card.addSystematic("CMS_scale_j","param",[0.0,0.012])
      card.addSystematic("CMS_res_j","param",[0.0,0.08])

      card.addSystematic("CMS_pdf","lnN",{'%s'%sig:1.01})    
      if correlate:
       card.addSystematic("CMS_lumi","lnN",{'%s'%sig:self.lumi_unc[dataset],"Wjets":self.lumi_unc[dataset],"Zjets":self.lumi_unc[dataset]})   
#       if category.find("VH")==-1 :
       card.addSystematic("CMS_VV_JJ_tau21_eff","lnN",{'%s'%sig:self.vtag_unc[category][dataset],"Wjets":self.vtag_unc[category][dataset],"Zjets":self.vtag_unc[category][dataset]})
       #else :
        # card.addSystematic("CMS_VH_JJ_tau21_eff","lnN",{'%s'%sig:self.vtag_unc[category][dataset],"Wjets":self.vtag_unc[category][dataset],"Zjets":self.vtag_unc[category][dataset]})
      else: 
       card.addSystematic("CMS_lumi","lnN",{'%s'%sig:self.lumi_unc[dataset]})
       print("Adding CMS_VV_JJ_tau21_eff for category ",category)
       print("dataset ",dataset)
       print("self.vtag_unc[category][dataset] ",self.vtag_unc[category][dataset]) 
       card.addSystematic("CMS_VV_JJ_tau21_eff","lnN",{'%s'%sig:self.vtag_unc[category][dataset]})
       #else:
        #card.addSystematic("CMS_VH_JJ_tau21_eff","lnN",{'%s'%sig:self.vtag_unc[category][dataset]})    
              
 def AddResBackgroundSystematics(self,card,category):
 
       card.addSystematic("CMS_VV_JJ_Wjets_norm","lnN",{'Wjets':1.2})
       card.addSystematic("CMS_VV_JJ_Zjets_norm","lnN",{'Zjets':1.2})    
       card.addSystematic("CMS_VV_JJ_Wjets_PTZ_"+category,"param",[0,0.1]) #0.333
       card.addSystematic("CMS_VV_JJ_Wjets_OPTZ_"+category,"param",[0,0.1]) #0.333
       card.addSystematic("CMS_VV_JJ_Zjets_PTZ_"+category,"param",[0,0.1]) #0.333
       card.addSystematic("CMS_VV_JJ_Zjets_OPTZ_"+category,"param",[0,0.1]) #0.333

 def AddNonResBackgroundSystematics(self,card,category):

      card.addSystematic("CMS_VV_JJ_nonRes_norm","lnN",{'nonRes':1.5})
      
      card.addSystematic("CMS_VV_JJ_nonRes_PT_"+category,"param",[0.0,0.333])
      card.addSystematic("CMS_VV_JJ_nonRes_OPT_"+category,"param",[0.0,0.333])
      card.addSystematic('CMS_VV_JJ_nonRes_altshape2_'+category,"param",[0.0,0.333])  
      card.addSystematic('CMS_VV_JJ_nonRes_altshape_'+category,"param",[0.0,0.333])
      card.addSystematic("CMS_VV_JJ_nonRes_OPT3_"+category,"param",[1.0,0.333])
