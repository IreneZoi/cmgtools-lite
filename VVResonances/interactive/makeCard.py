import sys,os
import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
from CMGTools.VVResonances.statistics.DataCardMaker import DataCardMaker
cmd='combineCards.py '


datasets=['2016']#,'2017']
addTT = False

pseudodata = "ZprimeZH"
outlabel="_sigonly_M2000"

lumi = {'2016':35900,'2017':41367}
lumi_unc = {'2016':1.025,'2017':1.023}

scales = {"2017" :[0.983,1.08], "2016":[1.014,1.086]}

vtag_unc = {'HPHP':{},'HPLP':{},'LPLP':{}}
vtag_unc['HPHP'] = {'2016':'1.232/0.792','2017':'1.269/0.763'}
vtag_unc['HPLP'] = {'2016':'0.882/1.12','2017':'0.866/1.136'}    
vtag_unc['LPLP'] = {'2016':'1.063','2017':'1.043'}

vtag_pt_dependence = {'HPHP':'((1+0.06*log(MH/2/300))*(1+0.06*log(MH/2/300)))','HPLP':'((1+0.06*log(MH/2/300))*(1+0.07*log(MH/2/300)))'}
  
purities= ['HPLP'] #['HPHP','HPLP']
#signals = ["BulkGWW", "BulkGZZ","ZprimeWW","WprimeWZ","VprimeWV"]
#signals = ['VprimeWV']
#signals = ['BulkGVV']
#signals = ['BulkGZZ']
signals = ['ZprimeZH']

for sig in signals:
  cmd ="combineCards.py"
  for dataset in datasets:
    cmd_combo="combineCards.py"
    for p in purities:

      ncontrib = 0
      
      cat='_'.join(['JJ',sig,p,'13TeV_'+dataset])
      card=DataCardMaker('',p,'13TeV_'+dataset,lumi[dataset],'JJ',cat)
      cmd=cmd+" "+cat.replace('_%s'%sig,'')+'=datacard_'+cat+'.txt '
      cmd_combo=cmd_combo+" "+cat.replace('_%s'%sig,'')+'=datacard_'+cat+'.txt '
      cardName='datacard_'+cat+'.txt'
      workspaceName='workspace_'+cat+outlabel+'.root'
      
      #SIGNAL
      if sig=='VprimeWV':
       sig1 = 'ZprimeWW'
       card.addMVVSignalParametricShape("%s_MVV"%sig1,"MJJ",dataset+"/JJ_%s_MVV.json"%sig1,{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig1,"MJ1" ,dataset+"/JJ_%s_MJl1_"%sig1+p+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig1,"MJ2" ,dataset+"/JJ_%s_MJl2_"%sig1+p+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       card.addParametricYieldHVTBR("%s"%sig1,ncontrib-1,dataset+"/JJ_%s_"%sig1+p+"_yield.json","../scripts/theoryXsec/HVTB.json","CX0(pb)","BRWW",1000.,'CMS_tau21_PtDependence',vtag_pt_dependence[p],1.0)
       card.product3D("%s"%sig1,"%s_Wqq1"%sig1,"%s_Wqq2"%sig1,"%s_MVV"%sig1)
       sig2 = 'WprimeWZ'
       card.addMVVSignalParametricShape("%s_MVV"%sig2,"MJJ",dataset+"/JJ_%s_MVV.json"%sig2,{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig2,"MJ1" ,dataset+"/JJ_%s_MJl1_"%sig2+p+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig2,"MJ2" ,dataset+"/JJ_%s_MJl2_"%sig2+p+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       card.addParametricYieldHVTBR("%s"%sig2,ncontrib,dataset+"/JJ_%s_"%sig2+p+"_yield.json","../scripts/theoryXsec/HVTB.json","CX+(pb),CX-(pb)","BRWZ",1000.,'CMS_tau21_PtDependence',vtag_pt_dependence[p],1.0)
       card.product3D("%s"%sig2,"%s_Wqq1"%sig2,"%s_Wqq2"%sig2,"%s_MVV"%sig2)
      elif sig=='BulkGVV':
       sig1 = 'BulkGWW'
       card.addMVVSignalParametricShape("%s_MVV"%sig1,"MJJ",dataset+"/JJ_%s_MVV.json"%sig1,{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig1,"MJ1" ,dataset+"/JJ_%s_MJl1_"%sig1+p+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig1,"MJ2" ,dataset+"/JJ_%s_MJl2_"%sig1+p+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       card.addParametricYieldHVTBR("%s"%sig1,ncontrib-1,dataset+"/JJ_%s_"%sig1+p+"_yield.json","../scripts/theoryXsec/BulkG.json","sigma","BRWW",1000.,'CMS_tau21_PtDependence',vtag_pt_dependence[p],1.0)
       card.product3D("%s"%sig1,"%s_Wqq1"%sig1,"%s_Wqq2"%sig1,"%s_MVV"%sig1)
       sig2 = 'BulkGZZ'
       card.addMVVSignalParametricShape("%s_MVV"%sig2,"MJJ",dataset+"/JJ_%s_MVV.json"%sig2,{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig2,"MJ1" ,dataset+"/JJ_%s_MJl1_"%sig2+p+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig2,"MJ2" ,dataset+"/JJ_%s_MJl2_"%sig2+p+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       card.addParametricYieldHVTBR("%s"%sig2,ncontrib,dataset+"/JJ_%s_"%sig2+p+"_yield.json","../scripts/theoryXsec/BulkG.json","sigma","BRZZ",1000.,'CMS_tau21_PtDependence',vtag_pt_dependence[p],1.0)
       card.product3D("%s"%sig2,"%s_Wqq1"%sig2,"%s_Wqq2"%sig2,"%s_MVV"%sig2)
      elif sig.find("H")==-1:
       card.addMVVSignalParametricShape("%s_MVV"%sig,"MJJ",dataset+"/JJ_%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1"%sig,"MJ1" ,dataset+"/JJ_%s_%s_MJrandom_VV_"%(sig,dataset)+p+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2"%sig,"MJ2" ,dataset+"/JJ_%s_%s_MJrandom_VV_"%(sig,dataset)+p+".json",{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       card.addParametricYieldWithUncertainty("%s"%sig,ncontrib,dataset+"/JJ_%s_%s_VV_%s_yield.json"%(sig,dataset,p),1,'CMS_tau21_PtDependence',vtag_pt_dependence[p],1.0)             
       card.product3D("%s"%sig,"%s_Wqq1"%sig,"%s_Wqq2"%sig,"%s_MVV"%sig)
      elif sig.find("H")!=-1:
    # only for test purposes put VH_all only as signal category!!!!!
       card.addMVVSignalParametricShape("%s_MVV_c1"%sig,"MJJ",dataset+"/JJ_%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeHiggs("%s_Wqq1_c1"%sig,"MJ1" ,dataset+"/JJ_Hjet_%s_%s_MJrandom_VH_all.json"%(sig,dataset),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq2_c1"%sig,"MJ2" ,dataset+"/JJ_Vjet_%s_%s_MJrandom_VH_all.json"%(sig,dataset),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       #card.product3D("%s_c1"%sig,"%s_Wqq1_c1"%sig,"%s_Wqq2_c1"%sig,"%s_MVV_c1"%sig) 
       card.conditionalProduct2("%s_c1"%sig,"%s_Wqq1_c1"%sig,"%s_Wqq2_c1"%sig,"%s_MVV_c1"%sig,"{MJ1,MJ2}")
       
       card.addMVVSignalParametricShape("%s_MVV_c2"%sig,"MJJ",dataset+"/JJ_%s_%s_MVV.json"%(sig,dataset),{'CMS_scale_j':1},{'CMS_res_j':1.0})
       card.addMJJSignalParametricShapeNOEXP("%s_Wqq1_c2"%sig,"MJ1" ,dataset+"/JJ_Vjet_%s_%s_MJrandom_VH_all.json"%(sig,dataset),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       card.addMJJSignalParametricShapeHiggs("%s_Wqq2_c2"%sig,"MJ2" ,dataset+"/JJ_Hjet_%s_%s_MJrandom_VH_all.json"%(sig,dataset),{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
       #card.product3D("%s_c2"%sig,"%s_Wqq1_c2"%sig,"%s_Wqq2_c2"%sig,"%s_MVV_c2"%sig) 
       card.conditionalProduct2("%s_c2"%sig,"%s_Wqq1_c2"%sig,"%s_Wqq2_c2"%sig,"%s_MVV_c2"%sig,"{MJ1,MJ2}")
       
       card.sum("%s"%sig,"%s_c1"%sig,"%s_c2"%sig,"0.5")
       
       if pseudodata=="":
          card.addParametricYieldWithUncertainty("%s"%sig,ncontrib,dataset+"/JJ_%s_%s_VV_%s_yield.json"%(sig,dataset,p),1,'CMS_tau21_PtDependence',vtag_pt_dependence[p],1.0)
       elif outlabel.find("sigonly"):
           card.addParametricYieldWithUncertainty("%s"%sig,ncontrib,dataset+"/JJ_%s_%s_VV_%s_yield.json"%(sig,dataset,p),1,'CMS_tau21_PtDependence',vtag_pt_dependence[p],100.)
          
      ncontrib+=1

      #---------------------------------------------------------------------------------
      #Vjets
      sys.path.append(dataset)
      from JJ_WJets_HPLP_gauss import Wjets_TTbar_nonRes_l1, Wjets_TTbar_Res_l1, Wjets_TTbar_nonRes_l2, Wjets_TTbar_Res_l2
      from JJ_WJets_HPLP_gauss import Zjets_Res_l1, Zjets_Res_l2, Zjets_nonRes_l1, Zjets_nonRes_l2  
         
            # begin W+jets background :
      
      # W+jets 
      rootFile = dataset+'/JJ_WJets_MVV_%s_new.root'%p
      card.addHistoShapeFromFile("Wjets_mjj_c1",["MJJ"],rootFile,"histo_nominal",['PT:CMS_VV_JJ_Wjets_PTZ_'+p,'OPT:CMS_VV_JJ_Wjets_OPTZ_'+p],False,0)
      card.addMJJSignalShapeNOEXP("Wjets_mjetRes_l1","MJ1","",Wjets_TTbar_Res_l1,{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
      card.addGaussianShape("Wjets_mjetNonRes_l2","MJ2",Wjets_TTbar_nonRes_l2)
      card.product3D("Wjets_c1","Wjets_mjetRes_l1","Wjets_mjetNonRes_l2","Wjets_mjj_c1")
      
      # jets + W
      rootFile = dataset+'/JJ_WJets_MVV_%s_new.root'%p
      card.addHistoShapeFromFile("Wjets_mjj_c2",["MJJ"],rootFile,"histo_nominal",['PT:CMS_VV_JJ_Wjets_PTZ_'+p,'OPT:CMS_VV_JJ_Wjets_OPTZ_'+p],False,0)
      card.addMJJSignalShapeNOEXP("Wjets_mjetRes_l2","MJ2","",Wjets_TTbar_Res_l2,{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
      card.addGaussianShape("Wjets_mjetNonRes_l1","MJ1",Wjets_TTbar_nonRes_l1)
      card.product3D("Wjets_c2","Wjets_mjetRes_l2","Wjets_mjetNonRes_l1","Wjets_mjj_c2")
      card.sumPdf('Wjets',"Wjets_c1","Wjets_c2","CMS_ratio_Wjets_"+p)
     
      if pseudodata=="":
          card.addFixedYieldFromFile('Wjets',ncontrib,dataset+"/JJ_WJets_%s.root"%p,"WJets")
      elif outlabel.find("sigonly"):
          card.addFixedYieldFromFile('Wjets',ncontrib,dataset+"/JJ_WJets_%s.root"%p,"WJets",0.000001)
      ncontrib+=1
            
      # begin Z+jets background :
      
      # Z+jets 
      rootFile = dataset+'/JJ_ZJets_MVV_%s_new.root'%p
      card.addHistoShapeFromFile("Zjets_mjj_c1",["MJJ"],rootFile,"histo_nominal",['PT:CMS_VV_JJ_Zjets_PTZ_'+p,'OPT:CMS_VV_JJ_Zjets_OPTZ_'+p],False,0)
      card.addMJJSignalShapeNOEXP("Zjets_mjetRes_l1","MJ1","",Zjets_Res_l1,{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
      card.addGaussianShape("Zjets_mjetNonRes_l2","MJ2",Zjets_nonRes_l2)
      card.product3D("Zjets_c1","Zjets_mjetRes_l1","Zjets_mjetNonRes_l2","Zjets_mjj_c1")
      
      
      # jets + Z
      rootFile = dataset+'/JJ_ZJets_MVV_%s_new.root'%p
      card.addHistoShapeFromFile("Zjets_mjj_c2",["MJJ"],rootFile,"histo_nominal",['PT:CMS_VV_JJ_Zjets_PTZ_'+p,'OPT:CMS_VV_JJ_Zjets_OPTZ_'+p],False,0)
      card.addMJJSignalShapeNOEXP("Zjets_mjetRes_l2","MJ2","",Zjets_Res_l2,{'CMS_scale_prunedj':1.},{'CMS_res_prunedj':1.},scales[dataset])
      card.addGaussianShape("Zjets_mjetNonRes_l1","MJ1",Zjets_nonRes_l1)
      card.product3D("Zjets_c2","Zjets_mjetRes_l2","Zjets_mjetNonRes_l1","Zjets_mjj_c2")
      card.sumPdf('Zjets',"Zjets_c1","Zjets_c2","CMS_ratio_Zjets_"+p)
      
      if pseudodata=="":
            card.addFixedYieldFromFile('Zjets',ncontrib,dataset+"/JJ_ZJets_%s.root"%p,"ZJets") 
      elif outlabel.find("sigonly"):
          card.addFixedYieldFromFile('Zjets',ncontrib,dataset+"/JJ_ZJets_%s.root"%p,"ZJets",0.000001)
      

      ncontrib+=1
                  
      #QCD      
      rootFile=dataset+"/save_new_shapes_pythia_"+p+"_3D.root"   
      card.addHistoShapeFromFile("nonRes",["MJ1","MJ2","MJJ"],rootFile,"histo",['PT:CMS_VV_JJ_nonRes_PT_'+p,'OPT:CMS_VV_JJ_nonRes_OPT_'+p,'OPT3:CMS_VV_JJ_nonRes_OPT3_'+p,'altshape:CMS_VV_JJ_nonRes_altshape_'+p,'altshape2:CMS_VV_JJ_nonRes_altshape2_'+p],False,0) ,    
      
      
      if pseudodata=="":
          card.addFixedYieldFromFile("nonRes",ncontrib,dataset+"/JJ_nonRes_"+p+".root","nonRes",0.8)
      if outlabel.find("sigonly")!=-1:
          card.addFixedYieldFromFile("nonRes",ncontrib,dataset+"/JJ_nonRes_"+p+".root","nonRes",0.0000000000001)

      #DATA
      if pseudodata=="":
        card.importBinnedData(dataset+"/JJ_"+p+".root","data",["MJ1","MJ2","MJJ"]) 
      if pseudodata=="ZprimeZH":
          card.importBinnedData(dataset+"/JJ_ZprimeZH_VH_all_M2000.root","data_obs",["MJ1","MJ2","MJJ"])
      if pseudodata=="WprimeWZ":
          card.importBinnedData(dataset+"/JJ_WprimeWZ_VV_HPLP_M2000.root","data_obs",["MJ1","MJ2","MJJ"])
      
      #SYSTEMATICS
      #luminosity
      card.addSystematic("CMS_lumi","lnN",{'%s'%sig:lumi_unc[dataset],"Wjets":lumi_unc[dataset],"Zjets":lumi_unc[dataset]})

      #PDF uncertainty for the signal
      card.addSystematic("CMS_pdf","lnN",{'%s'%sig:1.01})
    

      #background normalization
      card.addSystematic("CMS_VV_JJ_nonRes_norm","lnN",{'nonRes':1.5})
      card.addSystematic("CMS_VV_JJ_Wjets_norm","lnN",{'Wjets':1.2})
      card.addSystematic("CMS_VV_JJ_Zjets_norm","lnN",{'Zjets':1.2})
      
     
        
      #tau21 
      card.addSystematic("CMS_VV_JJ_tau21_eff","lnN",{'%s'%sig:vtag_unc[p][dataset],"Wjets":vtag_unc[p][dataset],"Zjets":vtag_unc[p][dataset]})
             
      #pruned mass scale  
      card.addSystematic("CMS_scale_prunedj","param",[0.0,0.02])
      card.addSystematic("CMS_res_prunedj","param",[0.0,0.08])
      card.addSystematic("CMS_scale_j","param",[0.0,0.012])
      card.addSystematic("CMS_res_j","param",[0.0,0.08])
    
      #systematics for dijet part of V+jets background
      card.addSystematic("CMS_VV_JJ_Wjets_PTZ_"+p,"param",[0,0.1]) #0.333
      card.addSystematic("CMS_VV_JJ_Wjets_OPTZ_"+p,"param",[0,0.1]) #0.333
      card.addSystematic("CMS_VV_JJ_Zjets_PTZ_"+p,"param",[0,0.1]) #0.333
      card.addSystematic("CMS_VV_JJ_Zjets_OPTZ_"+p,"param",[0,0.1]) #0.333
      
    
      #alternative shapes for QCD background
      card.addSystematic("CMS_VV_JJ_nonRes_PT_"+p,"param",[0.0,0.333])
      card.addSystematic("CMS_VV_JJ_nonRes_OPT_"+p,"param",[0.0,0.333])
      card.addSystematic('CMS_VV_JJ_nonRes_altshape2_'+p,"param",[0.0,0.333])  
      card.addSystematic('CMS_VV_JJ_nonRes_altshape_'+p,"param",[0.0,0.333])
      card.addSystematic("CMS_VV_JJ_nonRes_OPT3_"+p,"param",[1.0,0.333])
        
      card.makeCard()
      
      t2wcmd = "text2workspace.py %s -o %s"%(cardName,workspaceName)
      print t2wcmd
      os.system(t2wcmd)
    del card
    #make combined HPHP+HPLP card   
    combo_card = 'datacard_'+cat.replace("_HPHP","").replace("_HPLP","").replace("_LPLP","")+'.txt'
    combo_workspace = 'workspace_'+cat.replace("_HPHP","").replace("_HPLP","").replace("_LPLP","")+pseudodata+'.root'
    os.system('rm %s'%combo_card)
    cmd_combo+=' >> %s'%combo_card
    print cmd_combo
    os.system(cmd_combo)
    t2wcmd = "text2workspace.py %s -o %s"%(combo_card,combo_workspace)
    print t2wcmd
    os.system(t2wcmd)
  
  #make combine 2016+2017 card
  combo_card = 'datacard_'+cat.replace("_HPHP","").replace("_HPLP","").replace("_LPLP","").replace('_2016','').replace('_2017','')+'.txt'
  combo_workspace = 'workspace_'+cat.replace("_HPHP","").replace("_HPLP","").replace("_LPLP","").replace('_2016','').replace('_2017','')+pseudodata+'.root'
  os.system('rm %s'%combo_card)
  cmd+=' >> %s'%combo_card
  print cmd

  
  
  os.system(cmd)
  t2wcmd = "text2workspace.py %s -o %s"%(combo_card,combo_workspace)
  print t2wcmd
  os.system(t2wcmd)
  
  
  




