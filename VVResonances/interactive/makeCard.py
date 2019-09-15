from tools.DatacardTools import *
import sys,os
import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
from CMGTools.VVResonances.statistics.DataCardMaker import DataCardMaker
cmd='combineCards.py '

sf_qcd = 1.0
pseudodata = ""#"ZprimeZH"
outlabel = ""#"sigonly_ZprimeZH_M2000"

datasets=['2016']#,'2017']

resultsDir = {'2016':'results_2016','2017':'results_2017'}

lumi = {'2016':35900,'2017':41367}
lumi_unc = {'2016':1.025,'2017':1.023}

scales = {"2017" :[0.983,1.08], "2016":[1.014,1.086]}
scalesHiggs = {"2017" :[1.,1.], "2016":[1.,1.]}

#quick fix to add VH !!!
vtag_unc = {'VV_HPHP':{},'VV_HPLP':{},'VV_LPLP':{},'VH_HPHP':{},'VH_HPLP':{},'VH_LPHP':{}}
vtag_unc['VV_HPHP'] = {'2016':'1.232/0.792','2017':'1.269/0.763'}
vtag_unc['VV_HPLP'] = {'2016':'0.882/1.12','2017':'0.866/1.136'}    
vtag_unc['VV_LPLP'] = {'2016':'1.063','2017':'1.043'}
vtag_unc['VH_HPHP'] = {'2016':'1.232/0.792','2017':'1.269/0.763'}
vtag_unc['VH_HPLP'] = {'2016':'0.882/1.12','2017':'0.866/1.136'}    
vtag_unc['VH_LPHP'] = {'2016':'1.063','2017':'1.043'}

vtag_pt_dependence = {'VV_HPHP':'((1+0.06*log(MH/2/300))*(1+0.06*log(MH/2/300)))','VV_HPLP':'((1+0.06*log(MH/2/300))*(1+0.07*log(MH/2/300)))','VH_HPHP':'((1+0.06*log(MH/2/300))*(1+0.06*log(MH/2/300)))','VH_HPLP':'((1+0.06*log(MH/2/300))*(1+0.07*log(MH/2/300)))','VH_LPHP':'((1+0.06*log(MH/2/300))*(1+0.06*log(MH/2/300)))'}

'''
vtag_unc = {'VV_HPHP':{},'VV_HPLP':{},'VV_LPLP':{}}
vtag_unc['VV_HPHP'] = {'2016':'1.232/0.792','2017':'1.269/0.763'}
vtag_unc['VV_HPLP'] = {'2016':'0.882/1.12','2017':'0.866/1.136'}    
vtag_unc['VV_LPLP'] = {'2016':'1.063','2017':'1.043'}

vtag_pt_dependence = {'VV_HPHP':'((1+0.06*log(MH/2/300))*(1+0.06*log(MH/2/300)))','VV_HPLP':'((1+0.06*log(MH/2/300))*(1+0.07*log(MH/2/300)))'}
'''
  
#purities= ['VV_HPLP']
#purities= ['VV_HPLP','VV_HPHP']
#purities= ['VH_HPLP','VH_HPHP',VH_LPHP]
purities= ['VV_HPLP','VV_HPHP','VH_HPLP','VH_HPHP','VH_LPHP']
#signals = ["BulkGWW", "BulkGZZ","ZprimeWW","WprimeWZ","VprimeWV","'ZprimeZH'"]
signals = ["BulkGWW"]

Tools = DatacardTools(scales,scalesHiggs,vtag_pt_dependence,lumi_unc,vtag_unc,sf_qcd,pseudodata,outlabel)

for sig in signals:
  cmd ="combineCards.py"
  for dataset in datasets:
    cmd_combo="combineCards.py"
    cmd_combo_vv="combineCards.py"
    cmd_combo_vh="combineCards.py"

    for p in purities:

      ncontrib = 0
      
      cat='_'.join(['JJ',sig,p,'13TeV_'+dataset])
      card=DataCardMaker('',p,'13TeV_'+dataset,lumi[dataset],'JJ',cat)
      cmd=cmd+" "+cat.replace('_%s'%sig,'')+'=datacard_'+cat+'.txt '
      cmd_combo=cmd_combo+" "+cat.replace('_%s'%sig,'')+'=datacard_'+cat+'.txt '
      if p.find("VV") !=-1:
        cmd_combo_vv=cmd_combo_vv+" "+cat.replace('_%s'%sig,'')+'=datacard_'+cat+'.txt '
      if p.find("VH") !=-1:
        cmd_combo_vh=cmd_combo_vh+" "+cat.replace('_%s'%sig,'')+'=datacard_'+cat+'.txt '

      cardName='datacard_'+cat+'.txt'
      workspaceName='workspace_'+cat+outlabel+'.root'

      if sig.find("VV")!=-1 or sig.find("WV")!=-1 :            
        Tools.AddSignal(card,dataset,p,sig,resultsDir[dataset],ncontrib)
      else:
        Tools.AddSingleSignal(card,dataset,p,sig,resultsDir[dataset],ncontrib)                 
      ncontrib+=1

      '''
      #comment out W/Z jets as they need new modelling (and TTbar too) for VH
      rootFileMVV = resultsDir[dataset]+'/JJ_%s_WJets_MVV_'%dataset+p+'.root' 
      rootFileNorm = resultsDir[dataset]+'/JJ_%s_WJets_%s.root'%(dataset,p)
      Tools.AddWResBackground(card,dataset,p,rootFileMVV,rootFileNorm,resultsDir[dataset],ncontrib)
      ncontrib+=1
      
      rootFileMVV = resultsDir[dataset]+'/JJ_%s_ZJets_MVV_'%dataset+p+'.root'
      rootFileNorm = resultsDir[dataset]+"/JJ_%s_ZJets_%s.root"%(dataset,p)
      Tools.AddZResBackground(card,dataset,p,rootFileMVV,rootFileNorm,resultsDir[dataset],ncontrib)
      ncontrib+=1
      '''

#      rootFile3DPDF = resultsDir[dataset]+'/JJ_2016_nonRes_3D_VV_HPLP.root'
      rootFile3DPDF = resultsDir[dataset]+"/save_new_shapes_%s_pythia_"%dataset+p+"_3D.root"

#      rootFileNorm = resultsDir[dataset]+"/JJ_PDnoVjets_VV_HPLP.root"
      rootFileNorm = resultsDir[dataset]+"/JJ_%s_nonRes_"%dataset+p+".root"   
      print("AddNonResBackground")
      Tools.AddNonResBackground(card,dataset,p,rootFile3DPDF,rootFileNorm,ncontrib) 
      print(" Added NonResBackground")
#      rootFileData = resultsDir[dataset]+"/JJ_"+p+".root"
#      histName="data"
#      scaleData=1.0 #if you ru on real data
      #psseudodata
      rootFileData = resultsDir[dataset]+"/JJ_PDnoVjets_VV_HPLP.root"          
      histName="data" 
      scaleData=lumi[dataset]
#      rootFileData = "results_2016/JJ_2016_nonRes_VV_HPLP.root"   
#      histName="nonRes"
#      scaleData=lumi[dataset] #if you run on MC for transfer-kernel
      if pseudodata=="ZprimeZH":
       rootFileData = resultsDir[dataset]+"/JJ_ZprimeZH_VH_all_M2000.root"
       histName="data_obs"
       scaleData=1.0
      if pseudodata=="WprimeWZ":
       rootFileData = resultsDir[dataset]+"/JJ_WprimeWZ_VV_HPLP_M4500.root" 
       histName="data_obs"    
       scaleData=1.0
      print("Adding data")
      Tools.AddData(card,rootFileData,histName,scaleData)
      print("data added")

      Tools.AddSigSystematics(card,sig,dataset,p,1)
      Tools.AddResBackgroundSystematics(card,p)
      Tools.AddNonResBackgroundSystematics(card,p)
        
      card.makeCard()
      
      t2wcmd = "text2workspace.py %s -o %s"%(cardName,workspaceName)
      print(t2wcmd)
      os.system(t2wcmd)
    del card
    '''
    #make combined HPHP+HPLP card   
    combo_card = 'datacard_'+cat.replace("_HPHP","").replace("_HPLP","").replace("_LPLP","")+'.txt'
    combo_workspace = 'workspace_'+cat.replace("_HPHP","").replace("_HPLP","").replace("_LPLP","")+'.root'
    os.system('rm %s'%combo_card)
    cmd_combo+=' >> %s'%combo_card
    print cmd_combo
    os.system(cmd_combo)
    t2wcmd = "text2workspace.py %s -o %s"%(combo_card,combo_workspace)
    print t2wcmd
    os.system(t2wcmd)
    '''


    '''
    #make combined VV HPHP+HPLP card
    print "######### Make VV card ###########"
    combo_card = 'datacard_'+cat.replace("VV_HPHP","VV").replace("VV_HPLP","VV").replace("VV_LPLP","")+'.txt'
    print combo_card
    combo_workspace = 'workspace_'+cat.replace("VV_HPHP","VV").replace("VV_HPLP","VV").replace("VV_LPLP","")+'.root'
    print combo_workspace
    os.system('rm %s'%combo_card)
    cmd_combo+=' >> %s'%combo_card
    print cmd_combo
    os.system(cmd_combo)
    t2wcmd = "text2workspace.py %s -o %s"%(combo_card,combo_workspace)
    print t2wcmd
    os.system(t2wcmd)
    #make combined VH HPHP+HPLP+LPHP card
    print "######### Make VH card ###########"
    combo_card = 'datacard_'+cat.replace("VH_HPHP","VH").replace("VH_HPLP","VH").replace("VH_LPHP","")+'.txt'
    print combo_card
    combo_workspace = 'workspace_'+cat.replace("VH_HPHP","VH").replace("VH_HPLP","VH").replace("VH_LPHP","")+'.root'
    print combo_workspace
    os.system('rm %s'%combo_card)
    cmd_combo+=' >> %s'%combo_card
    print cmd_combo
    os.system(cmd_combo)
    t2wcmd = "text2workspace.py %s -o %s"%(combo_card,combo_workspace)
    print t2wcmd
    os.system(t2wcmd)
    '''

    '''
    print "######### Make VV card ###########"   
    combo_card_vv="datacard_JJ_BulkGWW_VV_13TeV_2016.txt"
    combo_workspace_vv = "workspace_JJ_BulkGWW_VV_13TeV_2016.root"
    os.system('rm %s'%combo_card_vv)
    cmd_combo_vv+=' >> %s'%combo_card_vv
    print cmd_combo_vv
    os.system(cmd_combo_vv)
    t2wcmd = "text2workspace.py %s -o %s"%(combo_card_vv,combo_workspace_vv)
    print t2wcmd
    os.system(t2wcmd)

    print "######### Make VH card ###########"   
    combo_card_vh="datacard_JJ_BulkGWW_VH_13TeV_2016.txt"
    combo_workspace_vh = "workspace_JJ_BulkGWW_VH_13TeV_2016.root"
    os.system('rm %s'%combo_card_vh)
    cmd_combo_vh+=' >> %s'%combo_card_vh
    print cmd_combo_vh
    os.system(cmd_combo_vh)
    t2wcmd = "text2workspace.py %s -o %s"%(combo_card_vh,combo_workspace_vh)
    print t2wcmd
    os.system(t2wcmd)

    print "######### Make VV + VH card ###########"     
    combo_card_tot="datacard_JJ_BulkGWW_TOT_13TeV_2016.txt"
    combo_workspace_tot = "workspace_JJ_BulkGWW_TOT_13TeV_2016.root"
    os.system('rm %s'%combo_card_tot)
#    cmd_combo_tot+=' >> %s'%combo_card_tot
    cmd_combo_tot="combineCards.py JJ_VV_13TeV_2016=datacard_JJ_BulkGWW_VV_13TeV_2016.txt JJ_VH_13TeV_2016=datacard_JJ_BulkGWW_VH_13TeV_2016.txt >> datacard_JJ_BulkGWW_TOT_13TeV_2016.txt"
    print cmd_combo_tot
    os.system(cmd_combo_tot)
    t2wcmd = "text2workspace.py %s -o %s"%(combo_card_tot,combo_workspace_tot)
    print t2wcmd
    os.system(t2wcmd)

    '''
    #make combined VV+VH card
    print("######### Make VV + VH card ###########")
    print(cat)
    combo_card = 'datacard_'+cat.replace("VV_HPHP","").replace("VV_HPLP","").replace("VV_LPLP","").replace("VH_HPHP","").replace("VH_HPLP","").replace("VH_LPHP","")+'.txt'
    print(combo_card)
    combo_workspace = 'workspace_'+cat.replace("VV_HPHP","").replace("VV_HPLP","").replace("VV_LPLP","").replace("VH_HPHP","").replace("VH_HPLP","").replace("VH_LPHP","")+'.root'
    print(combo_workspace)
    os.system('rm %s'%combo_card)
    cmd_combo+=' >> %s'%combo_card
    print(cmd_combo)
    os.system(cmd_combo)
    t2wcmd = "text2workspace.py %s -o %s"%(combo_card,combo_workspace)
    print(t2wcmd)
    os.system(t2wcmd)

  #make combine 2016+2017 card
  #combo_card = 'datacard_'+cat.replace("_HPHP","").replace("_HPLP","").replace("_LPLP","").replace('_2016','').replace('_2017','')+'.txt'
  #combo_workspace = 'workspace_'+cat.replace("_HPHP","").replace("_HPLP","").replace("_LPLP","").replace('_2016','').replace('_2017','')+'.root'
  #os.system('rm %s'%combo_card)
  #cmd+=' >> %s'%combo_card
  #print cmd

  
  
  #os.system(cmd)
  #t2wcmd = "text2workspace.py %s -o %s"%(combo_card,combo_workspace)
  #print t2wcmd
  #os.system(t2wcmd)


