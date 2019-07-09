from functions import *

period = 2016
samples= str(period)+"/" #for V+jets we use 2018 samples also for 2016 because the 2016 ones are buggy and they need to be processed before to add the NLO weights!

#samples= str(period)+"_new/"
sorting = 'random'
#sorting = 'btag'

submitToBatch = False #Set to true if you want to submit kernels + makeData to batch!
runParallel   = False #Set to true if you want to run all kernels in parallel! This will exit this script and you will have to run mergeKernelJobs when your jobs are done! TODO! Add waitForBatchJobs also here?
dijetBinning = True
useTriggerWeights = False

#scale factors to be updated!
HPSF = 0.937
LPSF = 1.006
if period == 2017:
    HPSF = 0.955
    LPSF = 1.003
    
addOption = ""
if useTriggerWeights: 
    addOption = "-t"
    
if dijetBinning:
    HCALbinsMVVSignal=" --binsMVV 1,3,6,10,16,23,31,40,50,61,74,88,103,119,137,156,176,197,220,244,270,296,325,354,386,419,453,489,526,565,606,649,693,740,788,838,890,944,1000,1058,1126,1181,1246,1313,1383,1455,1530,1607,1687,1770,1856,1945,2037,2132,2231,2332,2438,2546,2659,2775,2895,3019,3147,3279,3416,3558,3704,3854,4010,4171,4337,4509,4686,4869,5058,5253,5500,5663,5877,6099,6328,6564,6808"
    dijetbins = [1126,1181,1246,1313,1383,1455,1530,1607,1687,1770,1856,1945,2037,2132,2231,2332,2438,2546,2659,2775,2895,3019,3147,3279,3416,3558,3704,3854,4010,4171,4337,4509,4686,4869,5058,5253,5500] # ,7060,7320,7589]
    # dijetbins = [1126,1181,1246,1313,1383,1455,1530,1607,1687,1770,1856,1945,2037,2132,2231,2332,2438,2546,2659,2775,2895,3019,3147,3279,3416,3558,3704,3854,4010,4171,4337,4509,4686,4869,5058,5253,5455,5663,5877,6099,6328,6564,6808,7060,7320,7589]
    HCALbinsMVV  =" --binsMVV "
    HCALbinsMVV += ','.join(str(e) for e in dijetbins)
else:
    HCALbinsMVV=""
    HCALbinsMVVSignal=""

if period == 2018:
    lumi = 59690. #to be checked! https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis
if period == 2017:
    lumi = 41367.    
elif period == 2016:
    lumi = 35900.
    
catVtag = {}
catHtag = {}


# For retuned DDT tau 21, use this
catVtag['HP1'] = '(jj_l1_tau2/jj_l1_tau1+(0.080*TMath::Log((jj_l1_softDrop_mass*jj_l1_softDrop_mass)/jj_l1_pt)))<0.43'
catVtag['HP2'] = '(jj_l2_tau2/jj_l2_tau1+(0.080*TMath::Log((jj_l2_softDrop_mass*jj_l2_softDrop_mass)/jj_l2_pt)))<0.43'
catVtag['LP1'] = '(jj_l1_tau2/jj_l1_tau1+(0.080*TMath::Log((jj_l1_softDrop_mass*jj_l1_softDrop_mass)/jj_l1_pt)))>0.43&&(jj_l1_tau2/jj_l1_tau1+(0.080*TMath::Log((jj_l1_softDrop_mass*jj_l1_softDrop_mass)/jj_l1_pt)))<0.79'
catVtag['LP2'] = '(jj_l2_tau2/jj_l2_tau1+(0.080*TMath::Log((jj_l2_softDrop_mass*jj_l2_softDrop_mass)/jj_l2_pt)))>0.43&&(jj_l2_tau2/jj_l2_tau1+(0.080*TMath::Log((jj_l2_softDrop_mass*jj_l2_softDrop_mass)/jj_l2_pt)))<0.79'
catVtag['NP1'] = '(jj_l1_tau2/jj_l1_tau1+(0.080*TMath::Log((jj_l1_softDrop_mass*jj_l1_softDrop_mass)/jj_l1_pt)))>0.79'
catVtag['NP2'] = '(jj_l2_tau2/jj_l2_tau1+(0.080*TMath::Log((jj_l2_softDrop_mass*jj_l2_softDrop_mass)/jj_l2_pt)))>0.79'

catHtag['HP1'] = '(jj_l1_MassDecorrelatedDeepBoosted_ZHbbvsQCD>0.8)'
catHtag['HP2'] = '(jj_l2_MassDecorrelatedDeepBoosted_ZHbbvsQCD>0.8)'
catHtag['LP1'] = '(jj_l1_MassDecorrelatedDeepBoosted_ZHbbvsQCD>0.5&&jj_l1_MassDecorrelatedDeepBoosted_ZHbbvsQCD<0.8)'
catHtag['LP2'] = '(jj_l2_MassDecorrelatedDeepBoosted_ZHbbvsQCD>0.5&&jj_l2_MassDecorrelatedDeepBoosted_ZHbbvsQCD<0.8)'
catHtag['NP1'] = '(jj_l1_MassDecorrelatedDeepBoosted_ZHbbvsQCD<0.5)'
catHtag['NP2'] = '(jj_l2_MassDecorrelatedDeepBoosted_ZHbbvsQCD<0.5)'

cuts={}

cuts['common'] = '((HLT_JJ)*(run>500) + (run<500))*(passed_METfilters&&passed_PVfilter&&njj>0&&jj_LV_mass>700&&abs(jj_l1_eta-jj_l2_eta)<1.3&&jj_l1_softDrop_mass>0.&&jj_l2_softDrop_mass>0.)' #with rho

#signal regions
if sorting == 'random':
 print "Use random sorting!"
 cuts['VH_HPHP'] = '(' + '('+  '&&'.join([catVtag['HP1'],catHtag['HP2']]) + ')' + '||' + '(' + '&&'.join([catVtag['HP2'],catHtag['HP1']]) + ')' + ')'
 cuts['VH_HPLP'] = '(' + '('+  '&&'.join([catVtag['HP1'],catHtag['LP2']]) + ')' + '||' + '(' + '&&'.join([catVtag['HP2'],catHtag['LP1']]) + ')' + ')'
 cuts['VH_LPHP'] = '(' + '('+  '&&'.join([catVtag['LP1'],catHtag['HP2']]) + ')' + '||' + '(' + '&&'.join([catVtag['LP2'],catHtag['HP1']]) + ')' + ')'
 cuts['VH_LPLP'] = '(' + '('+  '&&'.join([catVtag['LP1'],catHtag['LP2']]) + ')' + '||' + '(' + '&&'.join([catVtag['LP2'],catHtag['LP1']]) + ')' + ')'
 cuts['VH_all'] =  '('+  '||'.join([cuts['VH_HPHP'],cuts['VH_HPLP'],cuts['VH_LPHP'],cuts['VH_LPLP']]) + ')'
 cuts['VV_HPHP'] = '(' + '!' + cuts['VH_all'] + '&&' + '(' + '&&'.join([catVtag['HP1'],catVtag['HP2']]) + ')' + ')'
 #cuts['VV_HPLP'] = '(' + '!' + cuts['VH_all'] + '&&' + '(' + '('+  '&&'.join([catVtag['HP1'],catVtag['LP2']]) + ')' + '||' + '(' + '&&'.join([catVtag['HP2'],catVtag['LP1']]) + ')' + ')' + ')'
 cuts['VV_HPLP'] = '(' + '('+  '&&'.join([catVtag['HP1'],catVtag['LP2']]) + ')' + '||' + '(' + '&&'.join([catVtag['HP2'],catVtag['LP1']]) + ')' + ')'

else:
 print "Use b-tagging sorting"
 cuts['VH_HPHP'] = '('+  '&&'.join([catHtag['HP1'],catVtag['HP2']]) + ')'
 cuts['VH_HPLP'] = '('+  '&&'.join([catHtag['HP1'],catVtag['LP2']]) + ')'
 cuts['VH_LPHP'] = '('+  '&&'.join([catHtag['LP1'],catVtag['HP2']]) + ')'
 cuts['VH_LPLP'] = '('+  '&&'.join([catHtag['LP1'],catVtag['LP2']]) + ')'
 cuts['VH_all'] =  '('+  '||'.join([cuts['VH_HPHP'],cuts['VH_HPLP'],cuts['VH_LPHP'],cuts['VH_LPLP']]) + ')'
 cuts['VV_HPHP'] = '(' + '!' + cuts['VH_all'] + '&&' + '(' + '&&'.join([catVtag['HP1'],catVtag['HP2']]) + ')' + ')'
 cuts['VV_HPLP'] = '(' + '!' + cuts['VH_all'] + '&&' + '(' + '('+  '&&'.join([catVtag['HP1'],catVtag['LP2']]) + ')' + '||' + '(' + '&&'.join([catVtag['HP2'],catVtag['LP1']]) + ')' + ')' + ')'

#validation regions --> we might need more of these here (control region of b-tagging?)
cuts['VV_LPLP'] = '(' + '&&'.join([catVtag['LP1'],catVtag['LP2']]) + ')'


#categories B2G18002
#cuts['HPHP'] = '('+cat['HP1']+'&&'+cat['HP2']+')'
#cuts['HPLP'] = '(('+cat['HP1']+'&&'+cat['LP2']+')||('+cat['LP1']+'&&'+cat['HP2']+'))'
#cuts['LPLP'] = '('+cat['LP1']+'&&'+cat['LP2']+')'
#cuts['NP'] = '(('+cat['LP1']+'&&'+cat['NP2']+')||('+cat['NP1']+'&&'+cat['LP2']+')||('+cat['NP1']+'&&'+cat['NP2']+'))'

#these will probably change to use mergedHTruth as well? do we need a mergedTopTruth?
cuts['nonres'] = '1'
cuts['res'] = '(jj_l1_mergedVTruth==1&&jj_l1_softDrop_mass>60&&jj_l1_softDrop_mass<120)'
cuts['resTT'] = '(jj_l1_mergedVTruth==1&&jj_l1_softDrop_mass>140&&jj_l1_softDrop_mass<200)'

#all categories
#categories=['VH_HPHP','VH_HPLP','VH_LPHP','VH_LPLP','VV_HPHP','VV_HPLP']
categories=['VV_HPLP']

#list of signal samples --> nb, radion and vbf samples to be added
BulkGravWWTemplate="BulkGravToWW"
BulkGravZZTemplate="BulkGravToZZToZhadZhad"
ZprimeWWTemplate= "ZprimeToWW"
ZprimeZHTemplate="ZprimeToZhToZhadhbb"
WprimeWZTemplate= "WprimeToWZToWhadZhad"
WprimeWHTemplate="WprimeToWhToWhadhbb"

# use arbitrary cross section 0.001 so limits converge better
BRZZ=1.*0.001*0.6991*0.6991
BRWW=1.*0.001 #ZprimeWW and GBulkWW are inclusive
BRZH=1.*0.001*0.6991*0.584
BRWZ=1.*0.001*0.6991*0.676
BRWH=1.*0.001*0.676*0.584

#data samples
dataTemplate="JetHT"

#background samples
nonResTemplate="QCD_Pt-" #low stat herwig
#nonResTemplate="QCD_HT" #medium stat madgraph+pythia
#nonResTemplate="QCD_Pt_" #high stat pythia8
if(period == 2016):
    TTemplate= "TT_Mtt-700to1000,TT_Mtt-1000toInf" #do we need a separate fit for ttbar?
else:
    TTemplate= "TTToHadronic" #do we need a separate fit for ttbar?
WresTemplate= "WJetsToQQ_HT400to600,WJetsToQQ_HT600to800,WJetsToQQ_HT800toInf,"+str(TTemplate)
ZresTemplate= "ZJetsToQQ_HT400to600,ZJetsToQQ_HT600to800,ZJetsToQQ_HT800toInf"
resTemplate= "ZJetsToQQ_HT400to600,ZJetsToQQ_HT600to800,ZJetsToQQ_HT800toInf,WJetsToQQ_HT400to600,WJetsToQQ_HT600to800,WJetsToQQ_HT800toInf,"+str(TTemplate)

#ranges and binning
minMJ=55.0
maxMJ=215.0
binsMJ=80

minMVV=838.0
maxMVV=6000.
binsMVV=100

minMX=1200.0
maxMX=4500.0
    
if dijetBinning:
    minMVV = float(dijetbins[0])
    maxMVV = float(dijetbins[-1])
    binsMVV= len(dijetbins)-1
      
cuts['acceptance']= "(jj_LV_mass>{minMVV}&&jj_LV_mass<{maxMVV}&&jj_l1_softDrop_mass>{minMJ}&&jj_l1_softDrop_mass<{maxMJ}&&jj_l2_softDrop_mass>{minMJ}&&jj_l2_softDrop_mass<{maxMJ})".format(minMVV=minMVV,maxMVV=maxMVV,minMJ=minMJ,maxMJ=maxMJ)
cuts['acceptanceMJ']= "(jj_l1_softDrop_mass>{minMJ}&&jj_l1_softDrop_mass<{maxMJ}&&jj_l2_softDrop_mass>{minMJ}&&jj_l2_softDrop_mass<{maxMJ})".format(minMJ=minMJ,maxMJ=maxMJ)
cuts['acceptanceMVV'] = "(jj_LV_mass>{minMVV}&&jj_LV_mass<{maxMVV})".format(minMVV=minMVV,maxMVV=maxMVV)
cuts['acceptanceGEN']='(jj_l1_gen_softDrop_mass>20.&&jj_l2_gen_softDrop_mass>20.&jj_l1_gen_softDrop_mass<300.&&jj_l2_gen_softDrop_mass<300.&&jj_gen_partialMass>800.&&jj_gen_partialMass<6000.&&TMath::Log(jj_l1_gen_softDrop_mass**2/jj_l1_gen_pt**2)<-1.5&&TMath::Log(jj_l2_gen_softDrop_mass**2/jj_l2_gen_pt**2)<-1.5)'
cuts['looseacceptanceMJ']= "(jj_l1_softDrop_mass>35&&jj_l1_softDrop_mass<300&&jj_l2_softDrop_mass>35&&jj_l2_softDrop_mass<300)"

#do not change the order here, add at the end instead
parameters = [cuts,minMVV,maxMVV,minMX,maxMX,binsMVV,HCALbinsMVV,samples,categories,minMJ,maxMJ,binsMJ,lumi,submitToBatch]   
f = AllFunctions(parameters)

#Fitting steps for one signal sample 
signal_inuse="BulkGravToWW"
signaltemplate_inuse=BulkGravWWTemplate
xsec_inuse=BRWW
f.makeSignalShapesMVV("JJ_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse) #nb, to be optimized
f.makeSignalShapesMJ("JJ_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,'l1')
f.makeSignalShapesMJ("JJ_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,'l2')
f.makeSignalYields("JJ_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,xsec_inuse,{'VH_HPHP':HPSF*HPSF,'VH_HPLP':HPSF*LPSF,'VH_LPHP':HPSF*LPSF,'VH_LPLP':LPSF*LPSF,'VV_HPHP':HPSF*HPSF,'VV_HPLP':HPSF*LPSF})

#Detector response
f.makeDetectorResponse("nonRes","JJ_"+str(period),nonResTemplate,cuts['nonres'])

# Make nonresonant QCD templates and normalization
if runParallel and submitToBatch:
  wait = False
  f.makeBackgroundShapesMVVKernel("nonRes","JJ_"+str(period),nonResTemplate,cuts['nonres'],"1D",wait)
  f.makeBackgroundShapesMVVConditional("nonRes","JJ_"+str(period),nonResTemplate,'l1',cuts['nonres'],"2Dl1",wait)
  f.makeBackgroundShapesMVVConditional("nonRes","JJ_"+str(period),nonResTemplate,'l2',cuts['nonres'],"2Dl2",wait)
  print "Exiting system! When all jobs are finished, please run mergeKernelJobs below"
  sys.exit()
  f.mergeKernelJobs("nonRes","JJ_"+str(period))
else:
  wait = True
  f.makeBackgroundShapesMVVKernel("nonRes","JJ_"+str(period),nonResTemplate,cuts['nonres'],"1D",wait)
  f.makeBackgroundShapesMVVConditional("nonRes","JJ_"+str(period),nonResTemplate,'l1',cuts['nonres'],"2Dl1",wait)
  f.makeBackgroundShapesMVVConditional("nonRes","JJ_"+str(period),nonResTemplate,'l2',cuts['nonres'],"2Dl2",wait)

f.mergeBackgroundShapes("nonRes","JJ_"+str(period))
f.makeNormalizations("nonRes","JJ_"+str(period),nonResTemplate,0,cuts['nonres'],"nRes")

#for V+jets
f.fitVJets("JJ_WJets",resTemplate,1.,1.)
f.makeBackgroundShapesMVVKernel("WJets","JJ_"+str(period),WresTemplate,cuts['nonres'],"1D",0,1.,1.)
f.makeBackgroundShapesMVVKernel("ZJets","JJ_"+str(period),ZresTemplate,cuts['nonres'],"1D",0,1.,1.)
f.makeNormalizations("WJets","JJ_"+str(period),WresTemplate,0,cuts['nonres'],"nRes","",HPSF,LPSF) #this contains also TTbar...
f.makeNormalizations("ZJets","JJ_"+str(period),ZresTemplate,0,cuts['nonres'],"nRes","",HPSF,LPSF)
#f.makeNormalizations("TTJets","JJ_"+str(period),TTemplate,0,cuts['nonres'],"nRes","") # ... so we do not need this 


## Do data or pseudodata
#f.makeNormalizations("data","JJ",dataTemplate,1,'1',"normD") #run on data. Currently run on pseudodata only (below)
#from modules.submitJobs import makePseudoData
#for p in purities: makePseudoData("JJ_nonRes_%s.root"%p,"JJ_nonRes_3D_%s.root"%p,"pythia","JJ_PDnoVjets_%s.root"%p,lumi)
#from modules.submitJobs import makePseudoDataVjets
#for p in purities: makePseudoDataVjets("/afs/cern.ch/user/t/thaarres/public/forJen/looseDDT/JJ_nonRes_%s.root"%p,"/afs/cern.ch/user/t/thaarres/public/forJen/looseDDT/JJ_nonRes_3D_%s.root"%p,"pythia","/afs/cern.ch/user/t/thaarres/public/forJen/looseDDT/JJ_PD_%s.root"%p,lumi,"/afs/cern.ch/user/t/thaarres/public/forJen/looseDDT/workspace_JJ_13TeV_2017.root",2017,p)

print " ########## I did everything I could! ###### "
