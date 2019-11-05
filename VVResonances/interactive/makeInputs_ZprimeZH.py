from functions import *
from optparse import OptionParser

# python makeInputs.py -p 2016 --run "signorm" --signal "ZprimeWW" --batch False 
# python makeInputs.py -p 2016 --run "vjets" --batch False                                                                                                                                      
# python makeInputs.py -p 2016 --run "qcdtemplates"
# python makeInputs.py -p 2016 --run "qcdkernel"
# python makeInputs.py -p 2016 --run "qcdnorm"
# python makeInputs.py -p 2016 --run "data"
# python makeInputs_doubleB.py -p 2016 --run "pseudoNOVJETS"
# python makeInputs_doubleB.py -p 2016 --run "pseudoVJETS"                                                                                                                                                                                                                   

parser = OptionParser()
parser.add_option("-p","--period",dest="period",type="int",default=2016,help="run period")
parser.add_option("-s","--sorting",dest="sorting",help="b-tag or random sorting",default='random')
parser.add_option("-b","--binning",action="store_false",dest="binning",help="use dijet binning or not",default=True)
parser.add_option("--batch",action="store_false",dest="batch",help="submit to batch or not ",default=True)
parser.add_option("--trigg",action="store_true",dest="trigg",help="add trigger weights or not ",default=False)
parser.add_option("--run",dest="run",help="decide which parts of the code should be run right now possible optoins are: all : run everything, sigmvv: run signal mvv fit sigmj: run signal mj fit, signorm: run signal norm, vjets: run vjets , qcdtemplates: run qcd templates, qcdkernel: run qcd kernel, qcdnorm: run qcd merge and norm, detector: run detector fit , data : run the data or pseudodata scripts ",default="all")
parser.add_option("--signal",dest="signal",default="BGWW",help="which signal do you want to run? options are BulkGWW, BulkGZZ, WprimeWZ, ZprimeWW, ZprimeZH")


(options,args) = parser.parse_args()

print options

period = options.period
samples= str(period)+"/" #for V+jets we use 2018 samples also for 2016 because the 2016 ones are buggy

sorting = options.sorting

submitToBatch = options.batch #Set to true if you want to submit kernels + makeData to batch!
runParallel   = True #Set to true if you want to run all kernels in parallel! This will exit this script and you will have to run mergeKernelJobs when your jobs are done! 

dijetBinning = options.binning
useTriggerWeights = options.trigg

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
print "################     you are using TAU21DDT !!!!! #########" 
catVtag['HP1'] = '(jj_l1_tau2/jj_l1_tau1+(0.080*TMath::Log((jj_l1_softDrop_mass*jj_l1_softDrop_mass)/jj_l1_pt)))<0.43'                                                    
catVtag['HP2'] = '(jj_l2_tau2/jj_l2_tau1+(0.080*TMath::Log((jj_l2_softDrop_mass*jj_l2_softDrop_mass)/jj_l2_pt)))<0.43'                                                                                                                                                         
catVtag['LP1'] = '(jj_l1_tau2/jj_l1_tau1+(0.080*TMath::Log((jj_l1_softDrop_mass*jj_l1_softDrop_mass)/jj_l1_pt)))>0.43&&(jj_l1_tau2/jj_l1_tau1+(0.080*TMath::Log((jj_l1_softDrop_mass*jj_l1_softDrop_mass)/jj_l1_pt)))<0.79'                                                   
catVtag['LP2'] = '(jj_l2_tau2/jj_l2_tau1+(0.080*TMath::Log((jj_l2_softDrop_mass*jj_l2_softDrop_mass)/jj_l2_pt)))>0.43&&(jj_l2_tau2/jj_l2_tau1+(0.080*TMath::Log((jj_l2_softDrop_mass*jj_l2_softDrop_mass)/jj_l2_pt)))<0.79'                                                    
catVtag['NP1'] = '(jj_l1_tau2/jj_l1_tau1+(0.080*TMath::Log((jj_l1_softDrop_mass*jj_l1_softDrop_mass)/jj_l1_pt)))>0.79'                                                                                                                                                         
catVtag['NP2'] = '(jj_l2_tau2/jj_l2_tau1+(0.080*TMath::Log((jj_l2_softDrop_mass*jj_l2_softDrop_mass)/jj_l2_pt)))>0.79'                                                                                                                                                         
'''
#MassDecorrelatedDeepBoosted_WvsQCD   
print "################     you are using DeepAK8 WvsQCD  !!!!! #########"                                                                                                                                                                                                     
catVtag['HP1'] = 'jj_l1_MassDecorrelatedDeepBoosted_WvsQCD>0.8'
catVtag['HP2'] = 'jj_l2_MassDecorrelatedDeepBoosted_WvsQCD>0.8'
catVtag['LP1'] = 'jj_l1_MassDecorrelatedDeepBoosted_WvsQCD<0.8&&jj_l1_MassDecorrelatedDeepBoosted_WvsQCD>0.5'
catVtag['LP2'] = 'jj_l2_MassDecorrelatedDeepBoosted_WvsQCD<0.8&&jj_l2_MassDecorrelatedDeepBoosted_WvsQCD>0.5'
catVtag['NP1'] = 'jj_l1_MassDecorrelatedDeepBoosted_WvsQCD<0.5'
catVtag['NP2'] = 'jj_l2_MassDecorrelatedDeepBoosted_WvsQCD<0.5'
'''

'''
print "################     you are using DeepAK8 ZHbbvsQCD  !!!!! #########"                                                                                                                                                                                             
catHtag['HP1'] = '(jj_l1_MassDecorrelatedDeepBoosted_ZHbbvsQCD>0.8)'
catHtag['HP2'] = '(jj_l2_MassDecorrelatedDeepBoosted_ZHbbvsQCD>0.8)'
catHtag['LP1'] = '(jj_l1_MassDecorrelatedDeepBoosted_ZHbbvsQCD>0.5&&jj_l1_MassDecorrelatedDeepBoosted_ZHbbvsQCD<0.8)'
catHtag['LP2'] = '(jj_l2_MassDecorrelatedDeepBoosted_ZHbbvsQCD>0.5&&jj_l2_MassDecorrelatedDeepBoosted_ZHbbvsQCD<0.8)'
catHtag['NP1'] = '(jj_l1_MassDecorrelatedDeepBoosted_ZHbbvsQCD<0.5)'
catHtag['NP2'] = '(jj_l2_MassDecorrelatedDeepBoosted_ZHbbvsQCD<0.5)'
'''

print "################     you are using double B  !!!!! #########"  
valueHP = 0.92
valueLP = 0.7
catHtag['HP1'] = '(jj_l1_MassIndependentDeepDoubleBvLJetprobHbb_>%.2f)'%valueHP
catHtag['HP2'] = '(jj_l2_MassIndependentDeepDoubleBvLJet_probHbb>%.2f)'%valueHP
catHtag['LP1'] = '(jj_l1_MassIndependentDeepDoubleBvLJetprobHbb_>%.2f&&jj_l1_MassIndependentDeepDoubleBvLJetprobHbb_<%.2f)'%(valueLP,valueHP)
catHtag['LP2'] = '(jj_l2_MassIndependentDeepDoubleBvLJet_probHbb>%.2f&&jj_l2_MassIndependentDeepDoubleBvLJet_probHbb<%.2f)'%(valueLP,valueHP)
catHtag['NP1'] = '(jj_l1_MassIndependentDeepDoubleBvLJetprobHbb_<%.2f)'%valueLP
catHtag['NP2'] = '(jj_l2_MassIndependentDeepDoubleBvLJet_probHbb<%.2f)'%valueLP

cuts={}

#cuts['common'] = '((HLT_JJ)*(run>500) + (run<500))*(passed_METfilters&&passed_PVfilter&&njj>0&&jj_LV_mass>700&&abs(jj_l1_eta-jj_l2_eta)<1.3&&jj_l1_softDrop_mass>0.&&jj_l2_softDrop_mass>0.)' #without rho
cuts['common'] = '((HLT_JJ)*(run>500) + (run<500))*(passed_METfilters&&passed_PVfilter&&njj>0&&jj_LV_mass>700&&abs(jj_l1_eta-jj_l2_eta)<1.3&&jj_l1_softDrop_mass>0.&&jj_l2_softDrop_mass>0.&&TMath::Log(jj_l1_softDrop_mass**2/jj_l1_pt**2)<-1.8&&TMath::Log(jj_l2_softDrop_mass**2/jj_l2_pt**2)<-1.8)'
#cuts['common_VV'] = '((HLT_JJ)*(run>500) + (run<500))*(passed_METfilters&&passed_PVfilter&&njj&&!njj_vbf&&jj_LV_mass>700&&abs(jj_l1_eta-jj_l2_eta)<1.3&&jj_l1_softDrop_mass>0.&&jj_l2_softDrop_mass>0.&&TMath::Log(jj_l1_softDrop_mass**2/jj_l1_pt**2)<-1.8&&TMath::Log(jj_l2_softDrop_mass**2/jj_l2_pt**2)<-1.8)'
#cuts['common_VBF'] = '((HLT_JJ)*(run>500) + (run<500))*(passed_METfilters&&passed_PVfilter&&njj_vbf&&jj_LV_mass>700&&abs(jj_l1_eta-jj_l2_eta)<1.3&&jj_l1_softDrop_mass>0.&&jj_l2_softDrop_mass>0.&&TMath::Log(jj_l1_softDrop_mass**2/jj_l1_pt**2)<-1.8&&TMath::Log(jj_l2_softDrop_mass**2/jj_l2_pt**2)<-1.8)'

#signal regions
if sorting == 'random':
 print "Use random sorting!"
 print "current VV + VH"
 cuts['VH_HPHP'] = '('+'('+'('+ '&&'.join([catVtag['HP1'],catHtag['HP2']])+')'+'||'+'('+'&&'.join([catHtag['HP1'],catHtag['HP2']]) +')'+')'+ '||' + '(' + '(' + '&&'.join([catVtag['HP2'],catHtag['HP1']]) +')'+ '||' +'(' + '&&'.join([catHtag['HP2'],catHtag['HP1']]) + ')' + ')' + ')' #anna                                                                                                                       
 cuts['VH_HPLP'] = '(' + '('+  '&&'.join([catVtag['HP1'],catHtag['LP2']]) + ')' + '||' + '(' + '&&'.join([catVtag['HP2'],catHtag['LP1']]) + ')' + ')'
 cuts['VH_LPHP'] = '('+ '('+  '&&'.join([catVtag['LP1'],catHtag['HP2']]) + ')' + '||' + '(' + '&&'.join([catVtag['LP2'],catHtag['HP1']]) + ')' + '||' + '(' + '&&'.join([catHtag['LP2'],catHtag['HP1']]) + ')'+'||'+ '(' + '&&'.join([catHtag['LP1'],catHtag['HP2']]) + ')'+')'
 cuts['VH_LPLP'] = '('+'('+ '&&'.join([catVtag['LP1'],catHtag['LP2']]) + ')'+ '||' +'(' + '&&'.join([catVtag['LP2'],catHtag['LP1']]) + ')' + '||' + '(' + '&&'.join([catVtag['LP2'],catHtag['LP1']]) + ')' +'||' + '(' + '&&'.join([catHtag['LP2'],catHtag['LP1']]) + ')' + ')'
 cuts['VH_all'] =  '('+  '||'.join([cuts['VH_HPHP'],cuts['VH_HPLP'],cuts['VH_LPHP'],cuts['VH_LPLP']]) + ')'
 print " ################ VV orthogonal  VH CATEGORY ############"
 cuts['VV_HPHP'] = '('+'!'+cuts['VH_all']+'&&'+'(' +  '&&'.join([catVtag['HP1'],catVtag['HP2']]) + ')' + ')'
 cuts['VV_HPLP'] = '(' + '!' + cuts['VH_all'] + '&&' + '(' + '('+  '&&'.join([catVtag['HP1'],catVtag['LP2']]) + ')' + '||' + '(' + '&&'.join([catVtag['HP2'],catVtag['LP1']]) + ')' + ')' + ')'

# print " ################ VV inclusive !! ! ############"
# cuts['VV_HPHP'] = '(' +  '&&'.join([catVtag['HP1'],catVtag['HP2']]) + ')' 
# cuts['VV_HPLP'] =  '(' + '('+  '&&'.join([catVtag['HP1'],catVtag['LP2']]) + ')' + '||' + '(' + '&&'.join([catVtag['HP2'],catVtag['LP1']]) + ')' + ')' 
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
#categories=['VH_HPHP','VH_HPLP','VH_LPHP','VV_HPHP','VV_HPLP'] #,'VBF_VV_HPHP','VBF_VV_HPLP']
categories=['VH_LPHP'] #,'VH_LPHP','VV_HPHP','VV_HPLP'] 
#categories=['VH_HPHP','VH_LPHP','VV_HPHP','VV_HPLP'] 
#categories=['VV_HPHP','VV_HPLP']
                                                                                                                                                                                   
#list of signal samples --> nb, radion and vbf samples to be added
BulkGravWWTemplate="BulkGravToWW_"
VBFBulkGravWWTemplate="VBF_BulkGravToWW_"
BulkGravZZTemplate="BulkGravToZZToZhadZhad_"
ZprimeWWTemplate= "ZprimeToWW_"
ZprimeZHTemplate="ZprimeToZhToZhadhbb_"
WprimeWZTemplate= "WprimeToWZToWhadZhad_"
WprimeWHTemplate="WprimeToWhToWhadhbb_"

# use arbitrary cross section 0.001 so limits converge better
BRZZ=1.*0.001*0.6991*0.6991
BRWW=1.*0.001 #ZprimeWW and GBulkWW are inclusive
BRZH=1.*0.001*0.6991*0.584
BRWZ=1.*0.001*0.6991*0.676
BRWH=1.*0.001*0.676*0.584

#data samples
dataTemplate="JetHT"

#background samples
#nonResTemplate="QCD_Pt-" #low stat herwig
#nonResTemplate="QCD_HT" #medium stat madgraph+pythia
nonResTemplate="QCD_Pt_" #high stat pythia8

print "TTbar is commented out!!! "
#commentin TT for the moment, will need better shaping
#if(period == 2016):
#    TTemplate= "TT_Mtt-700to1000,TT_Mtt-1000toInf" #do we need a separate fit for ttbar?
#else:
#    TTemplate= "TTToHadronic" #do we need a separate fit for ttbar?
#WresTemplate= "WJetsToQQ_HT400to600,WJetsToQQ_HT600to800,WJetsToQQ_HT800toInf,"+str(TTemplate)
WresTemplate= "WJetsToQQ_HT400to600,WJetsToQQ_HT600to800,WJetsToQQ_HT800toInf"
ZresTemplate= "ZJetsToQQ_HT400to600,ZJetsToQQ_HT600to800,ZJetsToQQ_HT800toInf"
#resTemplate= "ZJetsToQQ_HT400to600,ZJetsToQQ_HT600to800,ZJetsToQQ_HT800toInf,WJetsToQQ_HT400to600,WJetsToQQ_HT600to800,WJetsToQQ_HT800toInf,"+str(TTemplate)
resTemplate= "ZJetsToQQ_HT400to600,ZJetsToQQ_HT600to800,ZJetsToQQ_HT800toInf,WJetsToQQ_HT400to600,WJetsToQQ_HT600to800,WJetsToQQ_HT800toInf"


#ranges and binning
minMJ=55.0
maxMJ=215.0
binsMJ=20

minMVV=838.0
maxMVV=6000.
binsMVV=100

minMX=1200.0
maxMX=4000.0
    
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


#parser.add_option("--signal",dest="signal",default="BGWW",help="which signal do you want to run? options are BGWW, BGZZ, WprimeWZ, ZprimeWW, ZprimeZH")
if options.run.find("all")!=-1 or options.run.find("sig")!=-1:
    if options.signal.find("ZprimeZH")!=-1:
        signal_inuse="ZprimeZH"
        signaltemplate_inuse=ZprimeZHTemplate
        xsec_inuse=BRZH
    elif options.signal.find("BulkGWW")!=-1:
        signal_inuse="BulkGWW"
        signaltemplate_inuse=BulkGravWWTemplate
        xsec_inuse=BRWW
    elif options.signal.find("VBFBGWW")!=-1:
        signal_inuse="VBF_BulkGWW"
        signaltemplate_inuse=VBFBulkGravWWTemplate
        xsec_inuse=BRWW
    elif options.signal.find("BulkGZZ")!=-1:
        signal_inuse="BulkGZZ"
        signaltemplate_inuse=BulkGravZZTemplate
        xsec_inuse=BRZZ
    elif options.signal.find("ZprimeWW")!=-1:
        signal_inuse="ZprimeWW"
        signaltemplate_inuse=ZprimeWWTemplate
        xsec_inuse=BRWW
    elif options.signal.find("WprimeWZ")!=-1:
        signal_inuse="WprimeWZ"
        signaltemplate_inuse=WprimeWZTemplate
        xsec_inuse=BRWZ
    else:
        print "signal "+str(options.signal)+" not found!"
        sys.exit()


fixParsSig={"ZprimeZH":{
    "VV_HPLP": {"fixPars":"mean:91.5,n:1.83,n2:4.22,sigmaH:10.7,nH:130", "pol":"mean:pol0,sigma:pol5,alpha:pol5,n:pol0,alpha2:pol5,n2:pol0,meanH:pol4,sigmaH:pol0,alphaH:pol2,nH:pol3,alpha2H:pol3,n2H:pol4"}, 
    "VH_all": {"fixPars":"mean:91.5,n2:4.22,n:128,alphaH:0.51,nH:127","pol":"mean:pol0,sigma:pol5,alpha:pol5,n:pol0,alpha2:pol5,n2:pol0,meanH:pol5,sigmaH:pol7,alphaH:pol0,nH:pol3,alpha2H:pol3,n2H:pol4"}, 
    "VH_HPLP": {"fixPars":"mean:90.5,sigmaH:10,n:5,nH:5", "pol":"mean:pol0,sigma:pol5,alpha:pol5,n:pol0,alpha2:pol5,n2:pol3,meanH:pol5,sigmaH:pol6,alphaH:pol3,nH:pol3,alpha2H:pol5,n2H:pol4"},
    "VH_LPHP": {"fixPars":"mean:90.5,sigmaH:10,n:5,nH:5", "pol":"mean:pol0,sigma:pol5,alpha:pol5,n:pol0,alpha2:pol5,n2:pol3,meanH:pol5,sigmaH:pol6,alphaH:pol3,nH:pol3,alpha2H:pol5,n2H:pol4"},#irene
    "VV_HPHP": {"fixPars":"mean:90.9,alpha:1.1,n:1.83,n2:4.22,alphaH:0.5,nH:120", "pol":"mean:pol0,sigma:pol5,alpha:pol5,n:pol0,alpha2:pol5,n2:pol0,meanH:pol4,sigmaH:pol2,alphaH:pol0,nH:pol3,alpha2H:pol4,n2H:pol4"}, 
    "VH_HPHP": {"fixPars":"n:4.2,nH:132", "pol":"mean:pol3,sigma:pol5,alpha:pol5,n:pol0,alpha2:pol3,n2:pol3,meanH:pol5,sigmaH:pol6,alphaH:pol2,nH:pol0,alpha2H:pol3,n2H:pol4"} },
"BulkGWW":{ "VV_HPLP": {"fixPars":"alpha:1.125,n:2,n2:2","pol":"mean:pol4,sigma:pol3,alpha:pol3,n:pol0,alpha2:pol3,n2:pol3"},
            "VV_HPHP": {"fixPars":"alpha:1.08,n:6,n2:2", "pol":"mean:pol5,sigma:pol5,alpha:pol0,n:pol0,alpha2:pol5,n2:pol0"},
            "VH_HPLP": {"fixPars":"alpha:1.125,n:2,n2:2","pol":"mean:pol4,sigma:pol3,alpha:pol3,n:pol0,alpha2:pol3,n2:pol3"},
            "VH_HPHP": {"fixPars":"n:60,alpha:0.76", "pol":"mean:pol5,sigma:pol6,alpha:pol0,n:pol0,alpha2:pol5,n2:pol5"},
            "VH_LPHP": {"fixPars":"alpha:1.125,n:2,n2:2","pol":"mean:pol4,sigma:pol3,alpha:pol3,n:pol0,alpha2:pol3,n2:pol3"}}, #VH_LPHP irene
"BulkGZZ":{"VV_HPLP":{"fixPars":"alpha:1.024,n:3.25","pol":"mean:pol4,sigma:pol3,alpha:pol0,n:pol0,alpha2:pol3,n2:pol4"},
           "VV_HPHP":{"fixPars":"n2:4.8,n:2.8","pol":"mean:pol5,sigma:pol5,alpha:pol5,n:pol0,alpha2:pol3,n2:pol0"},
           "VH_HPLP":{"fixPars":"alpha:1.024,n:3.25","pol":"mean:pol4,sigma:pol3,alpha:pol0,n:pol0,alpha2:pol3,n2:pol4"},
           "VH_HPHP":{"fixPars":"n:64","pol":"mean:pol3,sigma:pol5,alpha:pol5,n:pol0,alpha2:pol3,n2:pol3"}},
"ZprimeWW":{"VV_HPLP": {"fixPars":"alpha:1.125","pol":"mean:pol5,sigma:pol5,alpha:pol0,n:pol3,alpha2:pol3,n2:pol3"},
            "VV_HPHP": {"fixPars":"alpha:1.083,n:3.5,n2:2.3","pol":"mean:pol5,sigma:pol4,alpha:pol0,n:pol0,alpha2:pol5,n2:pol0"},
            "VH_HPLP": {"fixPars":"n:5","pol":"mean:pol5,sigma:pol5,alpha:pol5,n:pol0,alpha2:pol5,n2:pol3"},
            "VH_HPHP": {"fixPars":"n:1.2,alpha:1.22", "pol":"mean:pol6,sigma:pol4,alpha:pol0,n:pol0,alpha2:pol5,n2:pol3"}},
"WprimeWZ":{"VV_HPLP":{"fixPars":"n:2.3","pol":"mean:pol3,sigma:pol3,alpha:pol3,n:pol0,alpha2:pol3,n2:pol1"},
            "VV_HPHP":{"fixPars":"n:2,n2:2,alpha:1.505", "pol":"mean:pol3,sigma:pol3,alpha:pol3,n:pol0,alpha2:pol3,n2:pol1"},
            "VH_HPLP":{"fixPars":"n:2.3","pol":"mean:pol3,sigma:pol3,alpha:pol3,n:pol0,alpha2:pol3,n2:pol1"},
            "VH_HPHP":{"fixPars":"n:0.24,alpha:1.6", "pol":"mean:pol3,sigma:pol5,alpha:pol3,n:pol0,alpha2:pol3,n2:pol3"}}}


#fixPars = "N1:1.61364,N2:4.6012"  
#fixParsSigMVV={"ZprimeZH":{"fixPars":"ALPHA2:2.42,N1:126.5", "pol":"MEAN:pol1,SIGMA:pol1,N1:pol0,ALPHA1:pol9,N2:pol3,ALPHA2:pol0"},"WprimeWZ":{"fixPars":"N1:7,N2:4","pol": "MEAN:pol1,SIGMA:pol3,N1:pol0,ALPHA1:pol7,N2:pol0,ALPHA2:pol5"}}
#fixParsSigMVV={"ZprimeZH":{"fixPars":"ALPHA2:2.42,N1:126.5", "pol":"MEAN:pol1,SIGMA:pol1,N1:pol0,ALPHA1:pol9,N2:pol3,ALPHA2:pol0"},"WprimeWZ":{"fixPars":"N1:7,N2:4","pol": "MEAN:pol1,SIGMA:pol3,N1:pol0,ALPHA1:pol7,N2:pol0,ALPHA2:pol5"},"BulkGWW":{"fixPars":"N1:1.61364,N2:4.6012","pol":"MEAN:pol1,SIGMA:pol6,ALPHA1:pol5,N1:pol0,ALPHA2:pol4,N2:pol0"},"BulkGZZ":{"fixPars":"N1:1.61364,N2:4.6012","pol":"MEAN:pol1,SIGMA:pol6,ALPHA1:pol5,N1:pol0,ALPHA2:pol4,N2:pol0"},"ZprimeWW":{"fixPars":"N1:1.61364,N2:4.6012","pol":"MEAN:pol1,SIGMA:pol6,ALPHA1:pol5,N1:pol0,ALPHA2:pol4,N2:pol0"}}

fixParsSigMVV={"ZprimeZH":{"fixPars":"ALPHA2:2.42,N1:126.5", "pol":"MEAN:pol1,SIGMA:pol1,N1:pol0,ALPHA1:pol5,N2:pol3,ALPHA2:pol0,corr_mean:pol1,corr_sigma:pol1"},              
"WprimeWZ":{"fixPars":"N1:7,N2:4","pol": "MEAN:pol1,SIGMA:pol3,N1:pol0,ALPHA1:pol7,N2:pol0,ALPHA2:pol5,corr_mean:pol1,corr_sigma:pol1"},
"BulkGWW":{"fixPars":"N1:1.61364,N2:4.6012","pol":"MEAN:pol1,SIGMA:pol6,ALPHA1:pol5,N1:pol0,ALPHA2:pol4,N2:pol0"},
"BulkGZZ":{"fixPars":"N1:1.61364,N2:4.6012","pol":"MEAN:pol1,SIGMA:pol6,ALPHA1:pol5,N1:pol0,ALPHA2:pol4,N2:pol0"},
"ZprimeWW":{"fixPars":"N1:1.61364,N2:4.6012","pol":"MEAN:pol1,SIGMA:pol6,ALPHA1:pol5,N1:pol0,ALPHA2:pol4,N2:pol0"}}


if options.run.find("all")!=-1 or options.run.find("sig")!=-1:
    print "run signal"
    if options.run.find("all")!=-1 or options.run.find("mj")!=-1:
        print "mj fit for signal "
        if sorting == "random":
            if signal_inuse.find("H")!=-1: 
                f.makeSignalShapesMJ("JJ_Vjet_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,'random', fixParsSig[signal_inuse],"jj_random_mergedVTruth==1")
                f.makeSignalShapesMJ("JJ_Hjet_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,'random',fixParsSig[signal_inuse],"jj_random_mergedHTruth==1")
            else:
                f.makeSignalShapesMJ("JJ_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,'random',fixParsSig[signal_inuse.replace('VBF_','')]) 
        else:
            if signal_inuse.find("H")!=-1: 
                f.makeSignalShapesMJ("JJ_Vjet_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,'l1',fixParsSig[signal_inuse],"jj_l1_mergedVTruth==1")
                f.makeSignalShapesMJ("JJ_Vjet_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,'l2',fixParsSig[signal_inuse],"jj_l2_mergedVTruth==1")
                f.makeSignalShapesMJ("JJ_Hjet_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,'l1',fixParsSig[signal_inuse],"jj_l1_mergedHTruth==1")
                f.makeSignalShapesMJ("JJ_Hjet_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,'l2',fixParsSig[signal_inuse],"jj_l2_mergedHTruth==1")
            else:
                f.makeSignalShapesMJ("JJ_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,'l1',fixParsSig[signal_inuse])
                f.makeSignalShapesMJ("JJ_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,'l2',fixParsSig[signal_inuse])
    if options.run.find("all")!=-1 or options.run.find("mvv")!=-1:
        print "mjj fit for signal ",signal_inuse
        if signal_inuse.find("H")!=-1:
            f.makeSignalShapesMVV("JJ_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,fixParsSigMVV[signal_inuse],"( jj_l1_softDrop_mass <= 150 && jj_l1_softDrop_mass > 105 && jj_l2_softDrop_mass <= 105 && jj_l2_softDrop_mass > 65) || (jj_l2_softDrop_mass <= 150 && jj_l2_softDrop_mass > 105 && jj_l1_softDrop_mass <= 105 && jj_l1_softDrop_mass > 65) ")
        elif signal_inuse.find("WZ")!=-1:
            f.makeSignalShapesMVV("JJ_j1"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,fixParsSigMVV[signal_inuse],"(jj_l1_softDrop_mass <= 105 && jj_l1_softDrop_mass > 85 && jj_l2_softDrop_mass <= 85 && jj_l2_softDrop_mass >= 65) || (jj_l2_softDrop_mass <= 105 && jj_l2_softDrop_mass > 85 && jj_l1_softDrop_mass <= 85 && jj_l1_softDrop_mass >= 65)")
        else:
            f.makeSignalShapesMVV("JJ_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,fixParsSigMVV[signal_inuse.replace('VBF_','')])
    

    if options.run.find("all")!=-1 or options.run.find("norm")!=-1:
        print "fit signal norm "
        f.makeSignalYields("JJ_"+str(signal_inuse)+"_"+str(period),signaltemplate_inuse,xsec_inuse,{'VH_HPHP':HPSF*HPSF,'VH_HPLP':HPSF*LPSF,'VH_LPHP':HPSF*LPSF,'VH_LPLP':LPSF*LPSF,'VV_HPHP':HPSF*HPSF,'VV_HPLP':HPSF*LPSF,'VH_all':HPSF*HPSF+HPSF*LPSF})
        #f.makeNormalizations(signal_inuse,"JJ_M2000_"+str(period),signaltemplate_inuse+"narrow_2000",0,cuts['nonres'],"nRes")

if options.run.find("all")!=-1 or options.run.find("detector")!=-1:
    print "make Detector response"
    f.makeDetectorResponse("nonRes","JJ_"+str(period),nonResTemplate,cuts['nonres'])

if options.run.find("all")!=-1 or options.run.find("qcd")!=-1:
    print "Make nonresonant QCD templates and normalization"
    if runParallel and submitToBatch:
        if options.run.find("all")!=-1 or options.run.find("templates")!=-1:
            wait = False
            f.makeBackgroundShapesMVVKernel("nonRes","JJ_"+str(period),nonResTemplate,cuts['nonres'],"1D",wait)
            f.makeBackgroundShapesMVVConditional("nonRes","JJ_"+str(period),nonResTemplate,'l1',cuts['nonres'],"2Dl1",wait)
            f.makeBackgroundShapesMVVConditional("nonRes","JJ_"+str(period),nonResTemplate,'l2',cuts['nonres'],"2Dl2",wait)
            print "Exiting system! When all jobs are finished, please run mergeKernelJobs below"
            sys.exit()
        elif options.run.find("all")!=-1 or options.run.find("kernel")!=-1:
            f.mergeKernelJobs("nonRes","JJ_"+str(period))
	    f.mergeBackgroundShapes("nonRes","JJ_"+str(period))
    else:
        if options.run.find("all")!=-1 or options.run.find("templates")!=-1:
            wait = True
            f.makeBackgroundShapesMVVKernel("nonRes","JJ_"+str(period),nonResTemplate,cuts['nonres'],"1D",wait)
            f.makeBackgroundShapesMVVConditional("nonRes","JJ_"+str(period),nonResTemplate,'l1',cuts['nonres'],"2Dl1",wait)
            f.makeBackgroundShapesMVVConditional("nonRes","JJ_"+str(period),nonResTemplate,'l2',cuts['nonres'],"2Dl2",wait)
            f.mergeBackgroundShapes("nonRes","JJ_"+str(period))
    if options.run.find("all")!=-1 or options.run.find("norm")!=-1:
        f.makeNormalizations("nonRes","JJ_"+str(period),nonResTemplate,0,cuts['nonres'],"nRes")

if options.run.find("all")!=-1 or options.run.find("vjets")!=-1:    
    print "for V+jets"
    print "first we fit"
    f.fitVJets("JJ_WJets",resTemplate,1.,1.)
    print "and then we make kernels"
    print " did you run Detector response  for this period? otherwise the kernels steps will not work!"
    print "first kernel W"
    f.makeBackgroundShapesMVVKernel("WJets","JJ_"+str(period),WresTemplate,cuts['nonres'],"1D",0,1.,1.)
    print "then kernel Z"
    f.makeBackgroundShapesMVVKernel("ZJets","JJ_"+str(period),ZresTemplate,cuts['nonres'],"1D",0,1.,1.)
    print "then norm W"
    f.makeNormalizations("WJets","JJ_"+str(period),WresTemplate,0,cuts['nonres'],"nRes","",HPSF,LPSF)
    print "then norm Z"
    f.makeNormalizations("ZJets","JJ_"+str(period),ZresTemplate,0,cuts['nonres'],"nRes","",HPSF,LPSF)
#    f.makeNormalizations("TTJets","JJ_"+str(period),TTemplate,0,cuts['nonres'],"nRes","") # ... so we do not need this


if options.run.find("all")!=-1 or options.run.find("data")!=-1:
    print " Do data "
    f.makeNormalizations("data","JJ_"+str(period),dataTemplate,1,'1',"normD") #run on data. Currently run on pseudodata only (below)
if options.run.find("all")!=-1 or options.run.find("pseudoNOVJETS")!=-1:
    print " Do pseudodata without vjets"
    from modules.submitJobs import makePseudoData
    for p in categories: makePseudoData("JJ_"+str(period)+"_nonRes_%s.root"%p,"save_new_shapes_"+str(period)+"_pythia_%s_3D.root"%p,"pythia","JJ_PDnoVjets_%s.root"%p,lumi)
if options.run.find("all")!=-1 or options.run.find("pseudoVJETS")!=-1:
    print " Do pseudodata with vjets: DID YOU PRODUCE THE WORKSPACE BEFORE???"
    from modules.submitJobs import makePseudoDataVjets
    for p in categories: makePseudoDataVjets("results_"+str(period)+"/JJ_"+str(period)+"_nonRes_%s.root"%p,"results_"+str(period)+"/save_new_shapes_"+str(period)+"_pythia_%s_3D.root"%p,"pythia","JJ_PDVjets_%s.root"%p,lumi,"results_"+str(period)+"/workspace_JJ_BulkGWW__13TeV_"+str(period)+"_prepVjPD.root",period,p)


print " ########## I did everything I could! ###### "
