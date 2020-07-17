from functions import *
from optparse import OptionParser
#from cuts import cuts, HPSF16, HPSF17, LPSF16, LPSF17, dijetbins, HCALbinsMVVSignal, minMJ,maxMJ,binsMJ, minMVV, maxMVV, binsMVV, minMX, maxMX, catVtag, catHtag
import cuts

## import cuts of the analysis from separate file

# python makeInputs.py -p 2016 --run "detector" --batch False
# python makeInputs.py -p 2016 --run "signorm" --signal "ZprimeWW" --batch False 
# python makeInputs.py -p 2016 --run "tt" --batch False 
# python makeInputs.py -p 2016 --run "vjets" --batch False   
# python makeInputs.py -p "2016,2017,2018"  --run "vjetsAll"
# python makeInputs.py -p "2016,2017,2018"  --run "vjetsfits"
# python makeInputs.py -p "2016,2017,2018"  --run "vjetskernel"
# python makeInputs.py -p "2016,2017,2018"  --run "vjetsnorm"
# python makeInputs.py -p 2016 --run "qcdtemplates"
# python makeInputs.py -p 2016 --run "qcdkernel"
# python makeInputs.py -p "2016,2017,2018" --run "qcdkernel"  --single True
# python makeInputs.py -p 2016 --run "qcdnorm"
# python makeInputs.py -p 2016 --run "data"
# python makeInputs.py -p 2016 --run "pseudoNOVJETS"
# python makeInputs.py -p 2016 --run "pseudoVJETS"                                                                                                                                                                                                                   

parser = OptionParser()
parser.add_option("-p","--period",dest="period",default="2016",help="run period")
parser.add_option("-s","--sorting",dest="sorting",help="b-tag or random sorting",default='random')
parser.add_option("-b","--binning",action="store_false",dest="binning",help="use dijet binning or not",default=True)
parser.add_option("--batch",action="store_false",dest="batch",help="submit to batch or not ",default=True)
parser.add_option("--trigg",action="store_true",dest="trigg",help="add trigger weights or not ",default=False)
parser.add_option("--run",dest="run",help="decide which parts of the code should be run right now possible optoins are: all : run everything, sigmvv: run signal mvv fit sigmj: run signal mj fit, signorm: run signal norm, vjets: run vjets, tt: run ttbar , qcdtemplates: run qcd templates, qcdkernel: run qcd kernel, qcdnorm: run qcd merge and norm, detector: run detector fit , data : run the data or pseudodata scripts ",default="all")
parser.add_option("--signal",dest="signal",default="BGWW",help="which signal do you want to run? options are BulkGWW, BulkGZZ, WprimeWZ, ZprimeWW, ZprimeZH")
parser.add_option("--single",dest="single",default=False,help="set to True to merge kernels also for single years when processing full run2 data")

(options,args) = parser.parse_args()

widerMVV=True


ctx  = cuts.cuts("init_VV_VH.json",options.period,options.sorting+"dijetbins",widerMVV)
if options.binning==False: ctx  = cuts.cuts("init_VV_VH.json",int(options.period),options.sorting,widerMVV)

basedir="deepAK8V2/"
print "options.period  ",options.period  
period = options.period
samples=""
filePeriod=period
rescale=False #for pseudodata: set to True if files comes from splitting Run2
if options.period.find(",")!=-1: 
    lumitot=0
    period = options.period.split(',') 
    filePeriod="Run2"
    rescale=True
    for year in period:
        print year
        if year==period[-1]: samples+=basedir+year+"/"
        else: samples+=basedir+year+"/,"
        lumitot+=ctx.lumi[year]
    ctx.lumi.update({"Run2":lumitot}) #this is used only to make pseudodata with full run2 lumi
else: samples=basedir+period+"/"
# NB to use the DDT decorrelation method, the ntuples in /eos/cms/store/cmst3/group/exovv/VVtuple/FullRun2VVVHNtuple/deepAK8V2/ should be used
#samples= str(period)+"trainingV2/" #for V+jets we use 2017 samples also for 2016 because the 2016 ones are buggy
#samples= str(period)+"/" #for V+jets we use 2017 samples also for 2016 because the 2016 ones are buggy
#samples="deepAK8V2/"+period
#if period=="201678": sample="deepAK8V2/201?/"

print "period ",period
print "sample ",samples

sorting = options.sorting

submitToBatch = options.batch #Set to true if you want to submit kernels + makeData to batch!
runParallel   = True #Set to true if you want to run all kernels in parallel! This will exit this script and you will have to run mergeKernelJobs when your jobs are done! 

dijetBinning = options.binning
useTriggerWeights = options.trigg


    
addOption = ""
if useTriggerWeights: 
    addOption = "-t"
    
#all categories
#categories = ['none']
categories=["VV_HPLP"]
#categories=["VH_HPLP"]
#categories =["NP"]
#categories =["VV_HPLP"]
       

#list of signal samples --> nb, radion and vbf samples to be added
BulkGravWWTemplate="BulkGravToWW_"
VBFBulkGravWWTemplate="VBF_BulkGravToWW_"
BulkGravZZTemplate="BulkGravToZZToZhadZhad_"
ZprimeWWTemplate= "ZprimeToWW_"
ZprimeZHTemplate="ZprimeToZhToZhadhbb_"
WprimeWZTemplate= "WprimeToWZToWhadZhad_"
WprimeWHTemplate="WprimeToWhToWhadhbb_" #"WprimeToWhToWhadhbb_"

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

#commentin TT for the moment, will need better shaping                                                                                                                                                                                       
if(period == "2016"):                                                                                                                                                                                                                        
    TTemplate= "TT_Mtt-700to1000,TT_Mtt-1000toInf" 
elif (filePeriod == "Run2"):
    TTemplate= "TT_Mtt-700to1000,TT_Mtt-1000toInf,TTToHadronic"
else:                                                                                                                                                                                                                                       
    TTemplate= "TTToHadronic" 
#WresTemplate= "WJetsToQQ_HT400to600,WJetsToQQ_HT600to800,WJetsToQQ_HT800toInf,"+str(TTemplate)                                                                                                                                              
WresTemplate= "WJetsToQQ_HT400to600,WJetsToQQ_HT600to800,WJetsToQQ_HT800toInf"
ZresTemplate= "ZJetsToQQ_HT400to600,ZJetsToQQ_HT600to800,ZJetsToQQ_HT800toInf"
resTemplate= "ZJetsToQQ_HT400to600,ZJetsToQQ_HT600to800,ZJetsToQQ_HT800toInf,WJetsToQQ_HT400to600,WJetsToQQ_HT600to800,WJetsToQQ_HT800toInf"

      

#do not change the order here, add at the end instead
parameters = [ctx.cuts,ctx.minMVV,ctx.maxMVV,ctx.minMX,ctx.maxMX,ctx.binsMVV,ctx.HCALbinsMVV,samples,categories,ctx.minMJ,ctx.maxMJ,ctx.binsMJ,ctx.lumi,submitToBatch]   
f = AllFunctions(parameters)


#parser.add_option("--signal",dest="signal",default="BGWW",help="which signal do you want to run? options are BGWW, BGZZ, WprimeWZ, ZprimeWW, ZprimeZH")
if options.run.find("all")!=-1 or options.run.find("sig")!=-1:
    if options.signal.find("ZprimeZH")!=-1:
        signal_inuse="ZprimeZH"
        signaltemplate_inuse=ZprimeZHTemplate
        xsec_inuse=BRZH
    elif options.signal.find("BGWW")!=-1 and not 'VBF' in options.signal:
        signal_inuse="BulkGWW"
        signaltemplate_inuse=BulkGravWWTemplate
        xsec_inuse=BRWW
    elif options.signal.find("VBFBGWW")!=-1:
        signal_inuse="VBF_BulkGWW"
        signaltemplate_inuse=VBFBulkGravWWTemplate
        xsec_inuse=BRWW
    elif options.signal.find("BGZZ")!=-1:
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
    elif options.signal.find("WprimeWH")!=-1:
        signal_inuse="WprimeWH"
        signaltemplate_inuse=WprimeWHTemplate
        xsec_inuse=BRWH
    else:
        print "signal "+str(options.signal)+" not found!"
        sys.exit()


fixParsSig=ctx.fixParsSig

fixParsSigMVV=ctx.fixParsSigMVV


if options.run.find("all")!=-1 or options.run.find("sig")!=-1:
    print "run signal"
    if options.run.find("all")!=-1 or options.run.find("mj")!=-1:
        print "mj fit for signal "
        if sorting == "random":
            if signal_inuse.find("H")!=-1: 
                f.makeSignalShapesMJ("JJ_Vjet_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,'random', fixParsSig[signal_inuse],"jj_random_mergedVTruth==1")
                f.makeSignalShapesMJ("JJ_Hjet_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,'random',fixParsSig[signal_inuse],"jj_random_mergedHTruth==1")
            else:
                f.makeSignalShapesMJ("JJ_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,'random',fixParsSig[signal_inuse.replace('VBF_','')]) 
        else:
            if signal_inuse.find("H")!=-1: 
                f.makeSignalShapesMJ("JJ_Vjet_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,'l1',fixParsSig[signal_inuse],"jj_l1_mergedVTruth==1")
                f.makeSignalShapesMJ("JJ_Vjet_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,'l2',fixParsSig[signal_inuse],"jj_l2_mergedVTruth==1")
                f.makeSignalShapesMJ("JJ_Hjet_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,'l1',fixParsSig[signal_inuse],"jj_l1_mergedHTruth==1")
                f.makeSignalShapesMJ("JJ_Hjet_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,'l2',fixParsSig[signal_inuse],"jj_l2_mergedHTruth==1")
            else:
                f.makeSignalShapesMJ("JJ_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,'l1',fixParsSig[signal_inuse])
                f.makeSignalShapesMJ("JJ_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,'l2',fixParsSig[signal_inuse])
    if options.run.find("all")!=-1 or options.run.find("mvv")!=-1:
        print "mjj fit for signal ",signal_inuse
        if signal_inuse.find("H")!=-1:
            f.makeSignalShapesMVV("JJ_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,fixParsSigMVV[signal_inuse],"( jj_l1_softDrop_mass <= 215 && jj_l1_softDrop_mass > 105 && jj_l2_softDrop_mass <= 105 && jj_l2_softDrop_mass > 55) || (jj_l2_softDrop_mass <= 215 && jj_l2_softDrop_mass > 105 && jj_l1_softDrop_mass <= 105 && jj_l1_softDrop_mass > 55) ")
        elif signal_inuse.find("WZ")!=-1:
            #f.makeSignalShapesMVV("JJ_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,fixParsSigMVV[signal_inuse],"(jj_l1_softDrop_mass <= 105 && jj_l1_softDrop_mass > 85 && jj_l2_softDrop_mass <= 85 && jj_l2_softDrop_mass >= 65) || (jj_l2_softDrop_mass <= 105 && jj_l2_softDrop_mass > 85 && jj_l1_softDrop_mass <= 85 && jj_l1_softDrop_mass >= 65)")
            f.makeSignalShapesMVV("JJ_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,fixParsSigMVV[signal_inuse],"1")
        else:
            f.makeSignalShapesMVV("JJ_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,fixParsSigMVV[signal_inuse.replace('VBF_','')],"1")
    

    if options.run.find("all")!=-1 or options.run.find("norm")!=-1:
        print "fit signal norm "
        f.makeSignalYields("JJ_"+str(signal_inuse)+"_"+filePeriod,signaltemplate_inuse,xsec_inuse,'"pol3"') #'"[0]*TMath::Log10(x)"')
        #f.makeNormalizations("sigonly_M2000","JJ_"+filePeriod+"_"+str(signal_inuse),signaltemplate_inuse+"narrow_2000",0,ctx.cuts['nonres'],"sig")
        #f.makeNormalizations("sigonly_M4000","JJ_"+filePeriod+"_"+str(signal_inuse),signaltemplate_inuse+"narrow_4000",0,ctx.cuts['nonres'],"sig")

if options.run.find("all")!=-1 or options.run.find("detector")!=-1:
    print "make Detector response"
    f.makeDetectorResponse("nonRes","JJ_"+filePeriod,nonResTemplate,ctx.cuts['nonres'])

if options.run.find("all")!=-1 or options.run.find("qcd")!=-1:
    print "Make nonresonant QCD templates and normalization"
    if runParallel and submitToBatch:
        if options.run.find("all")!=-1 or options.run.find("templates")!=-1:
            wait = False
            print "ctx.cuts['nonres'] = ",ctx.cuts['nonres']  
            f.makeBackgroundShapesMVVKernel("nonRes","JJ_"+filePeriod,nonResTemplate,ctx.cuts['nonres'],"1D",wait)
            f.makeBackgroundShapesMVVConditional("nonRes","JJ_"+filePeriod,nonResTemplate,'l1',ctx.cuts['nonres'],"2Dl1",wait)
            f.makeBackgroundShapesMVVConditional("nonRes","JJ_"+filePeriod,nonResTemplate,'l2',ctx.cuts['nonres'],"2Dl2",wait)
            print "Exiting system! When all jobs are finished, please run mergeKernelJobs below"
            sys.exit()
        elif options.run.find("all")!=-1 or options.run.find("kernel")!=-1:
            f.mergeKernelJobs("nonRes","JJ_"+filePeriod,options.single)
            print " calling merge Bckg shape"
	    f.mergeBackgroundShapes("nonRes","JJ_"+filePeriod,options.single)
    else:
        if options.run.find("all")!=-1 or options.run.find("templates")!=-1:
            wait = True
            f.makeBackgroundShapesMVVKernel("nonRes","JJ_"+filePeriod,nonResTemplate,ctx.cuts['nonres'],"1D",wait)
            f.makeBackgroundShapesMVVConditional("nonRes","JJ_"+filePeriod,nonResTemplate,'l1',ctx.cuts['nonres'],"2Dl1",wait)
            f.makeBackgroundShapesMVVConditional("nonRes","JJ_"+filePeriod,nonResTemplate,'l2',ctx.cuts['nonres'],"2Dl2",wait)
            f.mergeBackgroundShapes("nonRes","JJ_"+filePeriod)
    if options.run.find("all")!=-1 or options.run.find("norm")!=-1:
        f.makeNormalizations("nonRes","JJ_"+filePeriod,nonResTemplate,0,ctx.cuts['nonres'],"nRes",options.single)

if options.run.find("all")!=-1 or options.run.find("vjets")!=-1:    
    if options.run.find("all")!=-1 or options.run.find("fits")!=-1 or options.run.find("All")!=-1:
        print "for V+jets"
        print "first we fit"
        f.fitVJets("JJ_"+filePeriod+"_WJets",resTemplate,1.,1.)
    wait=False
    if options.batch == True : wait=True 
    if options.run.find("all")!=-1 or options.run.find("kernel")!=-1 or options.run.find("All")!=-1:
        print " did you run Detector response  for this period? otherwise the kernels steps will not work!"
        print "first kernel W"
        f.makeBackgroundShapesMVVKernel("WJets","JJ_"+filePeriod,WresTemplate,ctx.cuts['nonres'],"1DW",wait,1.,1.)
        print "then kernel Z"
        f.makeBackgroundShapesMVVKernel("ZJets","JJ_"+filePeriod,ZresTemplate,ctx.cuts['nonres'],"1DZ",wait,1.,1.)
    if options.run.find("all")!=-1 or options.run.find("vjetsnorm")!=-1 or options.run.find("All")!=-1:
        print "then norm W"
        f.makeNormalizations("WJets","JJ_"+filePeriod,WresTemplate,0,ctx.cuts['nonres'],"nResWJets",options.single,"1") #,HPSF_vtag,LPSF_vtag)
        print "then norm Z"
        f.makeNormalizations("ZJets","JJ_"+filePeriod,ZresTemplate,0,ctx.cuts['nonres'],"nResZJets",options.single,"1")



if options.run.find("all")!=-1 or options.run.find("tt")!=-1:
    f.fitTT   ("JJ_%s_TTJets"%(filePeriod),TTemplate,1.,)
    wait=False
    if options.batch == True : wait=True
    f.makeBackgroundShapesMVVKernel("TTJets","JJ_"+filePeriod,TTemplate,ctx.cuts['nonres'],"1D",wait,1.,1.)
    f.makeNormalizations("TTJets","JJ_"+filePeriod,TTemplate,0,ctx.cuts['nonres'],"nResTT",options.single,"1")
  
if options.run.find("all")!=-1 or options.run.find("data")!=-1:
    print " Do data "
    f.makeNormalizations("data","JJ_"+filePeriod,dataTemplate,1,'1',"normD") #run on data. Currently run on pseudodata only (below)
if options.run.find("all")!=-1 or options.run.find("pseudoNOVJETS")!=-1:
    print " Do pseudodata without vjets"
    from modules.submitJobs import makePseudoData
    #for p in categories: makePseudoData("results_"+filePeriod+"/JJ_"+filePeriod+"_nonRes_%s.root"%p,"results_"+filePeriod+"/save_new_shapes_"+filePeriod+"_pythia_%s_3D.root"%p,"pythia","JJ_%s_PDnoVjets_%s.root"%(filePeriod,p),ctx.lumi[filePeriod])
    for p in categories: makePseudoData("results_"+filePeriod+"/JJ_"+filePeriod+"_nonRes_%s.root"%p,"results_"+filePeriod+"/JJ_"+filePeriod+"_nonRes_3D_NP.root","pythia","JJ_%s_PDnoVjets_%s.root"%(filePeriod,p),ctx.lumi[filePeriod])
if options.run.find("all")!=-1 or options.run.find("pseudoVJETS")!=-1:
    print " Do pseudodata with vjets: DID YOU PRODUCE THE WORKSPACE BEFORE???"
    from modules.submitJobs import makePseudoDataVjets
    #for p in categories: makePseudoDataVjets("results_"+filePeriod+"/JJ_"+filePeriod+"_nonRes_%s.root"%p,"results_"+filePeriod+"/save_new_shapes_"+filePeriod+"_pythia_%s_3D.root"%p,"pythia","JJ_%s_PDVjets_%s.root"%(filePeriod,p),ctx.lumi[filePeriod],"results_"+filePeriod+"/workspace_JJ_BulkGWW_"+p+"_13TeV_"+filePeriod+"_VV_V_prepVjets.root",filePeriod,p)
    for p in categories: makePseudoDataVjets("results_"+filePeriod+"/JJ_"+filePeriod+"_nonRes_%s.root"%p,"results_"+filePeriod+"/JJ_"+filePeriod+"_nonRes_3D_NP.root","pythia","JJ_%s_PDVjets_%s.root"%(filePeriod,p),ctx.lumi[filePeriod],"results_"+filePeriod+"/workspace_JJ_BulkGWW_"+p+"_13TeV_"+filePeriod+"test_prepVjets.root",filePeriod,p,rescale)
if options.run.find("all")!=-1 or options.run.find("pseudoALL")!=-1:
    print " Do pseudodata. DID YOU PRODUCE THE WORKSPACE BEFORE???"
    from modules.submitJobs import makePseudoDataVjetsTT
    for p in categories: makePseudoDataVjetsTT("results_"+filePeriod+"/JJ_"+filePeriod+"_nonRes_%s.root"%p,"results_"+filePeriod+"/save_new_shapes_"+filePeriod+"_pythia_%s_3D.root"%p,"pythia","JJ_PD_%s.root"%p,lumi,"workspace_JJ_ZprimeZH_%s_13TeV_%sZprimeZH.root"%(p,filePeriod),period,p)


print " ########## I did everything I could! ###### "
