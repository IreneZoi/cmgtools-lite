import os,sys

class AllFunctions():

 def __init__(self,parameters):
  self.cuts = parameters[0]
  self.minMVV = parameters[1]
  self.maxMVV = parameters[2]
  self.minMX = parameters[3]
  self.maxMX = parameters[4]
  self.binsMVV = parameters[5]
  self.HCALbinsMVV = parameters[6]
  self.samples = parameters[7]
  self.categories = parameters[8]
  self.minMJ = parameters[9]
  self.maxMJ = parameters[10]
  self.binsMJ = parameters[11]
  self.submitToBatch = parameters[12]

  self.printAllParameters()
  
 def makeSignalShapesMVV(self,filename,template):
  print "*** makeSignalShapesMVV start ****" 
  cut='*'.join([self.cuts['common'],self.cuts['acceptanceMJ']])
   
  #the parameters to be fixed should be optimized
  rootFile=filename+"_MVV.root"
  fixPars = "N1:1.61364,N2:4.6012"  
  cmd='vvMakeSignalMVVShapes.py -s "{template}" -c "{cut}"  -o "{rootFile}" -V "jj_LV_mass" --fix "{fixPars}"   -m {minMVV} -M {maxMVV} --minMX {minMX} --maxMX {maxMX} {samples} '.format(template=template,cut=cut,rootFile=rootFile,minMVV=self.minMVV,maxMVV=self.maxMVV,minMX=self.minMX,maxMX=self.maxMX,fixPars=fixPars,samples=self.samples)
  os.system(cmd)
  
  jsonFile=filename+"_MVV.json"
  cmd='vvMakeJSON.py  -o "{jsonFile}" -g "MEAN:pol1,SIGMA:pol6,ALPHA1:pol5,N1:pol0,ALPHA2:pol4,N2:pol0" -m {minMX} -M {maxMX} {rootFile}  '.format(jsonFile=jsonFile,rootFile=rootFile,minMX=self.minMX,maxMX=self.maxMX)
  os.system(cmd)
  print "*** makeSignalShapesMVV start ****" 

 def makeSignalShapesMJ(self,filename,template,leg):
  print "*** makeSignalShapesMJ start ****" 
  for c in self.categories:
  
   cut='*'.join([self.cuts['common'],self.cuts[c]])
     
   #do not fix fit parameters for the moment --> to be optimized
   fixPars="" 
   rootFile=filename+"_MJ"+leg+"_"+c+".root"     
   cmd='vvMakeSignalMJShapes.py -s "{template}" -c "{cut}"  -o "{rootFile}" -V "jj_{leg}_softDrop_mass" -m {minMJ} -M {maxMJ} -f "{fixPars}" --minMX {minMX} --maxMX {maxMX} {samples} '.format(template=template,cut=cut,rootFile=rootFile,leg=leg,minMJ=self.minMJ,maxMJ=self.maxMJ,minMX=self.minMX,maxMX=self.maxMX,fixPars=fixPars,samples=self.samples)
   os.system(cmd)
   
   jsonFile=filename+"_MJ"+leg+"_"+c+".json"   
   cmd='vvMakeJSON.py  -o "{jsonFile}" -g "mean:pol3,sigma:pol3,alpha:pol3,n:pol3,alpha2:pol3,n2:pol3,slope:pol0,f:pol0,meanH:pol0,sigmaH:pol0,alphaH:pol0,nH:pol0,alpha2H:pol0,n2H:pol0,slopeH:pol0,fH:pol0" -m {minMX} -M {maxMX} {rootFile}  '.format(jsonFile=jsonFile,rootFile=rootFile,minMX=self.minMX,maxMX=self.maxMX)
   if filename.find('WH') != -1 or filename.find('ZH') != -1: cmd='vvMakeJSON.py  -o "{jsonFile}" -g "mean:pol3,sigma:pol3,alpha:pol3,n:pol3,alpha2:pol3,n2:pol3,slope:pol0,f:pol0,meanH:pol3,sigmaH:pol3,alphaH:pol3,nH:pol3,alpha2H:pol3,n2H:pol3,slopeH:pol0,fH:pol0" -m {minMX} -M {maxMX} {rootFile}  '.format(jsonFile=jsonFile,rootFile=rootFile,minMX=self.minMX,maxMX=self.maxMX)   
   os.system(cmd)
  print "*** makeSignalShapesMJ end ****" 

 def makeSignalYields(self,filename,template,branchingFraction,sfP = {'HPHP':1.0,'HPLP':1.0,'LPLP':1.0}):
  print "*** makeSignalYelds start ****" 
  print "using the following scalfactors:" ,sfP
  
  for c in self.categories:
   cut = "*".join([self.cuts[c],self.cuts['common'],self.cuts['acceptance'],str(sfP[c])])
   yieldFile=filename+"_"+c+"_yield"
   fnc = "pol7"
   cmd='vvMakeSignalYields.py -s {template} -c "{cut}" -o {output} -V "jj_LV_mass" -m {minMVV} -M {maxMVV} -f {fnc} -b {BR} --minMX {minMX} --maxMX {maxMX} {samples} '.format(template=template, cut=cut, output=yieldFile,minMVV=self.minMVV,maxMVV=self.maxMVV,fnc=fnc,BR=branchingFraction,minMX=self.minMX,maxMX=self.maxMX,samples=self.samples)
   os.system(cmd)
  print "*** makeSignalYelds end ****" 

 def makeDetectorResponse(self,name,filename,template,addCut="1",jobName="DetPar"):
   print "*** makeDetectorResponse start ****" 
   cut='*'.join([self.cuts['common'],addCut,self.cuts['acceptanceGEN'],self.cuts['looseacceptanceMJ']])
   resFile=filename+"_"+name+"_detectorResponse.root"	    
   print "Saving detector resolution to file: " ,resFile
   bins = "200,250,300,350,400,450,500,600,700,800,900,1000,1200,1500,1800,2200,2600,3000,3400,3800,5000"#4200,4600,5000,7000"#TODO: The last three bins are empty, remove next iteration!
   if self.submitToBatch:
    from modules.submitJobs import Make2DDetectorParam,merge2DDetectorParam 
    jobList, files = Make2DDetectorParam(resFile,template,cut,samples,jobName,bins)
    jobList = []
    files = []
    merge2DDetectorParam(resFile,bins,jobName)
   else:
    cmd='vvMake2DDetectorParam.py  -o "{rootFile}" -s "{template}" -c "{cut}"  -v "jj_LV_mass,jj_l1_softDrop_mass"  -g "jj_gen_partialMass,jj_l1_gen_softDrop_mass,jj_l1_gen_pt"  -b {bins}   {samples}'.format(rootFile=resFile,template=template,cut=cut,minMVV=self.minMVV,maxMVV=self.maxMVV,tag=name,bins=bins,samples=self.samples)
    os.system(cmd)
   
   print "Done with ",resFile
   print "*** makeDetectorResponse end ****" 

 def makeBackgroundShapesMVVKernel(self,name,filename,template,addCut="1",jobName="1DMVV",wait=True,corrFactorW=1,corrFactorZ=1):
  print "*** makeBackgroundShapesMVVKernel start ****" 
  pwd = os.getcwd()
  
  for c in self.categories:
  
   jobname = jobName+"_"+c
   print "Working on purity: ", c
   
   resFile  = pwd + "/"+ filename+"_"+name+"_detectorResponse.root"
   print "Reading " ,resFile

   rootFile = filename+"_"+name+"_MVV_"+c+".root"
   print "Saving to ",rootFile
   
   cut='*'.join([self.cuts['common'],self.cuts[c],addCut,self.cuts['acceptanceGEN'],self.cuts['looseacceptanceMJ']])
   smp = pwd +"/"+self.samples

   if self.submitToBatch:
    if name.find("Jets") == -1: template += ",QCD_Pt-,QCD_HT"
    from modules.submitJobs import Make1DMVVTemplateWithKernels,merge1DMVVTemplate
    jobList, files = Make1DMVVTemplateWithKernels(rootFile,template,cut,resFile,self.binsMVV,self.minMVV,self.maxMVV,smp,jobname,wait,self.HCALbinsMVV,addOption)
    if wait: merge1DMVVTemplate(jobList,files,jobname,c,self.binsMVV,self.minMVV,self.maxMVV,self.HCALbinsMVV)
   else:
    cmd='vvMake1DMVVTemplateWithKernels.py -H "x" -o "{rootFile}" -s "{template}" -c "{cut}"  -v "jj_gen_partialMass" -b {binsMVV}  -x {minMVV} -X {maxMVV} -r {res} {directory} --corrFactorW {corrFactorW} --corrFactorZ {corrFactorZ} '.format(rootFile=rootFile,template=template,cut=cut,res=resFile,binsMVV=self.binsMVV,minMVV=self.minMVV,maxMVV=self.maxMVV,corrFactorW=corrFactorW,corrFactorZ=corrFactorZ,directory=smp)
    cmd = cmd+self.HCALbinsMVV
    os.system(cmd)    
  print "*** makeBackgroundShapesMVVKernel end ****"
 
 def makeBackgroundShapesMVVConditional(self,name,filename,template,leg,addCut="",jobName="2DMVV",wait=True):
  print "*** makeBackgroundShapesMVVConditional start ****"  
  pwd = os.getcwd()  
  
  for c in self.categories:

   print " Working on purity: ", c

   jobname = jobName+"_"+c
   resFile=filename+"_"+name+"_detectorResponse.root"
   rootFile=filename+"_"+name+"_COND2D_"+c+"_"+leg+".root"       
   print "Reading " ,resFile
   print "Saving to ",rootFile
   
   cut='*'.join([self.cuts['common'],self.cuts[c],addCut])#,cuts['acceptanceGEN'],cuts['looseacceptanceMJ']])
   
   smp = pwd +"/"+self.samples 
  
   if self.submitToBatch:
    if name.find("VJets")== -1: template += ",QCD_Pt-,QCD_HT"
    from modules.submitJobs import Make2DTemplateWithKernels,merge2DTemplate
    jobList, files = Make2DTemplateWithKernels(rootFile,template,cut,leg,binsMVV,minMVV,maxMVV,resFile,binsMJ,minMJ,maxMJ,smp,jobname,wait,HCALbinsMVV,addOption)
    if wait: merge2DTemplate(jobList,files,jobname,p,leg,binsMVV,binsMJ,minMVV,maxMVV,minMJ,maxMJ,HCALbinsMVV)
   else:
      cmd='vvMake2DTemplateWithKernels.py -o "{rootFile}" -s "{template}" -c "{cut}"  -v "jj_{leg}_gen_softDrop_mass,jj_gen_partialMass"  -b {binsMJ} -B {binsMVV} -x {minMJ} -X {maxMJ} -y {minMVV} -Y {maxMVV}  -r {res} {samples}'.format(rootFile=rootFile,template=template,cut=cut,leg=leg,binsMVV=self.binsMVV,minMVV=self.minMVV,maxMVV=self.maxMVV,res=resFile,binsMJ=self.binsMJ,minMJ=self.minMJ,maxMJ=self.maxMJ,samples=smp)
      cmd=cmd+self.HCALbinsMVV
      os.system(cmd)
  print "*** makeBackgroundShapesMVVConditional end ****"  

 def mergeBackgroundShapes(self,name,filename):
  print "*** mergeBackgroundShapes start ****"  
  for c in self.categories:
  
   inputx=filename+"_"+name+"_COND2D_"+c+"_l1.root"
   inputy=filename+"_"+name+"_COND2D_"+c+"_l2.root"
   inputz=filename+"_"+name+"_MVV_"+c+".root"     
   rootFile=filename+"_"+name+"_3D_"+c+".root"
   
   print "Reading " ,inputx
   print "Reading " ,inputy
   print "Reading " ,inputz
   print "Saving to ",rootFile 
   
   cmd='vvMergeHistosToPDF3D.py -i "{inputx}" -I "{inputy}" -z "{inputz}" -o "{rootFile}"'.format(rootFile=rootFile,inputx=inputx,inputy=inputy,inputz=inputz)
   os.system(cmd)
   
   #print "Adding trigger shape uncertainties"
   #if useTriggerWeights: 
   #   cmd='vvMakeTriggerShapes.py -i "{rootFile}"'.format(rootFile=rootFile)
   #   os.system(cmd)
   print "*** mergeBackgroundShapes end ****"
  
 def makeNormalizations(self,name,filename,template,data=0,addCut='1',jobName="nR",factors="1",wait=True):
  print "*** makeNormalizations start ****"  
  pwd = os.getcwd()
  sam = pwd +"/"+self.samples
  print "Using files in" , sam
  
  for c in self.categories:
   
   #apply V/H tagging scale factors --> this will have to be updated
   if name.find("Jets")!=-1:
        if c.find("HPLP"): factors=factors+",sf:"+str(LPSF)
        else: factors=factors+",sf:"+str(HPSF)
      
   jobname = jobName+"_"+c
   rootFile=filename+"_"+name+"_"+c+".root"
   print "Saving to ",rootFile  
   
   cut='*'.join([self.cuts['common'],self.cuts[c],addCut,self.cuts['acceptance']])

   if self.submitToBatch:
       if name.find("nonRes")!= -1: template += ",QCD_Pt-,QCD_HT"
       from modules.submitJobs import makeData,mergeData
       jobList, files = makeData(template,cut,rootFile,binsMVV,binsMJ,minMVV,maxMVV,minMJ,maxMJ,factors,name,data,jobname,sam,wait,HCALbinsMVV,addOption)
       wait = True
       mergeData(jobname,p,rootFile)
   else:
        cmd='vvMakeData.py -s "{template}" -d {data} -c "{cut}"  -o "{rootFile}" -v "jj_l1_softDrop_mass,jj_l2_softDrop_mass,jj_LV_mass" -b "{bins},{bins},{BINS}" -m "{mini},{mini},{MINI}" -M "{maxi},{maxi},{MAXI}" -f {factors} -n "{name}" {samples}'.format(template=template,cut=cut,rootFile=rootFile,BINS=self.binsMVV,bins=self.binsMJ,MINI=self.minMVV,MAXI=self.maxMVV,mini=self.minMJ,maxi=self.maxMJ,factors=factors,name=name,data=data,samples=sam)
        cmd=cmd+self.HCALbinsMVV
        os.system(cmd)
  print "*** makeNormalizations end ****"  
 #this one I still have to fix and test, do not use submitToBatch yet
 def mergeKernelJobs(self):
    
    for p in purities:
        jobList = []
        files   = []
        with open("tmp1D_%s_joblist.txt"%p,'r') as infile:
            for line in infile:
                if line.startswith("job"):
                    for job in line.split("[")[1].split("]")[0].split(","):
                        jobList.append(job.replace("'","").replace(" ",""))
                if line.startswith("file"):
                    for job in line.split("[")[1].split("]")[0].split(","):
                        files.append(job.replace("'","").replace(" ",""))
        from modules.submitJobs import merge1DMVVTemplate
        merge1DMVVTemplate(jobList,files,"1D"+"_"+p,p,binsMVV,binsMJ,minMVV,maxMVV,minMJ,maxMJ,HCALbinsMVV)
      
        jobList = []
        files   = []
        with open("tmp2Dl1_%s_joblist.txt"%p,'r') as infile:
            for line in infile:
                if line.startswith("job"):
                    for job in line.split("[")[1].split("]")[0].split(","):
                        jobList.append(job.replace("'","").replace(" ",""))
                if line.startswith("file"):
                    for job in line.split("[")[1].split("]")[0].split(","):
                        files.append(job.replace("'","").replace(" ",""))

        from modules.submitJobs import merge2DTemplate
        merge2DTemplate(jobList,files,"2Dl1"+"_"+p,p,"l1",binsMVV,binsMJ,minMVV,maxMVV,minMJ,maxMJ,HCALbinsMVV)
        merge2DTemplate(jobList,files,"2Dl2"+"_"+p,p,"l2",binsMVV,binsMJ,minMVV,maxMVV,minMJ,maxMJ,HCALbinsMVV)
		            	    
 def printAllParameters(self):
 
  print "--------------------------------------------------------------------------"  
  print "Cuts:"
  print "--------------------------------------------------------------------------"
  print self.cuts
  print "--------------------------------------------------------------------------"
  print "mjj ranges: min ",self.minMVV,"max",self.maxMVV, "nbins",self.binsMVV
  print "--------------------------------------------------------------------------"
  print "reonance mX ranges: min ",self.minMX,"max",self.maxMX
  print "--------------------------------------------------------------------------" 
  print "mjj binning for bkg/data templates"
  print "--------------------------------------------------------------------------" 
  print self.HCALbinsMVV
  print "--------------------------------------------------------------------------" 
  print "Samples directory:",self.samples
  print "--------------------------------------------------------------------------" 
  print "Categories:",self.categories
  print "--------------------------------------------------------------------------" 
  print "mjet ranges: min ",self.minMJ,"max",self.maxMJ,"nbins",self.binsMJ
  print "--------------------------------------------------------------------------" 
  print "Submit to condor batch:",self.submitToBatch
  print "--------------------------------------------------------------------------" 
