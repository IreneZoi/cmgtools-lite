#!/usr/bin/env python
import os, re, copy
import commands
import math, time
import sys
import ROOT
from ROOT import *
import subprocess, thread
from array import array
ROOT.gROOT.SetBatch(True)

timeCheck = "30"
userName=os.environ['USER']


def getBinning(binsMVV):
    l=[]
    if binsMVV=="":
        return l
    else:
        s = binsMVV.split(",")
        for w in s:
            l.append(int(w))
    return l

useCondorBatch = True

def makeSubmitFileCondor(exe,jobname,jobflavour):
    print "make options file for condor job submission "
    submitfile = open("submit.sub","w")
    submitfile.write("executable  = "+exe+"\n")
    submitfile.write("arguments             = $(ClusterID) $(ProcId)\n")
    submitfile.write("output                = "+jobname+".$(ClusterId).$(ProcId).out\n")
    submitfile.write("error                 = "+jobname+".$(ClusterId).$(ProcId).err\n")
    submitfile.write("log                   = "+jobname+".$(ClusterId).log\n")
    submitfile.write('+JobFlavour           = "'+jobflavour+'"\n')
    submitfile.write("queue")
    submitfile.close()
    
def waitForBatchJobs( jobname, remainingjobs, listOfJobs, userName, timeCheck="30"):
    if listOfJobs-remainingjobs < listOfJobs:
        time.sleep(float(timeCheck))
        # nprocess = "bjobs -u %s | awk {'print $9'} | grep %s | wc -l" %(userName,jobname)
        if useCondorBatch:
                nprocess = "condor_q %s | grep %s | wc -l"%(userName,jobname)
        else:
                nprocess = "bjobs -u %s | grep %s | wc -l" %(userName,jobname)
        result = subprocess.Popen(nprocess, stdout=subprocess.PIPE, shell=True)
        runningJobs =  int(result.stdout.read())
        print "waiting for %d job(s) in the queue (out of total %d)" %(runningJobs,listOfJobs)
        waitForBatchJobs( jobname, runningJobs, listOfJobs, userName, timeCheck)  
    else:
        print "Jobs finished! Allow some time for files to be saved to your home directory"
        time.sleep(5)
        noutcmd = "ls res"+jobname+"/ | wc -l"
        result = subprocess.Popen(noutcmd, stdout=subprocess.PIPE, shell=True)
        nout =  int(result.stdout.read())
        if listOfJobs-nout > 0: 
            print "Uooh! Missing %i jobs. Resubmit when merging. Now return to main script" %int(listOfJobs-nout)
        else:
            print "Done! Have %i out of %i files in res%s directory. Return to main script"%(nout,listOfJobs,jobname)
        return
        
def submitJobs(minEv,maxEv,cmd,OutputFileNames,queue,jobname,path):
    joblist = []
    for k,v in minEv.iteritems():

      for j in range(len(v)):
       os.system("mkdir tmp"+jobname+"/"+str(k).replace(".root","")+"_"+str(j+1))
       os.chdir("tmp"+jobname+"/"+str(k).replace(".root","")+"_"+str(j+1))
      
       with open('job_%s_%i.sh'%(k.replace(".root",""),j+1), 'w') as fout:
          fout.write("#!/bin/sh\n")
          fout.write("echo\n")
          fout.write("echo\n")
          fout.write("echo 'START---------------'\n")
          fout.write("echo 'WORKDIR ' ${PWD}\n")
          fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
          fout.write("cd "+str(path)+"\n")
          fout.write("cmsenv\n")
          fout.write(cmd+" -o res"+jobname+"/"+OutputFileNames+"_"+str(j+1)+"_"+k+" -s "+k+" -e "+str(minEv[k][j])+" -E "+str(maxEv[k][j])+"\n")
          fout.write("echo 'STOP---------------'\n")
          fout.write("echo\n")
          fout.write("echo\n")

       if useCondorBatch:
               os.system("mv  job_*.sh "+jobname+".sh")
               makeSubmitFileCondor(jobname+".sh",jobname,"workday")
               os.system("condor_submit submit.sub")
       else:
               os.system("chmod 755 job_%s_%i.sh"%(k.replace(".root",""),j+1) )
               os.system("bsub -q "+queue+" -o logs job_%s_%i.sh -J %s"%(k.replace(".root",""),j+1,jobname))
       print "job nr " + str(j+1) + " file " + k + " being submitted"
       joblist.append("%s_%i"%(k.replace(".root",""),j+1))
       os.chdir("../..")
    return joblist     

def getEvents(template,samples):
    files = []
    sampleTypes = template.split(',')
    for f in os.listdir(samples):
        for t in sampleTypes:
            if f.find('.root') != -1 and f.find(t) != -1: files.append(f)

    minEv = {}
    maxEv = {}

    for f in files:
     inf = ROOT.TFile('%s/'%samples+f,'READ')
     print "opening file "+str(samples+f)
     minEv[f] = []
     maxEv[f] = []
     intree = inf.Get('AnalysisTree')
     nentries = intree.GetEntries()
     print f,nentries,nentries/500000
     if nentries/500000 == 0:
      minEv[f].append(0)
      maxEv[f].append(500000)
     else:
      for i in range(nentries/500000):
       minEv[f].append(i*500000)
       maxEv[f].append(499999)


    NumberOfJobs = 0
    for k in maxEv.keys():
     NumberOfJobs += len(maxEv[k])

    return minEv, maxEv, NumberOfJobs, files
    
def Make2DDetectorParam(rootFile,template,cut,samples,jobname="DetPar",bins="200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000",wait=True,leg="l1"):
   
    print 
    print 'START: Make2DDetectorParam with parameters:'
    print
    print "rootFile = %s" %rootFile  
    print "template = %s" %template  
    print "cut      = %s" %cut       
    print "samples  = %s" %samples   
    print "jobname  = %s" %jobname 
    
    if leg == "l1": cmd='vvMake2DDetectorParam.py  -c "{cut}"  -v "jj_LV_mass,jj_l1_softDrop_mass"  -g "jj_gen_partialMass,jj_l1_gen_softDrop_mass,jj_l1_gen_pt"  -b {bins}   {infolder}'.format(rootFile=rootFile,samples=template,cut=cut,bins=bins,infolder=samples)
    if leg == "l2": cmd='vvMake2DDetectorParam.py  -c "{cut}"  -v "jj_LV_mass,jj_l2_softDrop_mass"  -g "jj_gen_partialMass,jj_l2_gen_softDrop_mass,jj_l2_gen_pt"  -b {bins}   {infolder}'.format(rootFile=rootFile,samples=template,cut=cut,bins=bins,infolder=samples)
    OutputFileNames = rootFile.replace(".root","") # base of the output file name, they will be saved in res directory
    queue = "8nh" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 
    
    files = []
    sampleTypes = template.split(',')
    for f in os.listdir(samples):
        for t in sampleTypes:
            if f.find('.root') != -1 and f.find(t) != -1:
                files.append(f)
                print "file "+str(f)+" appended to jobs"
 
    NumberOfJobs= len(files) 
    print
    print "Submitting %i number of jobs "  ,NumberOfJobs
    print
    
    path = os.getcwd()
    try: os.system("rm -r tmp"+jobname)
    except: print "No tmp/ directory"
    os.system("mkdir tmp"+jobname)
    try: os.stat("res"+jobname) 
    except: os.mkdir("res"+jobname)

    #### Creating and sending jobs #####
    joblist = []
    ##### loop for creating and sending jobs #####
    for x in range(1, int(NumberOfJobs)+1):
     
       os.system("mkdir tmp"+jobname+"/"+str(files[x-1]).replace(".root",""))
       os.chdir("tmp"+jobname+"/"+str(files[x-1]).replace(".root",""))
       with open('job_%s.sh'%files[x-1].replace(".root",""), 'w') as fout:
          fout.write("#!/bin/sh\n")
          fout.write("echo\n")
          fout.write("echo\n")
          fout.write("echo 'START---------------'\n")
          fout.write("echo 'WORKDIR ' ${PWD}\n")
          fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
          fout.write("cd "+str(path)+"\n")
          fout.write("cmsenv\n")
          fout.write(cmd+" -o "+path+"/res"+jobname+"/"+OutputFileNames+"_"+files[x-1]+" -s "+files[x-1]+"\n")
          fout.write("echo 'STOP---------------'\n")
          fout.write("echo\n")
          fout.write("echo\n")
       if useCondorBatch:
         os.system("mv  job_*.sh "+jobname+".sh")
         makeSubmitFileCondor(jobname+".sh",jobname,"workday")
         os.system("condor_submit submit.sub")
       else:
         os.system("chmod 755 job_%s.sh"%(files[x-1].replace(".root","")) )
         os.system("bsub -q "+queue+" -o logs job_%s.sh -J %s"%(files[x-1].replace(".root",""),jobname))
       print "job nr " + str(x) + " submitted"
       joblist.append("%s"%(files[x-1].replace(".root","")))
       os.chdir("../..")
   
    print
    print "your jobs:"
    if useCondorBatch:
            os.system("condor_q")
    else:
            os.system("bjobs")
    userName=os.environ['USER']
    if wait: waitForBatchJobs(jobname,NumberOfJobs,NumberOfJobs, userName, timeCheck)
    
    print
    print 'END: Make2DDetectorParam'
    print
    return joblist, files   
    
def Make1DMVVTemplateWithKernels(rootFile,template,cut,resFile,binsMVV,minMVV,maxMVV,samples,jobName="1DMVV",wait=True,binning='',addOption=""):

    
    print 
    print 'START: Make1DMVVTemplateWithKernels with parameters:'
    print
    print "rootFile = %s" %rootFile  
    print "template = %s" %template  
    print "cut      = %s" %cut       
    print "resFile  = %s" %resFile   
    print "binsMVV  = %i" %binsMVV   
    print "minMVV   = %i" %minMVV    
    print "maxMVV   = %i" %maxMVV    
    print "samples  = %s" %samples   
    print "jobName  = %s" %jobName 
    print
    minEv, maxEv, NumberOfJobs, files = getEvents(template,samples) 
    print "Submitting %i number of jobs "  ,NumberOfJobs
    print

    cmd='vvMake1DMVVTemplateWithKernels.py -H "x" -c "{cut}"  -v "jj_gen_partialMass" {binning} -b {binsMVV}  -x {minMVV} -X {maxMVV} -r {res} {addOption} {infolder} '.format(rootFile=rootFile,cut=cut,res=resFile,binsMVV=binsMVV,minMVV=minMVV,maxMVV=maxMVV,infolder=samples,binning=binning,addOption=addOption)
    OutputFileNames = rootFile.replace(".root","") # base of the output file name, they will be saved in res directory
    queue = "8nh" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 
    
    path = os.getcwd()
    try: os.system("rm -r tmp"+jobName)
    except: print "No tmp/ directory"
    os.system("mkdir tmp"+jobName)
    try: os.stat("res"+jobName) 
    except: os.mkdir("res"+jobName)

    #### Creating and sending jobs #####
    joblist = submitJobs(minEv,maxEv,cmd,OutputFileNames,queue,jobName,path)
    with open('tmp'+jobName+'_joblist.txt','w') as outfile:
        outfile.write("jobList = %s\n" % joblist)
        outfile.write("files = %s\n" % files)
    outfile.close()
    print
    print "your jobs:"
    if useCondorBatch:
        os.system("condor_q")
    else:
        os.system("bjobs")
    if wait: waitForBatchJobs(jobName,NumberOfJobs,NumberOfJobs, userName, timeCheck)
    print
    print 'END: Make1DMVVTemplateWithKernels'
    print 

    return joblist, files 


def Make2DTemplateWithKernels(rootFile,template,cut,leg,binsMVV,minMVV,maxMVV,resFile,binsMJ,minMJ,maxMJ,samples,jobName="2DMVV",wait=True,binning='',addOption=''):

    
    print 
    print 'START: Make2DTemplateWithKernels'
    print
    print "rootFile  = %s" %rootFile   
    print "template  = %s" %template   
    print "cut       = %s" %cut       
    print "leg       = %s" %leg       
    print "binsMVV   = %i" %binsMVV    
    print "minMVV    = %i" %minMVV     
    print "maxMVV    = %i" %maxMVV     
    print "resFile   = %s" %resFile    
    print "binsMJ    = %s" %binsMJ    
    print "minMJ     = %s" %minMJ 
    print "maxMJ     = %s" %maxMJ 
    print "samples   = %s" %samples   
    print "jobName   = %s" %jobName   
    print
    minEv, maxEv, NumberOfJobs, files = getEvents(template,samples) 
    print "Submitting %i number of jobs "  ,NumberOfJobs
    print

    cmd='vvMake2DTemplateWithKernels.py -c "{cut}"  -v "jj_{leg}_gen_softDrop_mass,jj_gen_partialMass" {binning}  -b {binsMJ} -B {binsMVV} -x {minMJ} -X {maxMJ} -y {minMVV} -Y {maxMVV}  -r {res} {addOption} {infolder}'.format(rootFile=rootFile,samples=template,cut=cut,leg=leg,binsMVV=binsMVV,minMVV=minMVV,maxMVV=maxMVV,res=resFile,binsMJ=binsMJ,minMJ=minMJ,maxMJ=maxMJ,infolder=samples,binning=binning,addOption=addOption)
    OutputFileNames = rootFile.replace(".root","") # base of the output file name, they will be saved in res directory
    queue = "8nh" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 
    
    path = os.getcwd()
    try: os.system("rm -r tmp"+jobName)
    except: print "No tmp/ directory"
    os.system("mkdir tmp"+jobName)
    try: os.stat("res"+jobName) 
    except: os.mkdir("res"+jobName)

    #### Creating and sending jobs #####
    joblist = submitJobs(minEv,maxEv,cmd,OutputFileNames,queue,jobName,path)
    with open('tmp'+jobName+'_joblist.txt','w') as outfile:
        outfile.write("jobList = %s\n" % joblist)
        outfile.write("files = %s\n" % files)
    outfile.close()
    print
    print "your jobs:"
    if useCondorBatch:
        os.system("condor_q")
    else:
        os.system("bjobs")
    userName=os.environ['USER']
    if wait: waitForBatchJobs(jobName,NumberOfJobs,NumberOfJobs, userName, timeCheck)
    
    
      
    print
    print 'END: Make2DTemplateWithKernels'
    print

    return joblist, files

def unequalScale(histo,name,alpha,power=1,dim=1):
    newHistoU =copy.deepcopy(histo) 
    newHistoU.SetName(name+"Up")
    newHistoD =copy.deepcopy(histo) 
    newHistoD.SetName(name+"Down")
    if dim == 2:
	    maxFactor = max(pow(histo.GetXaxis().GetXmax(),power),pow(histo.GetXaxis().GetXmin(),power))
	    for i in range(1,histo.GetNbinsX()+1):
	        x= histo.GetXaxis().GetBinCenter(i)
	        for j in range(1,histo.GetNbinsY()+1):
	            nominal=histo.GetBinContent(i,j)
	            factor = 1+alpha*pow(x,power) 
	            newHistoU.SetBinContent(i,j,nominal*factor)
	            newHistoD.SetBinContent(i,j,nominal/factor)
	    if newHistoU.Integral()>0.0:        
	        newHistoU.Scale(1.0/newHistoU.Integral())        
	    if newHistoD.Integral()>0.0:        
	        newHistoD.Scale(1.0/newHistoD.Integral())        
    else:
	    for i in range(1,histo.GetNbinsX()+1):
	        x= histo.GetXaxis().GetBinCenter(i)
	        nominal=histo.GetBinContent(i) #ROOT.TMath.Log10(histo.GetBinContent(i))
		factor = 1+alpha*pow(x,power)
		#print i,x,power,alpha,factor,nominal,nominal*factor,nominal/factor
	        newHistoU.SetBinContent(i,nominal*factor)
	        if factor != 0: newHistoD.SetBinContent(i,nominal/factor)	
    return newHistoU,newHistoD 
    
def mirror(histo,histoNominal,name,dim=1):
    newHisto =copy.deepcopy(histoNominal) 
    newHisto.SetName(name)
    intNominal=histoNominal.Integral()
    intUp = histo.Integral()
    print intNominal,intUp
    if dim == 2:
		for i in range(1,histo.GetNbinsX()+1):
			for j in range(1,histo.GetNbinsY()+1):
				up=histo.GetBinContent(i,j)/intUp
				nominal=histoNominal.GetBinContent(i,j)/intNominal
				if up != 0: newHisto.SetBinContent(i,j,histoNominal.GetBinContent(i,j)*nominal/up)
    else:
		for i in range(1,histo.GetNbinsX()+1):
			up=histo.GetBinContent(i)/intUp
			nominal=histoNominal.GetBinContent(i)/intNominal
			newHisto.SetBinContent(i,histoNominal.GetBinContent(i)*nominal/up)	
    return newHisto       

def expandHisto(histo,suffix,binsMVV,binsMJ,minMVV,maxMVV,minMJ,maxMJ):
    histogram=ROOT.TH2F(histo.GetName()+suffix,"histo",binsMJ,minMJ,maxMJ,binsMVV,minMVV,maxMVV)
    for i in range(1,histo.GetNbinsX()+1):
        proje = histo.ProjectionY("q",i,i)
        graph=ROOT.TGraph(proje)
        for j in range(1,histogram.GetNbinsY()+1):
            x=histogram.GetYaxis().GetBinCenter(j)
            bin=histogram.GetBin(i,j)
            histogram.SetBinContent(bin,graph.Eval(x,0,"S"))
    return histogram

def expandHistoBinned(histo,suffix ,binsx,binsy):
    histogram=ROOT.TH2F(histo.GetName()+suffix,"histo",len(binsx)-1,array('f',binsx),len(binsy)-1,array('f',binsy))
    for i in range(1,histo.GetNbinsX()+1):
        proje = histo.ProjectionY("q",i,i)
        graph=ROOT.TGraph(proje)
        for j in range(1,histogram.GetNbinsY()+1):
            x=histogram.GetYaxis().GetBinCenter(j)
            bin=histogram.GetBin(i,j)
            histogram.SetBinContent(bin,graph.Eval(x,0,"S"))
    return histogram
        
def conditional(hist):
    for i in range(1,hist.GetNbinsY()+1):
        proj=hist.ProjectionX("q",i,i)
        integral=proj.Integral()
        if integral==0.0:
            print 'SLICE WITH NO EVENTS!!!!!!!!',hist.GetName()
            continue
        for j in range(1,hist.GetNbinsX()+1):
            hist.SetBinContent(j,i,hist.GetBinContent(j,i)/integral)

def getJobs(files,jobList,outdir,purity):
	resubmit = []
	jobsPerSample = {}
	exit_flag = False

	for s in files:
	 s = s.replace('.root','')
	 filelist = []
	 for t in jobList:
	  if t.find(s) == -1: continue
	  jobid = t.split("_")[-1]
	  found = False
	  for o in os.listdir(outdir):
	   if o.find(s) != -1 and o.find('_'+jobid+'_') != -1 and o.find(purity)!=-1:
	    found = True
	    filelist.append(outdir+"/"+o)
	    break
	  if not found:
	   print "SAMPLE ",s," JOBID ",jobid," NOT FOUND"
	   exit_flag = True
	   resubmit.append(s+"_"+jobid)
	 if len(filelist) > 0: jobsPerSample[s] = filelist
	return resubmit, jobsPerSample,exit_flag	 
	
def reSubmit(jobdir,resubmit,jobname):
 jobs = []
 for o in os.listdir(jobdir):
	 for jobs in resubmit:
		 if o.find(jobs) != -1: 
			 jobfolder = jobdir+"/"+jobs+"/"
			 os.chdir(jobfolder)
			 if useCondorBatch:
			    cmd = "condor_submit submit.sub"
			    script = jobname+".sh"
                         else:
                            script = "job_"+jobs+".sh"
                            cmd = "bsub -q 8nh -o logs %s -J %s"%(script,jobname)
			 print cmd
			 jobs += cmd
			 os.system("chmod 755 %s"%script)
			 os.system(cmd)
			 os.chdir("../..")
 return jobs

def merge2DDetectorParam(resFile,binsxStr,jobname,template="QCD_Pt_"):
    
    print "Merging 2D detector parametrization"

    outdir = 'res'+jobname
    jobdir = 'tmp'+jobname
    print "jobname "+str(jobname)
    print "outdir "+str(outdir)
    print "jobdir "+str(jobdir)
   
    filelist = os.listdir('./res'+jobname+'/')
    print "filelist "+str(filelist)

    pythia_files = []
    herwig_files = []
    mg_files = []
    pythia_template="QCD_Pt_"
    herwig_template="QCD_Pt-"
    mg_template="QCD_HT_"
    for f in filelist:
     if f.find(resFile.replace(".root",""))==-1: continue
     if f.find(pythia_template) != -1: pythia_files.append('./res'+jobname+'/'+f)
     elif f.find(mg_template) != -1: mg_files.append('./res'+jobname+'/'+f)
     elif f.find(herwig_template) != -1: herwig_files.append('./res'+jobname+'/'+f)
    

    #now hadd them
    tmp_files = []
    if len(pythia_files) > 0:
        print "doing pythia"
        cmd = 'hadd -f tmp_nominal.root '
        for f in pythia_files:
         cmd += f
         cmd += ' '
        print cmd
        os.system(cmd)
        tmp_files.append('tmp_nominal.root')
        
    if len(mg_files) > 0:
        print "doing mg"
        cmd = 'hadd -f tmp_altshape2.root '
        for f in mg_files:
         cmd += f
         cmd += ' '
        print cmd
        os.system(cmd)  
        tmp_files.append('tmp_altshape2.root')  

    if len(herwig_files) > 0:
        print "doing herwig"
        cmd = 'hadd -f tmp_altshapeUp.root '
        for f in herwig_files:
         cmd += f
         cmd += ' '
        print cmd
        os.system(cmd)  
        tmp_files.append('tmp_altshapeUp.root')
        
    #produce final det resolution files (one per sample, but at the end we use the pythia one in the following steps for all the samples)
    for f in tmp_files:
     print " file in tmp_files ",f 
    
     fin = ROOT.TFile.Open(f,'READ')
     print "WORKING ON FILE ", fin.GetName()
     label = fin.GetName().split('_')[1].split('.')[0]
     print "label ", label
     superHX = fin.Get("dataX")
     superHY = fin.Get("dataY")
     # superHNsubj = fin.Get("dataNsubj")
     
     binsx=[]
     for b in binsxStr.split(','):
         binsx.append(float(b))
     outname = resFile.replace(".root","")+f.split('_')[1]
     print "outname "+str(outname)
     fout = ROOT.TFile(outname,"RECREATE")
     
     scalexHisto=ROOT.TH1F("scalexHisto","scaleHisto",len(binsx)-1,array('d',binsx))
     resxHisto=ROOT.TH1F("resxHisto","resHisto",len(binsx)-1,array('d',binsx))
     scaleyHisto=ROOT.TH1F("scaleyHisto","scaleHisto",len(binsx)-1,array('d',binsx))
     resyHisto=ROOT.TH1F("resyHisto","resHisto",len(binsx)-1,array('d',binsx))
     #scaleNsubjHisto=ROOT.TH1F("scaleNsubjHisto","scaleHisto",len(binsx)-1,array('d',binsx))
     #resNsubjHisto=ROOT.TH1F("resNsubjHisto","resHisto",len(binsx)-1,array('d',binsx))
     
     for bin in range(1,superHX.GetNbinsX()+1):
        tmp=superHX.ProjectionY("q",bin,bin)
        if bin==1: 
                 scalexHisto.SetBinContent(bin,tmp.GetMean())
                 scalexHisto.SetBinError(bin,tmp.GetMeanError())
                 resxHisto.SetBinContent(bin,tmp.GetRMS())
                 resxHisto.SetBinError(bin,tmp.GetRMSError())
                 continue        
        startbin   = 0.
        maxcontent = 0.
        for b in range(tmp.GetXaxis().GetNbins()):
          if tmp.GetXaxis().GetBinCenter(b+1) > startbin and tmp.GetBinContent(b+1)>maxcontent:
            maxbin = b
            maxcontent = tmp.GetBinContent(b+1)
        tmpmean = tmp.GetXaxis().GetBinCenter(maxbin)
        tmpwidth = 0.5
        g1 = ROOT.TF1("g1","gaus", tmpmean-tmpwidth,tmpmean+tmpwidth)
        tmp.Fit(g1, "SR")
        c1 =ROOT.TCanvas("c","",800,800)
        tmp.Draw()
        c1.SaveAs("debug_%s_fit1_mvvres_%i.png"%(label,bin))
        tmpmean = g1.GetParameter(1)
        tmpwidth = g1.GetParameter(2)
        g1 = ROOT.TF1("g1","gaus", tmpmean-(tmpwidth*2),tmpmean+(tmpwidth*2))
        tmp.Fit(g1, "SR")
        c1 =ROOT.TCanvas("c","",800,800)
        tmp.Draw()
        c1.SaveAs("debug_%s_fit2_mvvres_%i.png"%(label,bin))
        tmpmean = g1.GetParameter(1)
        tmpmeanErr = g1.GetParError(1)
        tmpwidth = g1.GetParameter(2)
        tmpwidthErr = g1.GetParError(2)
        scalexHisto.SetBinContent(bin,tmpmean)
        scalexHisto.SetBinError  (bin,tmpmeanErr)
        resxHisto.SetBinContent  (bin,tmpwidth)
        resxHisto.SetBinError    (bin,tmpwidthErr)
     for bin in range(1,superHY.GetNbinsX()+1): 
        tmp=superHY.ProjectionY("q",bin,bin)
        if bin==1:
                 scaleyHisto.SetBinContent(bin,tmp.GetMean())
                 scaleyHisto.SetBinError(bin,tmp.GetMeanError())
                 resyHisto.SetBinContent(bin,tmp.GetRMS())
                 resyHisto.SetBinError(bin,tmp.GetRMSError())       
                 continue       
        startbin   = 0.
        maxcontent = 0.
        for b in range(tmp.GetXaxis().GetNbins()):
          if tmp.GetXaxis().GetBinCenter(b+1) > startbin and tmp.GetBinContent(b+1)>maxcontent:
            maxbin = b
            maxcontent = tmp.GetBinContent(b+1)
        tmpmean = tmp.GetXaxis().GetBinCenter(maxbin)
        tmpwidth = 0.3
        g1 = ROOT.TF1("g1","gaus", tmpmean-tmpwidth,tmpmean+tmpwidth)
        tmp.Fit(g1, "SR")
        c1 =ROOT.TCanvas("c","",800,800)
        tmp.Draw()
        c1.SaveAs("debug_%s_fit1_mjres_%i.png"%(label,bin))
        tmpmean = g1.GetParameter(1)
        tmpwidth = g1.GetParameter(2)
        g1 = ROOT.TF1("g1","gaus", tmpmean-(tmpwidth*1.1),tmpmean+(tmpwidth*1.1))
        tmp.Fit(g1, "SR")
        c1 =ROOT.TCanvas("c","",800,800)
        tmp.Draw()
        c1.SaveAs("debug_%s_fit2_mjres_%i.png"%(label,bin))
        tmpmean = g1.GetParameter(1)
        tmpmeanErr = g1.GetParError(1)
        tmpwidth = g1.GetParameter(2)
        tmpwidthErr = g1.GetParError(2)
        scaleyHisto.SetBinContent(bin,tmpmean)
        scaleyHisto.SetBinError  (bin,tmpmeanErr)
        resyHisto.SetBinContent  (bin,tmpwidth)
        resyHisto.SetBinError    (bin,tmpwidthErr)
         

         # tmp=superHX.ProjectionY("q",bin,bin)
#        scalexHisto.SetBinContent(bin,tmp.GetMean())
#        scalexHisto.SetBinError(bin,tmp.GetMeanError())
#        resxHisto.SetBinContent(bin,tmp.GetRMS())
#        resxHisto.SetBinError(bin,tmp.GetRMSError())
#
#        tmp=superHY.ProjectionY("q",bin,bin)
#        scaleyHisto.SetBinContent(bin,tmp.GetMean())
#        scaleyHisto.SetBinError(bin,tmp.GetMeanError())
#        resyHisto.SetBinContent(bin,tmp.GetRMS())
#        resyHisto.SetBinError(bin,tmp.GetRMSError())

         #tmp=superHNsubj.ProjectionY("q",bin,bin)
         #scaleNsubjHisto.SetBinContent(bin,tmp.GetMean())
         #scaleNsubjHisto.SetBinError(bin,tmp.GetMeanError())
         #resNsubjHisto.SetBinContent(bin,tmp.GetRMS())
         #resNsubjHisto.SetBinError(bin,tmp.GetRMSError())
         
     scalexHisto.Write()
     scaleyHisto.Write()
     #scaleNsubjHisto.Write()
     resxHisto.Write()
     resyHisto.Write()
     #resNsubjHisto.Write()
     superHX.Write("dataX")
     superHY.Write("dataY")
     #superHNsubj.Write("dataNsubj")

     fout.Close()
     fin.Close()
     
     os.system('rm '+f)
    
    #use the pythia det resolution for all the sample in the following steps
    if outname.find("nominal") != -1 and template==pythia_template :
        print " outname used for copy: "+str(outname)
        os.system( 'cp %s %s'%(outname,resFile) )
    elif outname.find("altshapeUp") != -1 and template==herwig_template:
        print " outname used for copy: "+str(outname)
        os.system( 'cp %s %s'%(outname,resFile) )
    elif outname.find("altshape2") != -1 and template==mg_template:
        print " outname used for copy: "+str(outname)
        os.system( 'cp %s %s'%(outname,resFile) )

def merge1DMVVTemplate(jobList,files,jobname,purity,binsMVV,minMVV,maxMVV,HCALbinsMVV,name,filename):
	print "Merging 1D templates"
	print
	print "Jobs to merge :   " ,jobList
	print "Files ran over:   " ,files
	
	outdir = 'res'+jobname
	jobdir = 'tmp'+jobname
	
	resubmit, jobsPerSample,exit_flag = getJobs(files,jobList,outdir,purity)
	
	if exit_flag:
	 submit = raw_input("The following files are missing: %s. Do you  want to resubmit the jobs to the batch system before merging? [y/n] "%resubmit)
	 if submit == 'y' or submit=='Y':
		 print "Resubmitting jobs:"
		 jobs = reSubmit(jobdir,resubmit,jobname)
		 waitForBatchJobs(jobname,len(resubmit),len(resubmit), userName, timeCheck)
		 resubmit, jobsPerSample,exit_flag = getJobs(files,jobList,outdir,purity)
		 if exit_flag: 
			 print "Job crashed again! Please resubmit manually before attempting to merge again"
			 for j in jobs: print j 
			 sys.exit()
	 else:
		 submit = raw_input("Some files are missing. [y] == Exit without merging, [n] == continue ? ")
		 if submit == 'y' or submit=='Y':
			 print "Exit without merging!"
			 sys.exit()
		 else:
			 print "Continuing merge!"
 
	try: 
		os.stat(outdir+'_out') 
		os.system('rm -r '+outdir+'_out')
		os.mkdir(outdir+'_out')
	except: os.mkdir(outdir+'_out')

	for s in jobsPerSample.keys():
	 factor = 1./float(len(jobsPerSample[s]))
	 print "sample: ", s,"number of files:",len(jobsPerSample[s]),"adding histo with scale factor:",factor
 
	 outf = ROOT.TFile.Open(outdir+'_out/%s_%s_MVV_%s_%s.root'%(filename,name,s,purity),'RECREATE')
  
	 finalHistos = {}
	 finalHistos['histo_nominal'] = ROOT.TH1F("histo_nominal_out","histo_nominal_out",binsMVV,minMVV,maxMVV)
	 finalHistos['mvv_nominal'] = ROOT.TH1F("mvv_nominal_out","mvv_nominal_out",binsMVV,minMVV,maxMVV)
	 if HCALbinsMVV!="":
             a,b,bins = HCALbinsMVV.split(" ")
             binning = getBinning(bins)
             binning = array("f",binning)
             finalHistos['histo_nominal'] = ROOT.TH1F("histo_nominal_out","histo_nominal_out",len(binning)-1,binning)
             finalHistos['mvv_nominal'] = ROOT.TH1F("mvv_nominal_out","mvv_nominal_out",len(binning)-1,binning)
    
	 for f in jobsPerSample[s]:
    
	  inf = ROOT.TFile.Open(f,'READ')
    
	  for h in inf.GetListOfKeys():
  
	   if (h.GetName() == 'histo_nominal' and h.GetTitle() == 'histo_nominal') or (h.GetName() == 'mvv_nominal' and h.GetTitle() == 'mvv_nominal'):

	    histo = ROOT.TH1F()
	    histo = inf.Get(h.GetName())

	    finalHistos[h.GetName()].Add(histo,factor)
      
   
	 print "Write file: ",outdir+'_out/%s_%s_MVV_%s_%s.root'%(filename,name,s,purity)
   
	 outf.cd()  
	 finalHistos['histo_nominal'].Write('histo_nominal')
	 finalHistos['mvv_nominal'].Write('mvv_nominal')
   
	 outf.Close()
	 outf.Delete()

	# read out files
	filelist = os.listdir('./'+outdir+'_out'+'/')

	mg_files = []
	pythia_files = []
	herwig_files = []
	dijet_files = []

	for f in filelist:
	 if f.find('QCD_HT') != -1: mg_files.append('./'+outdir+'_out'+'/'+f)
	 elif f.find('QCD_Pt_') != -1: pythia_files.append('./'+outdir+'_out'+'/'+f)
	 elif f.find('QCD_Pt-') != -1: herwig_files.append('./'+outdir+'_out'+'/'+f)
	 else: dijet_files.append('./'+outdir+'_out'+'/'+f)
	 
	doMadGraph = False
	doHerwig   = False
	doPythia   = False
	doDijet    = False
	
	#now hadd them
	if len(mg_files) > 0:
		cmd = 'hadd -f %s_%s_MVV_%s_altshape2.root '%(filename,name,purity)
		for f in mg_files:
		 cmd += f
		 cmd += ' '
		print cmd
		os.system(cmd)
		
		fhadd_madgraph = ROOT.TFile.Open('%s_%s_MVV_%s_altshape2.root'%(filename,name,purity),'READ')
		mvv_altshape2 = fhadd_madgraph.Get('mvv_nominal')
		mvv_altshape2.SetName('mvv_altshape2')
		mvv_altshape2.SetTitle('mvv_altshape2')
		histo_altshape2Up = fhadd_madgraph.Get('histo_nominal')
		histo_altshape2Up.SetName('histo_altshape2Up')
		histo_altshape2Up.SetTitle('histo_altshape2Up')
		
		doMadGraph = True
		
	if len(herwig_files) > 0:
		cmd = 'hadd -f %s_%s_MVV_%s_altshapeUp.root '%(filename,name,purity)
		for f in herwig_files:
		 cmd += f
		 cmd += ' '
		print cmd
		os.system(cmd)
		
		fhadd_herwig = ROOT.TFile.Open('%s_%s_MVV_%s_altshapeUp.root'%(filename,name,purity),'READ')
		mvv_altshapeUp = fhadd_herwig.Get('mvv_nominal')
		mvv_altshapeUp.SetName('mvv_altshapeUp')
		mvv_altshapeUp.SetTitle('mvv_altshapeUp')
		histo_altshapeUp = fhadd_herwig.Get('histo_nominal')
		histo_altshapeUp.SetName('histo_altshapeUp')
		histo_altshapeUp.SetTitle('histo_altshapeUp')
		
		doHerwig = True
 	
	if len(pythia_files) > 0:
		cmd = 'hadd -f %s_%s_MVV_%s_nominal.root '%(filename,name,purity)
		for f in pythia_files:
		 cmd += f
		 cmd += ' '
		print cmd
		os.system(cmd)
		
		fhadd_pythia = ROOT.TFile.Open('%s_%s_MVV_%s_nominal.root'%(filename,name,purity),'READ')
		mvv_nominal = fhadd_pythia.Get('mvv_nominal')
		mvv_nominal.SetName('mvv_nominal')
		mvv_nominal.SetTitle('mvv_nominal')
		histo_nominal = fhadd_pythia.Get('histo_nominal')
		histo_nominal.SetName('histo_nominal')
		histo_nominal.SetTitle('histo_nominal')
		
		doPythia = True

	if len(dijet_files) > 0:
		cmd = 'hadd -f %s_%s_MVV_%s_NLO.root '%(filename,name,purity)
		for f in dijet_files:
		 cmd += f
		 cmd += ' '
		print cmd
		os.system(cmd)
		
		fhadd_dijet = ROOT.TFile.Open('%s_%s_MVV_%s_NLO.root'%(purity),'READ')
		mvv_NLO = fhadd_dijet.Get('mvv_nominal')
		mvv_NLO.SetName('mvv_NLO')
		mvv_NLO.SetTitle('mvv_NLO')
		histo_NLO = fhadd_dijet.Get('histo_nominal')
		histo_NLO.SetName('histo_NLO')
		histo_NLO.SetTitle('histo_NLO')
		
		doDijet = True
		
	outf = ROOT.TFile.Open('%s_%s_MVV_%s.root'%(filename,name,purity),'RECREATE') 
    
	if doPythia:
		print "doing Pythia"
                mvv_nominal.Write('mvv_nominal')
		histo_nominal.Scale(1./histo_nominal.Integral())
		histo_nominal.Write('histo_nominal')
			
		print "Now pT"
		alpha=1.5/float(maxMVV)
		histogram_pt_up,histogram_pt_down=unequalScale(histo_nominal,"histo_nominal_PT",alpha)
		histogram_pt_down.SetName('histo_nominal_PTDown')
		histogram_pt_down.SetTitle('histo_nominal_PTDown')
		histogram_pt_down.Write('histo_nominal_PTDown')
		histogram_pt_up.SetName('histo_nominal_PTUp')
		histogram_pt_up.SetTitle('histo_nominal_PTUp')
		histogram_pt_up.Write('histo_nominal_PTUp')

                print "Now OPT"
		alpha=1.5*float(minMVV)
		histogram_opt_up,histogram_opt_down=unequalScale(histo_nominal,"histo_nominal_OPT",alpha,-1)
		histogram_opt_down.SetName('histo_nominal_OPTDown')
		histogram_opt_down.SetTitle('histo_nominal_OPTDown')
		histogram_opt_down.Write('histo_nominal_OPTDown')
		histogram_opt_up.SetName('histo_nominal_OPTUp')
		histogram_opt_up.SetTitle('histo_nominal_OPTUp')
		histogram_opt_up.Write('histo_nominal_OPTUp')
				
	if doHerwig:
                print "doing Herwig"
		mvv_altshapeUp.Write('mvv_altshapeUp')
		histo_altshapeUp.Write('histo_altshapeUp')
    
		print "Now pT"
		alpha=1.5/float(maxMVV)
		histogram_altshapeUp_pt_up,histogram_altshapeUp_pt_down=unequalScale(histo_nominal,"histo_altshapeUp_PT",alpha)
		histogram_altshapeUp_pt_down.SetName('histo_altshapeUp_PTDown')
		histogram_altshapeUp_pt_down.SetTitle('histo_altshapeUp_PTDown')
		histogram_altshapeUp_pt_down.Write('histo_altshapeUp_PTDown')
		histogram_altshapeUp_pt_up.SetName('histo_altshapeUp_PTUp')
		histogram_altshapeUp_pt_up.SetTitle('histo_altshapeUp_PTUp')
		histogram_altshapeUp_pt_up.Write('histo_altshapeUp_PTUp')

                print "Now OPT"
		alpha=1.5*float(minMVV)
		histogram_altshapeUp_opt_up,histogram_altshapeUp_opt_down=unequalScale(histo_altshapeUp,"histo_altshapeUp_OPT",alpha,-1)
		histogram_altshapeUp_opt_down.SetName('histo_altshapeUp_OPTDown')
		histogram_altshapeUp_opt_down.SetTitle('histo_altshapeUp_OPTDown')
		histogram_altshapeUp_opt_down.Write('histo_altshapeUp_OPTDown')
		histogram_altshapeUp_opt_up.SetName('histo_altshapeUp_OPTUp')
		histogram_altshapeUp_opt_up.SetTitle('histo_altshapeUp_OPTUp')
		histogram_altshapeUp_opt_up.Write('histo_altshapeUp_OPTUp')

		if doPythia:
			histogram_altshapeDown=mirror(histo_altshapeUp,histo_nominal,"histo_altshapeDown")
			histogram_altshapeDown.SetName('histo_altshapeDown')
			histogram_altshapeDown.SetTitle('histo_altshapeDown')
			histogram_altshapeDown.Write('histo_altshapeDown')
		
	if doMadGraph:
                print "doing Madgraph"
		mvv_altshape2.Write('mvv_altshape2')
		histo_altshape2Up.Write('histo_altshape2Up')
    
		print "Now pT"
		alpha=1.5/float(maxMVV)
                print " #################   careful!!! newt line are commented out because histo_nominal is missing at the moment!!! ####################"
	        histogram_altshape2_pt_up,histogram_altshape2_pt_down=unequalScale(histo_nominal,"histo_altshape2_PT",alpha)
	        histogram_altshape2_pt_down.SetName('histo_altshape2_PTDown')
	        histogram_altshape2_pt_down.SetTitle('histo_altshape2_PTDown')
	        histogram_altshape2_pt_down.Write('histo_altshape2_PTDown')
	        histogram_altshape2_pt_up.SetName('histo_altshape2_PTUp')
	        histogram_altshape2_pt_up.SetTitle('histo_altshape2_PTUp')
	        histogram_altshape2_pt_up.Write('histo_altshape2_PTUp')

                print "Now OPT"
		alpha=1.5*float(minMVV)
		histogram_altshape2_opt_up,histogram_altshape2_opt_down=unequalScale(histo_altshape2Up,"histo_altshape2_OPT",alpha,-1)
		histogram_altshape2_opt_down.SetName('histo_altshape2_OPTDown')
		histogram_altshape2_opt_down.SetTitle('histo_altshape2_OPTDown')
		histogram_altshape2_opt_down.Write('histo_altshape2_OPTDown')
		histogram_altshape2_opt_up.SetName('histo_altshape2_OPTUp')
		histogram_altshape2_opt_up.SetTitle('histo_altshape2_OPTUp')
		histogram_altshape2_opt_up.Write('histo_altshape2_OPTUp')
    
		if doPythia:
			histogram_altshape2Down=mirror(histo_altshape2Up,histo_nominal,"histo_altshape2Down")
			histogram_altshape2Down.SetName('histo_altshape2Down')
			histogram_altshape2Down.SetTitle('histo_altshape2Down')
			histogram_altshape2Down.Write('histo_altshape2Down')
			
	if doDijet:
		mvv_NLO.Write('mvv_NLO')
		histo_NLO.Write('histo_NLO')
		if doPythia:
			histogram_NLODown=mirror(histo_NLO,histo_nominal,"histo_NLODown")
			histogram_NLODown.SetName('histo_NLODown')
			histogram_NLODown.SetTitle('histo_NLODown')
			histogram_NLODown.Write('histo_NLODown')
						
	os.system('rm -rf '+outdir+'_out/')

def merge2DTemplate(jobList,files,jobname,purity,leg,binsMVV,binsMJ,minMVV,maxMVV,minMJ,maxMJ,HCALbinsMVV,name,filename):  
  
  print "Merging 2D templates"
  print
  print "Jobs to merge :   " ,jobList
  print "Files ran over:   " ,files
  
  outdir = 'res'+jobname
  jobdir = 'tmp'+jobname
  
  resubmit, jobsPerSample,exit_flag = getJobs(files,jobList,outdir,purity)
  
  if exit_flag:
   submit = raw_input("The following files are missing: %s. Do you  want to resubmit the jobs to the batch system before merging? [y/n] "%resubmit)
   if submit == 'y' or submit=='Y':
     print "Resubmitting jobs:"
     jobs = reSubmit(jobdir,resubmit,jobname)
     waitForBatchJobs(jobname,len(resubmit),len(resubmit), userName, timeCheck)
     resubmit, jobsPerSample,exit_flag = getJobs(files,jobList,outdir,purity)
     if exit_flag: 
       print "Job crashed again! Please resubmit manually before attempting to merge again"
       for j in jobs: print j 
       sys.exit()
   else:
     submit = raw_input("Some files are missing. [y] == Exit without merging, [n] == continue ? ")
     if submit == 'y' or submit=='Y':
       print "Exit without merging!"
       sys.exit()
     else:
       print "Continuing merge!"
  
 
  try: 
    os.stat(outdir+'_out') 
    os.system('rm -r '+outdir+'_out')
    os.mkdir(outdir+'_out')
  except: os.mkdir(outdir+'_out')

  for s in jobsPerSample.keys():

   factor = 1./float(len(jobsPerSample[s]))
   print "sample: ", s,"number of files:",len(jobsPerSample[s]),"adding histo with scale factor:",factor
 
   outf = ROOT.TFile.Open(outdir+'_out/%s_%s_COND2D_%s_%s_%s.root'%(filename,name,s,leg,purity),'RECREATE')
  
   finalHistos = {}
   finalHistos['histo_nominal_coarse'] = ROOT.TH2F("histo_nominal_coarse_out","histo_nominal_coarse_out",binsMJ,minMJ,maxMJ,binsMVV,minMVV,maxMVV)
   finalHistos['mjet_mvv_nominal'] = ROOT.TH2F("mjet_mvv_nominal_out","mjet_mvv_nominal_out",binsMJ,minMJ,maxMJ,binsMVV,minMVV,maxMVV)
   finalHistos['mjet_mvv_nominal_3D'] = ROOT.TH3F("mjet_mvv_nominal_3D_out","mjet_mvv_nominal_3D_out",binsMJ,minMJ,maxMJ,binsMJ,minMJ,maxMJ,binsMVV,minMVV,maxMVV)
   if HCALbinsMVV!="":
       a,b,bins = HCALbinsMVV.split(" ")
       binning = getBinning(bins)
       binning = array("d",binning)
       xbins = []
       for i in range(0,binsMJ+1):
          xbins.append(minMJ + i* (maxMJ - minMJ)/binsMJ)
       xbins = array("d",xbins)
       finalHistos['histo_nominal_coarse'] = ROOT.TH2F("histo_nominal_coarse_out","histo_nominal_coarse_out",len(xbins)-1,xbins,len(binning)-1,binning)
       finalHistos['mjet_mvv_nominal'] = ROOT.TH2F("mjet_mvv_nominal_out","mjet_mvv_nominal_out",len(xbins)-1,xbins,len(binning)-1,binning)
       finalHistos['mjet_mvv_nominal_3D'] = ROOT.TH3F("mjet_mvv_nominal_3D_out","mjet_mvv_nominal_3D_out",len(xbins)-1,xbins,len(xbins)-1,xbins,len(binning)-1,binning)
       print "use binning "+str(binning)
   print binsMVV
   print finalHistos['histo_nominal_coarse'].GetNbinsY()
   for f in jobsPerSample[s]:

    inf = ROOT.TFile.Open(f,'READ')
    #print "open file "+f
    
    for h in inf.GetListOfKeys():
  
     for k in finalHistos.keys():
        
      if h.GetName() == k:
       histo = ROOT.TH1F()
       histo = inf.Get(h.GetName())

       finalHistos[h.GetName()].Add(histo,factor)
   
   print "Write file: ",outdir+'_out/%s_%s_COND2D_%s_%s_%s.root'%(filename,name,s,leg,purity)
   
   outf.cd()  
 
   for k in finalHistos.keys():
    finalHistos[k].SetTitle(k)
    finalHistos[k].Write(k)
   
   outf.Close()
   outf.Delete()


  # read out files
  filelist = os.listdir('./'+outdir+'_out'+'/')

  mg_files = []
  pythia_files = []
  herwig_files = []
  dijet_files = []

  for f in filelist:
   if f.find('COND2D') == -1: continue
   if f.find('QCD_HT') != -1: mg_files.append('./'+outdir+'_out'+'/'+f)
   elif f.find('QCD_Pt_') != -1: pythia_files.append('./'+outdir+'_out'+'/'+f)
   elif f.find('QCD_Pt-') != -1: herwig_files.append('./'+outdir+'_out'+'/'+f)
   else: dijet_files.append('./'+outdir+'_out'+'/'+f)
   
   
  doMadGraph = False
  doHerwig   = False
  doPythia   = False
  doDijet    = False
  
  #now hadd them
  if len(mg_files) > 0:
    print "doing MadGraph "
    cmd = 'hadd -f %s_%s_COND2D_%s_%s_altshape2.root '%(filename,name,purity,leg)
    for f in mg_files:
     cmd += f
     cmd += ' '
    print cmd
    os.system(cmd)


    fhadd_madgraph = ROOT.TFile.Open('%s_%s_COND2D_%s_%s_altshape2.root'%(filename,name,purity,leg),'READ')
    mjet_mvv_altshape2_3D = fhadd_madgraph.Get('mjet_mvv_nominal_3D') 
    mjet_mvv_altshape2_3D.SetName('mjet_mvv_altshape2_3D')
    mjet_mvv_altshape2_3D.SetTitle('mjet_mvv_altshape2_3D')
    mjet_mvv_altshape2 = fhadd_madgraph.Get('mjet_mvv_nominal')
    mjet_mvv_altshape2.SetName('mjet_mvv_altshape2')
    mjet_mvv_altshape2.SetTitle('mjet_mvv_altshape2')
    histo_altshape2Up = fhadd_madgraph.Get('histo_nominal_coarse')
    histo_altshape2Up.SetName('histo_altshape2_coarse')
    histo_altshape2Up.SetTitle('histo_altshape2_coarse')
    #histo_altshape2 = fhadd_madgraph.Get('histo_nominal')
    #histo_altshape2.SetName('histo_altshape2')
    #histo_altshape2.SetTitle('histo_altshape2')

    doMadGraph = True


  if len(herwig_files) > 0:
    print "doing Herwig"
    cmd = 'hadd -f %s_%s_COND2D_%s_%s_altshapeUp.root '%(filename,name,purity,leg)
    for f in herwig_files:
     cmd += f
     cmd += ' '
    print cmd
    os.system(cmd)

    fhadd_herwig = ROOT.TFile.Open('%s_%s_COND2D_%s_%s_altshapeUp.root'%(filename,name,purity,leg),'READ')
    mjet_mvv_altshapeUp_3D = fhadd_herwig.Get('mjet_mvv_nominal_3D') 
    mjet_mvv_altshapeUp_3D.SetName('mjet_mvv_altshapeUp_3D')
    mjet_mvv_altshapeUp_3D.SetTitle('mjet_mvv_altshapeUp_3D')
    mjet_mvv_altshapeUp = fhadd_herwig.Get('mjet_mvv_nominal')
    mjet_mvv_altshapeUp.SetName('mjet_mvv_altshapeUp')
    mjet_mvv_altshapeUp.SetTitle('mjet_mvv_altshapeUp')
    histo_altshapeUp = fhadd_herwig.Get('histo_nominal_coarse')
    histo_altshapeUp.SetName('histo_altshapeUp_coarse')
    histo_altshapeUp.SetTitle('histo_altshapeUp_coarse')
    #histo_altshapeUp = fhadd_herwig.Get('histo_nominal')
    #histo_altshapeUp.SetName('histo_altshapeUp')
    #histo_altshapeUp.SetTitle('histo_altshapeUp')
    doHerwig = True

  if len(pythia_files) > 0:
    print "doing pythia"
    cmd = 'hadd -f %s_%s_COND2D_%s_%s_nominal.root '%(filename,name,purity,leg)
    for f in pythia_files:
     cmd += f
     cmd += ' '
    print cmd
    os.system(cmd)
    
    fhadd_pythia = ROOT.TFile.Open('%s_%s_COND2D_%s_%s_nominal.root'%(filename,name,purity,leg),'READ')
    mjet_mvv_nominal_3D = fhadd_pythia.Get('mjet_mvv_nominal_3D') 
    mjet_mvv_nominal_3D.SetName('mjet_mvv_nominal_3D')
    mjet_mvv_nominal_3D.SetTitle('mjet_mvv_nominal_3D')
    mjet_mvv_nominal = fhadd_pythia.Get('mjet_mvv_nominal')
    mjet_mvv_nominal.SetName('mjet_mvv_nominal')
    mjet_mvv_nominal.SetTitle('mjet_mvv_nominal')
    fhadd_pythia.Print()
    histo_nominal = fhadd_pythia.Get('histo_nominal_coarse')
    histo_nominal.SetName('histo_nominal_coarse')
    histo_nominal.SetTitle('histo_nominal_coarse')
    #histo_nominal = fhadd_pythia.Get('histo_nominal')
    #histo_nominal.SetName('histo_nominal')
    #histo_nominal.SetTitle('histo_nominal')
    
    doPythia = True

  if len(dijet_files) > 0:
    cmd = 'hadd -f %s_%s_COND2D_%s_%s_NLO.root '%(filename,name,purity,leg)
    for f in dijet_files:
     cmd += f
     cmd += ' '
    print cmd
    os.system(cmd)

    fhadd_dijet = ROOT.TFile.Open('%s_%s_COND2D_%s_%s_NLO.root'%(filename,name,purity,leg),'READ')
    mjet_mvv_NLO_3D = fhadd_dijet.Get('mjet_mvv_nominal_3D') 
    mjet_mvv_NLO_3D.SetName('mjet_mvv_NLO_3D')
    mjet_mvv_NLO_3D.SetTitle('mjet_mvv_NLO_3D')
    mjet_mvv_NLO = fhadd_dijet.Get('mjet_mvv_nominal')
    mjet_mvv_NLO.SetName('mjet_mvv_NLO')
    mjet_mvv_NLO.SetTitle('mjet_mvv_NLO')
    histo_NLO = fhadd_dijet.Get('histo_nominal_coarse')
    histo_NLO.SetName('histo_NLO_coarse')
    histo_NLO.SetTitle('histo_NLO_coarse')
    doDijet = True
    
    
  outf = ROOT.TFile.Open('%s_%s_COND2D_%s_%s.root'%(filename,name,purity,leg),'RECREATE') 
  finalHistograms = {}
  
  if doPythia:
    mjet_mvv_nominal.Write('mjet_mvv_nominal')
    mjet_mvv_nominal_3D.Write('mjet_mvv_nominal_3D')

    histo_nominal.Write('histo_nominal_coarse')
    print "make conditional histogram"
    conditional(histo_nominal)
    
    if HCALbinsMVV!="":
      expanded=expandHistoBinned(histo_nominal,"",xbins,binning)
    else:
      expanded=expandHisto(histo_nominal,"",binsMVV,binsMJ,minMVV,maxMVV,minMJ,maxMJ)
    conditional(expanded)
    expanded.SetName('histo_nominal')
    expanded.SetTitle('histo_nominal')
    expanded.Write('histo_nominal')
    finalHistograms['histo_nominal'] = histo_nominal
    
    alpha=1.5/float(maxMJ)
    histogram_pt_up,histogram_pt_down=unequalScale(finalHistograms['histo_nominal'],"histo_nominal_PT",alpha,1,2)
    conditional(histogram_pt_down)
    histogram_pt_down.SetName('histo_nominal_PTDown')
    histogram_pt_down.SetTitle('histo_nominal_PTDown')
    histogram_pt_down.Write('histo_nominal_PTDown')
    conditional(histogram_pt_up)
    histogram_pt_up.SetName('histo_nominal_PTUp')
    histogram_pt_up.SetTitle('histo_nominal_PTUp')
    histogram_pt_up.Write('histo_nominal_PTUp')

    alpha=1.5*float(minMJ)
    h1,h2=unequalScale(finalHistograms['histo_nominal'],"histo_nominal_OPT",alpha,-1,2)
    conditional(h1)
    h1.SetName('histo_nominal_OPTUp')
    h1.SetTitle('histo_nominal_OPTUp')
    h1.Write('histo_nominal_OPTUp')
    conditional(h2)
    h2.SetName('histo_nominal_OPTDown')
    h2.SetTitle('histo_nominal_OPTDown')
    h2.Write('histo_nominal_OPTDown')

    alpha=float(maxMJ)*float(maxMJ)
    histogram_pt2_up,histogram_pt2_down=unequalScale(finalHistograms['histo_nominal'],"histo_nominal_PT2",alpha,2,2)
    conditional(histogram_pt2_down)
    histogram_pt2_down.SetName('histo_nominal_PT2Down')
    histogram_pt2_down.SetTitle('histo_nominal_PT2Down')
    histogram_pt2_down.Write('histo_nominal_PT2Down')
    conditional(histogram_pt2_up)
    histogram_pt2_up.SetName('histo_nominal_PT2Up')
    histogram_pt2_up.SetTitle('histo_nominal_PT2Up')
    histogram_pt2_up.Write('histo_nominal_PT2Up')

    alpha=float(minMJ)*float(minMJ)
    histogram_opt2_up,histogram_opt2_down=unequalScale(finalHistograms['histo_nominal'],"histo_nominal_OPT2",alpha,-2,2)
    conditional(histogram_opt2_down)
    histogram_opt2_down.SetName('histo_nominal_OPT2Down')
    histogram_opt2_down.SetTitle('histo_nominal_OPT2Down')
    histogram_opt2_down.Write('histo_nominal_OPT2Down')
    conditional(histogram_opt2_up)
    histogram_opt2_up.SetName('histo_nominal_OPT2Up')
    histogram_opt2_up.SetTitle('histo_nominal_OPT2Up')
    histogram_opt2_up.Write('histo_nominal_OPT2Up')
        
  if doHerwig:
    histo_altshapeUp.Write('histo_altshapeUp_coarse')
    conditional(histo_altshapeUp)
    expanded=expandHisto(histo_altshapeUp,"herwig",binsMVV,binsMJ,minMVV,maxMVV,minMJ,maxMJ)
    if HCALbinsMVV!="":
      expanded=expandHistoBinned(histo_altshapeUp,"",xbins,binning)
    conditional(expanded)
    expanded.SetName('histo_altshapeUp')
    expanded.SetTitle('histo_altshapeUp')
    expanded.Write('histo_altshapeUp')
    finalHistograms['histo_altshapeUp'] = expanded
    
    
    alpha=1.5/float(maxMJ)
    histogram_pt_up,histogram_pt_down=unequalScale(finalHistograms['histo_altshapeUp'],"histo_altshapeUp_PT",alpha,1,2)
    conditional(histogram_pt_down)
    histogram_pt_down.SetName('histo_altshapeUp_PTDown')
    histogram_pt_down.SetTitle('histo_altshapeUp_PTDown')
    histogram_pt_down.Write('histo_altshapeUp_PTDown')
    conditional(histogram_pt_up)
    histogram_pt_up.SetName('histo_altshapeUp_PTUp')
    histogram_pt_up.SetTitle('histo_altshapeUp_PTUp')
    histogram_pt_up.Write('histo_altshapeUp_PTUp')

    alpha=1.5*float(minMJ)
    h1,h2=unequalScale(finalHistograms['histo_altshapeUp'],"histo_altshapeUp_OPT",alpha,-1,2)
    conditional(h1)
    h1.SetName('histo_altshapeUp_OPTUp')
    h1.SetTitle('histo_altshapeUp_OPTUp')
    h1.Write('histo_altshapeUp_OPTUp')
    conditional(h2)
    h2.SetName('histo_altshapeUp_OPTDown')
    h2.SetTitle('histo_altshapeUp_OPTDown')
    h2.Write('histo_altshapeUp_OPTDown')

    alpha=float(maxMJ)*float(maxMJ)
    histogram_pt2_up,histogram_pt2_down=unequalScale(finalHistograms['histo_altshapeUp'],"histo_altshapeUp_PT2",alpha,2,2)
    conditional(histogram_pt2_down)
    histogram_pt2_down.SetName('histo_altshapeUp_PT2Down')
    histogram_pt2_down.SetTitle('histo_altshapeUp_PT2Down')
    histogram_pt2_down.Write('histo_altshapeUp_PT2Down')
    conditional(histogram_pt2_up)
    histogram_pt2_up.SetName('histo_altshapeUp_PT2Up')
    histogram_pt2_up.SetTitle('histo_altshapeUp_PT2Up')
    histogram_pt2_up.Write('histo_altshapeUp_PT2Up')

    alpha=float(minMJ)*float(minMJ)
    histogram_opt2_up,histogram_opt2_down=unequalScale(finalHistograms['histo_altshapeUp'],"histo_altshapeUp_OPT2",alpha,-2,2)
    conditional(histogram_opt2_down)
    histogram_opt2_down.SetName('histo_altshapeUp_OPT2Down')
    histogram_opt2_down.SetTitle('histo_altshapeUp_OPT2Down')
    histogram_opt2_down.Write('histo_altshapeUp_OPT2Down')
    conditional(histogram_opt2_up)
    histogram_opt2_up.SetName('histo_altshapeUp_OPT2Up')
    histogram_opt2_up.SetTitle('histo_altshapeUp_OPT2Up')
    histogram_opt2_up.Write('histo_altshapeUp_OPT2Up')
    
    if doPythia:
      histogram_altshapeDown=mirror(finalHistograms['histo_altshapeUp'],finalHistograms['histo_nominal'],"histo_altshapeDown",2)
      conditional(histogram_altshapeDown)
      histogram_altshapeDown.SetName('histo_altshapeDown')
      histogram_altshapeDown.SetTitle('histo_altshapeDown')
      histogram_altshapeDown.Write()

  if doMadGraph:
    histo_altshape2Up.Write('histo_altshape2_coarse')
    conditional(histo_altshape2Up)
    expanded=expandHisto(histo_altshape2Up,"madgraph",binsMVV,binsMJ,minMVV,maxMVV,minMJ,maxMJ)
    if HCALbinsMVV!="":
      expanded=expandHistoBinned(histo_altshape2Up,"",xbins,binning)
    conditional(expanded)
    expanded.SetName('histo_altshape2Up')
    expanded.SetTitle('histo_altshape2Up')
    expanded.Write('histo_altshape2Up')
    finalHistograms['histo_altshape2Up'] = expanded
    
    alpha=1.5/float(maxMJ)
    histogram_pt_up,histogram_pt_down=unequalScale(finalHistograms['histo_altshape2Up'],"histo_altshape2_PT",alpha,1,2)
    conditional(histogram_pt_down)
    histogram_pt_down.SetName('histo_altshape2_PTDown')
    histogram_pt_down.SetTitle('histo_altshape2_PTDown')
    histogram_pt_down.Write('histo_altshape2_PTDown')
    conditional(histogram_pt_up)
    histogram_pt_up.SetName('histo_altshape2_PTUp')
    histogram_pt_up.SetTitle('histo_altshape2_PTUp')
    histogram_pt_up.Write('histo_altshape2_PTUp')

    alpha=1.5*float(minMJ)
    h1,h2=unequalScale(finalHistograms['histo_altshape2Up'],"histo_altshape2_OPT",alpha,-1,2)
    conditional(h1)
    h1.SetName('histo_altshape2_OPTUp')
    h1.SetTitle('histo_altshape2_OPTUp')
    h1.Write('histo_altshape2_OPTUp')
    conditional(h2)
    h2.SetName('histo_altshape2_OPTDown')
    h2.SetTitle('histo_altshape2_OPTDown')
    h2.Write('histo_altshape2_OPTDown')

    alpha=float(maxMJ)*float(maxMJ)
    histogram_pt2_up,histogram_pt2_down=unequalScale(finalHistograms['histo_altshape2Up'],"histo_altshape2_PT2",alpha,2,2)
    conditional(histogram_pt2_down)
    histogram_pt2_down.SetName('histo_altshape2_PT2Down')
    histogram_pt2_down.SetTitle('histo_altshape2_PT2Down')
    histogram_pt2_down.Write('histo_altshape2_PT2Down')
    conditional(histogram_pt2_up)
    histogram_pt2_up.SetName('histo_altshape2_PT2Up')
    histogram_pt2_up.SetTitle('histo_altshape2_PT2Up')
    histogram_pt2_up.Write('histo_altshape2_PT2Up')

    alpha=float(minMJ)*float(minMJ)
    histogram_opt2_up,histogram_opt2_down=unequalScale(finalHistograms['histo_altshape2Up'],"histo_altshape2_OPT2",alpha,-2,2)
    conditional(histogram_opt2_down)
    histogram_opt2_down.SetName('histo_altshape2_OPT2Down')
    histogram_opt2_down.SetTitle('histo_altshape2_OPT2Down')
    histogram_opt2_down.Write('histo_altshape2_OPT2Down')
    conditional(histogram_opt2_up)
    histogram_opt2_up.SetName('histo_altshape2_OPT2Up')
    histogram_opt2_up.SetTitle('histo_altshape2_OPT2Up')
    histogram_opt2_up.Write('histo_altshape2_OPT2Up')
    
    if doPythia:
      histogram_altshape2Down=mirror(finalHistograms['histo_altshape2Up'],finalHistograms['histo_nominal'],"histo_altshape2Down",2)
      conditional(histogram_altshape2Down)
      histogram_altshape2Down.SetName('histo_altshape2Down')
      histogram_altshape2Down.SetTitle('histo_altshape2Down')
      histogram_altshape2Down.Write()
      
  if doDijet:
    histo_NLO.Write('histo_NLO_coarse')
    conditional(histo_NLO)
    expanded=expandHisto(histo_NLO,"NLO",binsMVV,binsMJ,minMVV,maxMVV,minMJ,maxMJ)
    if HCALbinsMVV!="":
                    expanded=expandHistoBinned(histo_NLO,"",xbins,binning)
    conditional(expanded)
    expanded.SetName('histo_NLO')
    expanded.SetTitle('histo_NLO')
    expanded.Write('histo_NLO')
    finalHistograms['histo_NLO'] = expanded
    if doPythia:
      histogram_NLODown=mirror(finalHistograms['histo_NLO'],finalHistograms['histo_nominal'],"histo_NLODown",2)
      conditional(histogram_NLODown)
      histogram_NLODown.SetName('histo_NLODown')
      histogram_NLODown.SetTitle('histo_NLODown')
      histogram_NLODown.Write()
          
  os.system('rm -r '+outdir+'_out')
  # os.system('rm -r '+outdir)
  
def makeData(template,cut,rootFile,binsMVV,binsMJ,minMVV,maxMVV,minMJ,maxMJ,factors,name,data,jobname,samples,wait=True,binning='',addOption=""):
    print 
    print 'START: makeData'
    print "template = ",template
    print "cut      = ",cut
    print "rootFile = ",rootFile
    print "binsMVV  = ",binsMVV
    print "binsMJ   = ",binsMJ
    print "minMVV   = ",minMVV
    print "maxMVV   = ",maxMVV
    print "minMJ    = ",minMJ
    print "maxMJ    = ",maxMJ
    print "factor   = ",factors
    print "name     = ",name
    print "data     = ",data
    print "jobname  = ",jobname
    print "samples  = ",samples
    print 
    files = []
    sampleTypes = template.split(',')
    for f in os.listdir(samples):
        for t in sampleTypes:
            if f.find('.root') != -1 and f.find(t) != -1: files.append(f)
 
    NumberOfJobs= len(files) 
    OutputFileNames = rootFile.replace(".root","")
    cmd='vvMakeData.py -d {data} -c "{cut}"  -v "jj_l1_softDrop_mass,jj_l2_softDrop_mass,jj_LV_mass" {binning} -b "{bins},{bins},{BINS}" -m "{mini},{mini},{MINI}" -M "{maxi},{maxi},{MAXI}" -f {factors} -n "{name}" {addOption} {infolder} '.format(cut=cut,BINS=binsMVV,bins=binsMJ,MINI=minMVV,MAXI=maxMVV,mini=minMJ,maxi=maxMJ,factors=factors,name=name,data=data,infolder=samples,binning=binning,addOption=addOption)  
    queue = "1nd" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 
    path = os.getcwd()
    try: os.system("rm -r tmp"+jobname)
    except: print "No tmp/ directory"
    os.system("mkdir tmp"+jobname)
    try: os.stat("res"+jobname) 
    except: os.mkdir("res"+jobname)
        
    
    ##### Creating and sending jobs #####
    joblist = []
    ###### loop for creating and sending jobs #####
    for x in range(1, int(NumberOfJobs)+1):
       os.system("mkdir tmp"+jobname+"/"+str(files[x-1]).replace(".root",""))
       os.chdir("tmp"+jobname+"/"+str(files[x-1]).replace(".root",""))
       #os.system("mkdir tmp"+jobname+"/"+str(files[x-1]).replace(".root",""))
       #os.chdir(path+"tmp"+jobname+"/"+str(files[x-1]).replace(".root",""))
     
       with open('job_%s.sh'%files[x-1].replace(".root",""), 'w') as fout:
          fout.write("#!/bin/sh\n")
          fout.write("echo\n")
          fout.write("echo\n")
          fout.write("echo 'START---------------'\n")
          fout.write("echo 'WORKDIR ' ${PWD}\n")
          fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
          fout.write("cd "+str(path)+"\n")
          fout.write("cmsenv\n")
          fout.write(cmd+" -o "+path+"/res"+jobname+"/"+OutputFileNames+"_"+files[x-1]+" -s "+files[x-1]+"\n")
          fout.write("echo 'STOP---------------'\n")
          fout.write("echo\n")
          fout.write("echo\n")
       os.system("chmod 755 job_%s.sh"%(files[x-1].replace(".root","")) )

       if useCondorBatch:
           os.system("mv  job_*.sh "+jobname+".sh")
           makeSubmitFileCondor(jobname+".sh",jobname,"workday")
           os.system("condor_submit submit.sub")
       else:
           os.system("bsub -q "+queue+" -o logs job_%s.sh -J %s"%(files[x-1].replace(".root",""),jobname))
       print "job nr " + str(x) + " submitted"
       joblist.append("%s"%(files[x-1].replace(".root","")))
       os.chdir("../..")
   
    print
    print "your jobs:"
    if useCondorBatch:
        os.system("condor_q")
    else:
        os.system("bjobs")
    userName=os.environ['USER']
    if wait: waitForBatchJobs(jobname,NumberOfJobs,NumberOfJobs, userName, timeCheck)
    
    print
    print 'END: makeData'
    print
    return joblist, files

def mergeData(jobname,purity,rootFile,filename,name):
    
    print "Merging data from job " ,jobname
    print "Purity is " ,purity
    # read out files
    filelist = os.listdir('./res'+jobname+'/')

    mg_files     = []
    pythia_files = []
    herwig_files = []
    data_files   = []

    for f in filelist:
     #if f.find('COND2D') == -1: continue
     if f.find(purity)==-1:
             continue
     if f.find('QCD_HT')    != -1: mg_files.append('./res'+jobname+'/'+f)
     elif f.find('QCD_Pt_') != -1: pythia_files.append('./res'+jobname+'/'+f)
     elif f.find('JetHT')   != -1: data_files.append('./res'+jobname+'/'+f)
     else: herwig_files.append('./res'+jobname+'/'+f)

    #now hadd them
    if len(mg_files) > 0:
        cmd = 'hadd -f %s_%s_%s_altshape2.root '%(filename,name,purity)
        for f in mg_files:
         cmd += f
         cmd += ' '
        print cmd
        os.system(cmd)

    if len(herwig_files) > 0:
        cmd = 'hadd -f %s_%s_%s_altshapeUp.root '%(filename,name,purity)
        for f in herwig_files:
         cmd += f
         cmd += ' '
        print cmd
        os.system(cmd)

    if len(pythia_files) > 0:
        cmd = 'hadd -f %s '%rootFile
        for f in pythia_files:
         cmd += f
         cmd += ' '
        print cmd
        os.system(cmd)
    
    if len(data_files) > 0:
        cmd = 'hadd -f JJ_%s.root '%purity
        for f in data_files:
         cmd += f
         cmd += ' '
        print cmd
        os.system(cmd)
    print "Done merging data!"

def getListOfBinsLowEdge(hist,dim):
    axis =0
    N = 0
    if dim =="x":
        axis= hist.GetXaxis()
        N = hist.GetNbinsX()
    if dim =="y":
        axis = hist.GetYaxis()
        N = hist.GetNbinsY()
    if dim =="z":
        axis = hist.GetZaxis()
        N = hist.GetNbinsZ()
    if axis==0:
        return {}
    
    mmin = axis.GetXmin()
    mmax = axis.GetXmax()
    r =[]
    for i in range(1,N+2):
        #v = mmin + i * (mmax-mmin)/float(N)
        r.append(axis.GetBinLowEdge(i))
    return r

def makePseudoData(input="JJ_nonRes_LPLP.root",kernel="JJ_nonRes_3D_LPLP.root",mc="pythia",output="JJ_LPLP.root",lumi=35900):

 pwd = os.getcwd()
 
 ROOT.gRandom.SetSeed(0)
 
 finmc = ROOT.TFile.Open(pwd+'/'+input,'READ')
 hmcin = finmc.Get('nonRes')

 findata = ROOT.TFile.Open(pwd+'/'+kernel,'READ')
 #findata.ls()
 hdata = ROOT.TH3F()
 
 if   mc == 'pythia': hdata = findata.Get('histo')
 elif mc == 'herwig': hdata = findata.Get('histo_altshapeUp')
 elif mc == 'madgraph': hdata = findata.Get('histo_altshape2')
 elif mc == 'powheg': hdata = findata.Get('histo_NLO')
 
 fout = ROOT.TFile.Open(output,'RECREATE')
 #hmcin.Scale(10.)
 hmcin.Write('nonRes')

 xbins = array("f",getListOfBinsLowEdge(hmcin,"x"))
 zbins = array("f",getListOfBinsLowEdge(hmcin,"z"))
 hout = ROOT.TH3F('data','data',len(xbins)-1,xbins,len(xbins)-1,xbins,len(zbins)-1,zbins)
 hout.FillRandom(hdata,int(hmcin.Integral()*lumi))
 hout.Write('data')
 print "Writing histograms nonRes and data to file ", output

 finmc.Close()
 findata.Close()
 fout.Close()
 

def submitCPs(samples,template,wait,jobname="CPs",rootFile="controlplots_2017.root"):
  print 
  print 'START: submitCPs'
  print "template = ",template
  print "jobname  = ",jobname
  print "samples  = ",samples
  print 
  files = []
  sampleTypes = template.split(',')
  for f in os.listdir(samples):
    for t in sampleTypes:
      if f.find(t) == -1: continue 
      if f.startswith('.'): continue
      if f.find('.root') != -1 and f.find('rawPUMC') == -1: 
        print f
        files.append(f)
  
  NumberOfJobs= len(files)
  OutputFileNames = rootFile.replace(".root","")
  cmd = "python submit_CP.py"
  queue = "8nh" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw

  try: os.system("rm -r tmp"+jobname)
  except: print "No tmp/ directory"
  os.system("mkdir tmp"+jobname)
  try: os.stat("res"+jobname)
  except: os.mkdir("res"+jobname)


  ##### Creating and sending jobs #####
  joblist = []
  ###### loop for creating and sending jobs #####
  path = os.getcwd()
  for x in range(1, int(NumberOfJobs)+1):
     os.system("mkdir tmp"+jobname+"/"+str(files[x-1]).replace(".root",""))
     os.chdir("tmp"+jobname+"/"+str(files[x-1]).replace(".root",""))
     #os.system("mkdir tmp"+jobname+"/"+str(files[x-1]).replace(".root",""))
     os.chdir(path+"/tmp"+jobname+"/"+str(files[x-1]).replace(".root",""))

     with open('job_%s.sh'%files[x-1].replace(".root",""), 'w') as fout:
        fout.write("#!/bin/sh\n")
        fout.write("echo\n")
        fout.write("echo\n")
        fout.write("echo 'START---------------'\n")
        fout.write("echo 'WORKDIR ' ${PWD}\n")
        fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
        fout.write("cd "+str(path)+"\n")
        fout.write("cmsenv\n")
        fout.write(cmd+" "+files[x-1]+" "+path+"/res"+jobname+"/"+OutputFileNames+"_"+files[x-1]+" "+samples+"\n")
        print "EXECUTING: ",cmd+" "+files[x-1]+" "+path+"/res"+jobname+"/"+OutputFileNames+"_"+files[x-1]+" "+samples+"\n"
        fout.write("echo 'STOP---------------'\n")
        fout.write("echo\n")
        fout.write("echo\n")
     os.system("chmod 755 job_%s.sh"%(files[x-1].replace(".root","")) )

     if useCondorBatch:
       os.system("mv  job_*.sh "+jobname+".sh")
       makeSubmitFileCondor(jobname+".sh",jobname,"workday")
       os.system("condor_submit submit.sub")
     else:
       os.system("bsub -q "+queue+" -o logs job_%s.sh -J %s"%(files[x-1].replace(".root",""),jobname))
     print "job nr " + str(x) + " submitted: " + files[x-1].replace(".root","")
     joblist.append("%s"%(files[x-1].replace(".root","")))
     os.chdir("../..")

  print
  print "your jobs:"
  if useCondorBatch: os.system("condor_q")
  else: os.system("bjobs")
  userName=os.environ['USER']
  if wait: waitForBatchJobs(jobname,NumberOfJobs,NumberOfJobs, userName, timeCheck)

  print
  print 'END: makeData'
  print
  return joblist, files

# def mergeCPs(template,jobname="CPs"):
#
#   dir = "res"+jobname+"/"
#
#   print "Merging data from job " ,jobname
#   # read out files
#   resDir = 'res'+jobname+'/'
#   filelist = os.listdir('./'+resDir)
#
#   samples = {}
#   sampleTypes = template.split(',')
#   samples = {}
#   sampleTypes = template.split(',')
#   for t in sampleTypes:
#     print "For sample: ",t
#     files = []
#     for f in os.listdir(dir):
#       if f.find('.root') != -1 and f.find(t) != -1 and f.find("hadd") == -1:
#          print "Adding file = ",dir+f
#          files.append(dir+f)
#       samples[t] = files
#
#   vars =[]
#   ftest = ROOT.TFile(samples[sampleTypes[0]][0],"READ")
#   for h in ftest.GetListOfKeys():
#     vars.append(h.GetName())
#   ftest.Close()
#
#   print "Doing histograms for the following variables: "
#   print [v for v in vars]; print "" ;
#   # vars = ["looseSel_Dijet_invariant_mass"]
#
#   for t in sampleTypes:
#     print "For sample: " ,t
#     histlist = []
#     for v in vars:
#       print "For variable: " ,v
#       hZero = ROOT.TH1D()
#       for i,j in enumerate(samples[t]):
#         f1 = ROOT.TFile(j,"READ")
#         print "Opened file ", f1.GetName()
#         h1 = ROOT.TH1D(f1.Get(v))
#         h1.SetFillColor(36)
#         if i==0: hZero = copy.deepcopy(h1); continue
#         ROOT.gROOT.cd()
#         hnew = h1.Clone()
#         hZero.Add(hnew)
#       hClone   = copy.deepcopy(hZero)
#       histlist.append(hClone)
#
#     outf = ROOT.TFile.Open(dir+"/out_"+t+'.root','RECREATE')
#     outf.cd()
#     for h in histlist:
#       h.Write()
#     outf.Write()
#     outf.Close()

def mergeCPs(template,jobname="CPs"):
  
  dir = "res"+jobname+"/"

  print "Merging data from job " ,jobname
  # read out files
  resDir = 'res'+jobname+'/'

  sampleTypes = template.split(',')
  for t in sampleTypes:
    cmd = "hadd %s/%s.root %s/controlplots_*%s*.root" %(resDir,t,resDir,t)
    os.system(cmd)

def makePseudoDataVjets(input,kernel,mc,output,lumi,workspace,year,purity):
 
 pwd = os.getcwd()
 pwd = "/"
 ROOT.gRandom.SetSeed(0)
 
 finmc = ROOT.TFile.Open(pwd+'/'+input,'READ')
 hmcin = finmc.Get('nonRes')

 findata = ROOT.TFile.Open(pwd+'/'+kernel,'READ')
 #findata.ls()
 hdata = ROOT.TH3F()
 
 if   mc == 'pythia': hdata = findata.Get('histo')
 elif mc == 'herwig': hdata = findata.Get('histo_altshapeUp')
 elif mc == 'madgraph': hdata = findata.Get('histo_altshape2')
 elif mc == 'powheg': hdata = findata.Get('histo_NLO')
 
 fout = ROOT.TFile.Open(output,'RECREATE')
 #hmcin.Scale(10.)
 hmcin.Write('nonRes')
 xbins = array("f",getListOfBinsLowEdge(hmcin,"x"))
 zbins = array("f",getListOfBinsLowEdge(hmcin,"z"))
 hout = ROOT.TH3F('data','data',len(xbins)-1,xbins,len(xbins)-1,xbins,len(zbins)-1,zbins)
 
 hout.FillRandom(hdata,int(hmcin.Integral()*lumi))
 
 ws_file = ROOT.TFile.Open(workspace,'READ')
 ws = ws_file.Get('w')
 ws_file.Close()
 ws.Print()

 model = ws.pdf('shapeBkg_Vjet_JJ_%s_13TeV_%i'%(purity,year))
 category = ws.obj("CMS_channel==CMS_channel::JJ_"+purity+"_13TeV_%i"%year)

 MJ1= ws.var("MJ1");
 MJ2= ws.var("MJ2");
 MJJ= ws.var("MJJ");
 args = ROOT.RooArgSet(MJ1,MJ2,MJJ)
 print "n_exp_binJJ_"+purity+"_13TeV_%i_proc_Vjet"%year
 o_norm_vjets = ws.obj("n_exp_binJJ_"+purity+"_13TeV_%i_proc_Vjet"%year)
 hout_vjets = ROOT.TH3F('vjets','vjets',len(xbins)-1,xbins,len(xbins)-1,xbins,len(zbins)-1,zbins)
 
 nEvents = o_norm_vjets.getVal()
 print "Expected V+jets events: ",nEvents
 vjets = model.generate(args,int(nEvents))
 if vjets!=None:
  #print signal.sumEntries()
  for i in range(0,int(vjets.sumEntries())):
   a = vjets.get(i)
   it = a.createIterator()
   var = it.Next()
   x=[]
   while var:
       x.append(var.getVal())
       var = it.Next()
   #print x
   hout_vjets.Fill(x[0],x[1],x[2])
      
 hout.Add(hout_vjets)
 
 fout.cd()
 hout.Write('data')
 
 print "input    ", input
 print "kernel   ", kernel
 print "mc       ", mc
 print "output   ", output
 print "lumi     ", lumi
 print "workspace", workspace
 print "year     ", year
 print "purity   ", purity
 print "Expected V+jets events: ",nEvents
 print "Expected QCD events: ",int(hmcin.Integral()*lumi)
 print "Writing histograms nonRes and data to file ", output

 finmc.Close()
 findata.Close()
 fout.Close()    
