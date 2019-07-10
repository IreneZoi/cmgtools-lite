import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
import json
import sys


class DataCardMaker:
    def __init__(self,finalstate,category,period,luminosity=1.0,physics="LJ",cat=""):
        self.physics=physics
        self.finalstate=finalstate
        self.category=category
        self.period=period
        self.contributions=[]
        self.systematics=[]
        self.tag=self.physics+"_"+category+"_"+period
	self.cat=cat
	if self.cat!="": self.rootFile = ROOT.TFile("datacardInputs_"+cat+".root","RECREATE")
        else: self.rootFile = ROOT.TFile("datacardInputs_"+self.tag+".root","RECREATE")
        self.rootFile.cd()
        self.w=ROOT.RooWorkspace("w","w")
        self.luminosity=luminosity
        self.w.factory(self.physics+"_"+period+"_lumi["+str(luminosity)+"]")
        if period=='8TeV':
            self.sqrt_s=8000.0
        if period=='13TeV':
            self.sqrt_s=13000.0


    def addSystematic(self,name,kind,values,addPar = ""):
        self.systematics.append({'name':name,'kind':kind,'values':values })


    def addMVVSignalParametricShape(self,name,variable,jsonFile,scale ={},resolution={}):
        self.w.factory("MH[2000]")
        self.w.var("MH").setConstant(1)
       
        scaleStr='0'
        resolutionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
        for syst,factor in scale.iteritems():
            self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.iteritems():
            self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            resolutionSysts.append(syst)
       
        MVV=variable            
        self.w.factory(variable+"[0,13000]")
        self.w.factory("MJ1[0,13000]")
        self.w.factory("MJ2[0,13000]")

        
        f=open(jsonFile)
        info=json.load(f)

        SCALEVar="_".join(["MEAN",name,self.tag])
        SIGMAVar="_".join(["SIGMA",name,self.tag])
        if name.find("H")!=-1:
            self.w.factory("expr::"+SCALEVar+"('("+info['MEAN']+")*(1+"+scaleStr+")',{MH,MJ1,MJ2},"+','.join(scaleSysts)+")")
            self.w.factory("expr::"+SIGMAVar+"('("+info['SIGMA']+")*(1+"+resolutionStr+")',{MH,MJ1,MJ2},"+','.join(resolutionSysts)+")") 
        else:    
            self.w.factory("expr::"+SCALEVar+"('("+info['MEAN']+")*(1+"+scaleStr+")',MH,"+','.join(scaleSysts)+")")
            self.w.factory("expr::"+SIGMAVar+"('("+info['SIGMA']+")*(1+"+resolutionStr+")',MH,"+','.join(resolutionSysts)+")")
        
        

        ALPHA1Var="_".join(["ALPHA1",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=ALPHA1Var,param=info['ALPHA1']))

        ALPHA2Var="_".join(["ALPHA2",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=ALPHA2Var,param=info['ALPHA2']))

        N1Var="_".join(["N1",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=N1Var,param=info['N1']))

        N2Var="_".join(["N2",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=N2Var,param=info['N2']))        

        pdfName="_".join([name,self.tag])
        vvMass = ROOT.RooDoubleCB(pdfName,pdfName,self.w.var(MVV),self.w.function(SCALEVar),self.w.function(SIGMAVar),self.w.function(ALPHA1Var),self.w.function(N1Var),self.w.function(ALPHA2Var),self.w.function(N2Var))
 	getattr(self.w,'import')(vvMass,ROOT.RooFit.Rename(pdfName))
        f.close()
	
    def addMVVSignalParametricShape2(self,name,variable,jsonFile,scale ={},resolution={}):

        if self.w.var("MH") == None: self.w.factory("MH[2000]")
	self.w.var("MH").setVal(2000.)
        self.w.var("MH").setConstant(1)
       
        scaleStr='0'
        resolutionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
       
        for syst,factor in scale.iteritems():
            self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.iteritems():
            #self.w.factory(syst+"[0,-0.8,0.8]")
            self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
    
            resolutionSysts.append(syst)

        MVV=variable    
	if self.w.var(MVV) == None: self.w.factory(MVV+"[0,13000]")

        f=open(jsonFile)
        info=json.load(f)

        SCALEVar="_".join(["MEAN",name,self.tag])
        self.w.factory("expr::{name}('({param})*(1+{vv_syst})',MH,{vv_systs})".format(name=SCALEVar,param=info['MEAN'],vv_syst=scaleStr,vv_systs=','.join(scaleSysts)))

        SIGMAVar="_".join(["SIGMA",name,self.tag])
        self.w.factory("expr::{name}('({param})*(1+{vv_syst})',MH,{vv_systs})".format(name=SIGMAVar,param=info['SIGMA'],vv_syst=resolutionStr,vv_systs=','.join(resolutionSysts)))
        print SIGMAVar

        ALPHAVar="_".join(["ALPHA",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=ALPHAVar,param=info['ALPHA']))

        SCALESIGMAVar="_".join(["SCALESIGMA",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=SCALESIGMAVar,param=info['SCALESIGMA']))

        NVar="_".join(["N",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=NVar,param=info['N']))

        fVar="_".join(["f",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=fVar,param=info['f']))        

        GSIGMAVar="_".join(["GSIGMA",name,self.tag])
        self.w.factory("prod::{name}({param1},{param2})".format(name=GSIGMAVar,param1=SIGMAVar,param2=SCALESIGMAVar))        

        self.w.var(MVV).Print()
	self.w.function(SCALEVar).Print()
	self.w.function(SIGMAVar).Print()
	self.w.function(GSIGMAVar).Print()

        pdfName="_".join([name,self.tag])
	gaus = ROOT.RooGaussian(pdfName+'_Gaus',pdfName+'_Gaus',self.w.var(MVV),self.w.function(SCALEVar),self.w.function(GSIGMAVar))
	cb = ROOT.RooCBShape(pdfName+'_CB',pdfName+'_CB',self.w.var(MVV),self.w.function(SCALEVar),self.w.function(SIGMAVar),self.w.function(ALPHAVar),self.w.function(NVar))
        vvMass = ROOT.RooAddPdf(pdfName,pdfName,gaus,cb,self.w.function(fVar))
        getattr(self.w,'import')(vvMass,ROOT.RooFit.Rename(pdfName))
        f.close()
        
        
        
    def addGaussianShape(self,name,variable,parameters,newTag=""):
        print parameters
        Mjet=variable
        if self.w.var(Mjet) == None: self.w.factory(Mjet+"[0,10000]")
        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag
        print 'mean'
        mean="_".join(["mean",tag])
        if "mean" in parameters.keys():
            val = parameters['mean']['val']
            ranges = val
        else:
            val = 15.0
            print "attention set value to default in addMjetBackgroundShapeVJets"
        self.w.factory("{name}[{val},{ranges},{ranges}]".format(name=mean,val=val,ranges=ranges))
	self.w.var(mean).setConstant(1)
	print 'sigma'
	sigma="_".join(["sigma",tag])
        if "sigma" in parameters.keys():
            val = parameters['sigma']['val']
            ranges = val
        else:
            val = 15.0
            print "attention set value to default in addMjetBackgroundShapeVJets"
        self.w.factory("{name}[{val},{ranges},{ranges}]".format(name=sigma,val=val,ranges=ranges))
	self.w.var(sigma).setConstant(1)
        print 'function '
        pdfName="_".join([name,self.tag])
	gaus = ROOT.RooGaussian(pdfName,pdfName,self.w.var(Mjet),self.w.var(mean),self.w.var(sigma))
        getattr(self.w,'import')(gaus,ROOT.RooFit.Rename(pdfName))





    def addMJJSignalParametricShape(self,name,variable,jsonFile,scale ={},resolution={},scales=[1,1],varToReplace="MH"):
        self.w.factory("MH[2000]")
        self.w.var("MH").setConstant(1)
       
        scaleStr='0'
        resolutionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
        for syst,factor in scale.iteritems():
            self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.iteritems():
            self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            resolutionSysts.append(syst)
       
        MJJ=variable            
        self.w.factory(variable+"[0,1000]")

        
        f=open(jsonFile)
        info=json.load(f)

        SCALEVar="_".join(["mean",name,self.tag])
        self.w.factory("expr::{name}('({param}*{scale})*(1+{vv_syst})',MH,{vv_systs})".format(name=SCALEVar,param=info['mean'],scale=str(scales[0]),vv_syst=scaleStr,vv_systs=','.join(scaleSysts)).replace("MH",varToReplace))

        SIGMAVar="_".join(["sigma",name,self.tag])
        self.w.factory("expr::{name}('({param}*{res})*(1+{vv_syst})',MH,{vv_systs})".format(name=SIGMAVar,param=info['sigma'],res=str(scales[1]),vv_syst=resolutionStr,vv_systs=','.join(resolutionSysts)).replace("MH",varToReplace))

        ALPHAVar="_".join(["alpha",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=ALPHAVar,param=info['alpha']).replace("MH",varToReplace))

        NVar="_".join(["n",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=NVar,param=info['n']).replace("MH",varToReplace))

        ALPHAVar2="_".join(["alpha2",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=ALPHAVar2,param=info['alpha2']).replace("MH",varToReplace))

        NVar2="_".join(["n2",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=NVar2,param=info['n2']).replace("MH",varToReplace))

        pdfName="_".join([name+"peak",self.tag])
        vvMass = ROOT.RooDoubleCB(pdfName,pdfName,self.w.var(MJJ),self.w.function(SCALEVar),self.w.function(SIGMAVar),self.w.function(ALPHAVar),self.w.function(NVar),self.w.function(ALPHAVar2),self.w.function(NVar2))
        getattr(self.w,'import')(vvMass,ROOT.RooFit.Rename(pdfName))


        SLOPEVar="_".join(["slope",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=SLOPEVar,param=info['slope']).replace("MH",varToReplace))

        FVar="_".join(["f",name,self.tag])
        self.w.factory("expr::{name}('min(MH*0+{param},1)',MH)".format(name=FVar,param=info['f']).replace("MH",varToReplace))

        pdfName2="_".join([name+"bkg",self.tag])
        
        self.w.factory("RooExponential::{name}({var},{SLOPE})".format(name=pdfName2,var=MJJ,SLOPE=SLOPEVar).replace("MH",varToReplace))

        pdfName3="_".join([name,self.tag])
        self.w.factory("SUM::{name}({f}*{PDF1},{PDF2})".format(name=pdfName3,f=FVar,PDF1=pdfName,PDF2=pdfName2))

        f.close()




    def addMJJTopMergedParametricShape(self,name,variable,jsonFile,scale ={},resolution={},fraction={},varToReplace="MH"):
        self.w.factory("MH[2000]")
        self.w.var("MH").setConstant(1)
       
        scaleStr='0'
        resolutionStr='0'
        fractionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
        fractionSysts=[]
        for syst,factor in scale.iteritems():
            self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.iteritems():
            self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            resolutionSysts.append(syst)
        for syst,factor in fraction.iteritems():
            self.w.factory(syst+"[0,-50,50]")
            fractionStr=fractionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            fractionSysts.append(syst)
       
        MJJ=variable            
        self.w.factory(variable+"[0,1000]")

        
        f=open(jsonFile)
        info=json.load(f)

        MEANWVar="_".join(["meanW",name,self.tag])
        self.w.factory("expr::{name}('({param})*(1+{vv_syst})',MH,{vv_systs})".format(name=MEANWVar,param=info['meanW'],vv_syst=scaleStr,vv_systs=','.join(scaleSysts)).replace("MH",varToReplace))

        MEANTOPVar="_".join(["meanTop",name,self.tag])
        self.w.factory("expr::{name}('({param})*(1+{vv_syst})',MH,{vv_systs})".format(name=MEANTOPVar,param=info['meanTop'],vv_syst=scaleStr,vv_systs=','.join(scaleSysts)).replace("MH",varToReplace))

        SIGMAWVar="_".join(["sigmaW",name,self.tag])
        self.w.factory("expr::{name}('({param})*(1+{vv_syst})',MH,{vv_systs})".format(name=SIGMAWVar,param=info['sigmaW'],vv_syst=resolutionStr,vv_systs=','.join(resolutionSysts)).replace("MH",varToReplace))

        SIGMATOPVar="_".join(["sigmaTop",name,self.tag])
        self.w.factory("expr::{name}('({param})*(1+{vv_syst})',MH,{vv_systs})".format(name=SIGMATOPVar,param=info['sigmaTop'],vv_syst=resolutionStr,vv_systs=','.join(resolutionSysts)).replace("MH",varToReplace))

        ALPHAWVar="_".join(["alphaW",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=ALPHAWVar,param=info['alphaW']).replace("MH",varToReplace))

        ALPHAWVar2="_".join(["alphaW2",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=ALPHAWVar2,param=info['alphaW2']).replace("MH",varToReplace))

        ALPHATOPVar="_".join(["alphaTop",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=ALPHATOPVar,param=info['alphaTop']).replace("MH",varToReplace))

        ALPHATOPVar2="_".join(["alphaTop2",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=ALPHATOPVar2,param=info['alphaTop2']).replace("MH",varToReplace))

        NVar="_".join(["n",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=NVar,param=info['n']).replace("MH",varToReplace))


        FVar="_".join(["f",name,self.tag])
        self.w.factory("expr::{name}('min((MH*0+{param})*(1+{vv_syst}),1.0)',MH,{vv_systs})".format(name=FVar,param=info['f'],vv_syst=fractionStr,vv_systs=','.join(fractionSysts)).replace("MH",varToReplace))


        pdfNameW="_".join([name+"PeakW",self.tag])
        vvMass = ROOT.RooDoubleCB(pdfNameW,pdfNameW,self.w.var(MJJ),self.w.function(MEANWVar),self.w.function(SIGMAWVar),self.w.function(ALPHAWVar),self.w.function(NVar),self.w.function(ALPHAWVar2),self.w.function(NVar))
        getattr(self.w,'import')(vvMass,ROOT.RooFit.Rename(pdfNameW))

        pdfNameTop="_".join([name+"PeakTop",self.tag])
        vvMass2 = ROOT.RooDoubleCB(pdfNameTop,pdfNameTop,self.w.var(MJJ),self.w.function(MEANTOPVar),self.w.function(SIGMATOPVar),self.w.function(ALPHATOPVar),self.w.function(NVar),self.w.function(ALPHATOPVar2),self.w.function(NVar))
        getattr(self.w,'import')(vvMass2,ROOT.RooFit.Rename(pdfNameTop))

        pdfNamePeak="_".join([name+"Peak",self.tag])
        self.w.factory("SUM::{name}({f}*{PDF1},{PDF2})".format(name=pdfNamePeak,f=FVar,PDF1=pdfNameW,PDF2=pdfNameTop))


        SLOPEVar="_".join(["slope",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=SLOPEVar,param=info['slope']).replace("MH",varToReplace))

        pdfNameBKG="_".join([name+"bkg",self.tag])
        self.w.factory("RooExponential::{name}({var},{SLOPE})".format(name=pdfNameBKG,var=MJJ,SLOPE=SLOPEVar).replace("MH",varToReplace))

        FVar2="_".join(["f2",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=FVar2,param=info['f2']).replace("MH",varToReplace))

        pdfName="_".join([name,self.tag])
        self.w.factory("SUM::{name}({f}*{PDF1},{PDF2})".format(name=pdfName,f=FVar2,PDF1=pdfNamePeak,PDF2=pdfNameBKG))
        f.close()



   
    def addMJJSignalParametricShapeNOEXP(self,name,variable,jsonFile,scale ={},resolution={},scales=[1,1],varToReplace="MH"):

        if self.w.var("MH") == None: self.w.factory("MH[2000]")
	self.w.var("MH").setVal(2000.)
        self.w.var("MH").setConstant(1)
       
        scaleStr='0'
        resolutionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
        for syst,factor in scale.iteritems():
            if self.w.var(syst) != None: continue
            self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.iteritems():
            if self.w.var(syst) != None: continue
            self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)

            resolutionSysts.append(syst)
       
        MJJ=variable            
        if self.w.var(MJJ) == None: self.w.factory(MJJ+"[0,1000]")

        f=open(jsonFile)
        info=json.load(f)

        SCALEVar="_".join(["mean",name,self.tag])
        self.w.factory("expr::{name}('({param}*{sc})*(1+{vv_syst})',MH,{vv_systs})".format(name=SCALEVar,param=info['mean'],sc=scales[0],vv_syst=scaleStr,vv_systs=','.join(scaleSysts)).replace("MH",varToReplace))

        
        SIGMAVar="_".join(["sigma",name,self.tag])
        self.w.factory("expr::{name}('({param}*{res})*(1+{vv_syst})',MH,{vv_systs})".format(name=SIGMAVar,param=info['sigma'],res=scales[1],vv_syst=resolutionStr,vv_systs=','.join(resolutionSysts)).replace("MH",varToReplace))

        ALPHAVar="_".join(["alpha",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=ALPHAVar,param=info['alpha']).replace("MH",varToReplace))

        NVar="_".join(["n",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=NVar,param=info['n']).replace("MH",varToReplace))

        ALPHAVar2="_".join(["alpha2",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=ALPHAVar2,param=info['alpha2']).replace("MH",varToReplace))

        NVar2="_".join(["n2",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=NVar2,param=info['n2']).replace("MH",varToReplace))

        pdfName="_".join([name,self.tag])
        vvMass = ROOT.RooDoubleCB(pdfName,pdfName,self.w.var(MJJ),self.w.function(SCALEVar),self.w.function(SIGMAVar),self.w.function(ALPHAVar),self.w.function(NVar),self.w.function(ALPHAVar2),self.w.function(NVar2))
        getattr(self.w,'import')(vvMass,ROOT.RooFit.Rename(pdfName))

        f.close()
        
    def addMJJSignalParametricShapeHiggs(self,name,variable,jsonFile,scale ={},resolution={},scales=[1,1],varToReplace="MH"):

        if self.w.var("MH") == None: self.w.factory("MH[2000]")
	self.w.var("MH").setVal(2000.)
        self.w.var("MH").setConstant(1)
       
        scaleStr='0'
        resolutionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
        for syst,factor in scale.iteritems():
            if self.w.var(syst) != None: continue
            self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.iteritems():
            if self.w.var(syst) != None: continue
            self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)

            resolutionSysts.append(syst)
       
        MJJ=variable            
        if self.w.var(MJJ) == None: self.w.factory(MJJ+"[0,1000]")

        f=open(jsonFile)
        info=json.load(f)
        print "attention : "
        print "meanH "+str(info['meanH'])
        print variable
        print jsonFile
        SCALEVar="_".join(["meanH",name,self.tag])
        self.w.factory("expr::{name}('({param}*{sc})*(1+{vv_syst})',MH,{vv_systs})".format(name=SCALEVar,param=info['meanH'],sc=scales[0],vv_syst=scaleStr,vv_systs=','.join(scaleSysts)).replace("MH",varToReplace))

        
        SIGMAVar="_".join(["sigmaH",name,self.tag])
        self.w.factory("expr::{name}('({param}*{res})*(1+{vv_syst})',MH,{vv_systs})".format(name=SIGMAVar,param=info['sigmaH'],res=scales[1],vv_syst=resolutionStr,vv_systs=','.join(resolutionSysts)).replace("MH",varToReplace))

        ALPHAVar="_".join(["alphaH",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=ALPHAVar,param=info['alphaH']).replace("MH",varToReplace))

        NVar="_".join(["nH",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=NVar,param=info['nH']).replace("MH",varToReplace))

        ALPHAVar2="_".join(["alpha2H",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=ALPHAVar2,param=info['alpha2H']).replace("MH",varToReplace))

        NVar2="_".join(["n2H",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(name=NVar2,param=info['n2H']).replace("MH",varToReplace))

        pdfName="_".join([name,self.tag])
        vvMass = ROOT.RooDoubleCB(pdfName,pdfName,self.w.var(MJJ),self.w.function(SCALEVar),self.w.function(SIGMAVar),self.w.function(ALPHAVar),self.w.function(NVar),self.w.function(ALPHAVar2),self.w.function(NVar2))
        getattr(self.w,'import')(vvMass,ROOT.RooFit.Rename(pdfName))

        f.close()    

    def addMJJSignalShapeNOEXP(self,name,variable,newTag="",preconstrains={},scale ={},resolution={},scales=[1,1]):
       
        scaleStr='0'
        resolutionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
        for syst,factor in scale.iteritems():
            #if self.w.var(syst) != None: continue
            self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.iteritems():
            #if self.w.var(syst) != None: continue
            self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            resolutionSysts.append(syst)
       
        MJJ=variable            
        if self.w.var(MJJ) == None: self.w.factory(MJJ+"[0,1000]")

        SCALEVar="_".join(["mean",name,self.tag])
        self.w.factory("expr::{name}('({param}*{sc})*(1+{vv_syst})',{vv_systs})".format(name=SCALEVar,param=preconstrains['mean']['val'],sc=scales[0],vv_syst=scaleStr,vv_systs=','.join(scaleSysts)))

        SIGMAVar="_".join(["sigma",name,self.tag])
        self.w.factory("expr::{name}('({param}*{res})*(1+{vv_syst})',{vv_systs})".format(name=SIGMAVar,param=preconstrains['sigma']['val'],res=scales[1],vv_syst=resolutionStr,vv_systs=','.join(resolutionSysts)))

        ALPHAVar="_".join(["alpha",name,self.tag])
        self.w.factory("{name}[{param}]".format(name=ALPHAVar,param=preconstrains['alpha']['val']))

        NVar="_".join(["n",name,self.tag])
        self.w.factory("{name}[{param}]".format(name=NVar,param=preconstrains['n']['val']))

        ALPHAVar2="_".join(["alpha2",name,self.tag])
        self.w.factory("{name}[{param}]".format(name=ALPHAVar2,param=preconstrains['alpha2']['val']))

        NVar2="_".join(["n2",name,self.tag])
        self.w.factory("{name}[{param}]".format(name=NVar2,param=preconstrains['n2']['val']))

        pdfName="_".join([name,self.tag])
        vvMass = ROOT.RooDoubleCB(pdfName,pdfName,self.w.var(MJJ),self.w.function(SCALEVar),self.w.function(SIGMAVar),self.w.function(ALPHAVar),self.w.function(NVar),self.w.function(ALPHAVar2),self.w.function(NVar2))
        getattr(self.w,'import')(vvMass,ROOT.RooFit.Rename(pdfName))

    def addHistoShapeFromFile(self,name,observables,filename,histoname,systematics=[],conditional = False,order=0,newTag=""):     
        varset=ROOT.RooArgSet()
        varlist=ROOT.RooArgList()
        varPointers=[]
        for var in observables:
            if self.w.var(var) == None: self.w.factory(var+"[0,10000]")
            varPointers.append(self.w.var(var))
            varset.add(self.w.var(var))
            varlist.add(self.w.var(var))

        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag
        print "########################################"
        print tag
        print self.tag
        print "#########################################"
        FR=ROOT.TFile(filename)

        #Load PDF
        histo=FR.Get(histoname)
	print histo.GetName()


        if len(systematics)>0:
            histName="_".join([name+"NominalHIST",tag])
            pdfName="_".join([name+"Nominal",self.tag])
        else:
            histName="_".join([name+"HIST",tag])
            pdfName="_".join([name,self.tag])

        roohist = ROOT.RooDataHist(histName,histName,varlist,histo)      
	varset.Print()
	varlist.Print()
        pdf=ROOT.RooHistPdf(pdfName,pdfName,varset,roohist,order)
        if self.w.data(histName) == None: getattr(self.w,'import')(roohist,ROOT.RooFit.Rename(histName))
        if self.w.pdf(pdfName) == None: getattr(self.w,'import')(pdf,ROOT.RooFit.Rename(pdfName))
        #Load SYstematics
        coeffList=ROOT.RooArgList()
        pdfList=ROOT.RooArgList(self.w.pdf(pdfName))

        for systval in systematics:
            splitted=systval.split(':')
            systName=splitted[1]
            syst=splitted[0]
            if self.w.var(systName) == None: self.w.factory(systName+"[-1,1]")
            coeffList.add(self.w.var(systName))

            for variation in ["Up","Down"]:
                histo=FR.Get(histoname+"_"+syst+variation)
                print 'loaded',histoname+"_"+syst+variation
                histName="_".join([name+"_"+syst+variation+"HIST",tag])
                roohist = ROOT.RooDataHist(histName,histName,varlist,histo)
       
                pdfName="_".join([name+"_"+syst+variation,self.tag])
                pdf=ROOT.RooHistPdf(pdfName,pdfName,varset,roohist,order)

                if self.w.data(histName) == None: getattr(self.w,'import')(roohist,ROOT.RooFit.Rename(histName))
                if self.w.pdf(pdfName) == None: getattr(self.w,'import')(pdf,ROOT.RooFit.Rename(pdfName))
                pdfList.add(self.w.pdf(pdfName))

        pdfName="_".join([name,self.tag])
        if len(systematics)>0:
            if len(observables)==1:
                total=ROOT.FastVerticalInterpHistPdf(pdfName,pdfName,self.w.var(observables[0]),pdfList, coeffList)
            elif len(observables)==2:
                total=ROOT.FastVerticalInterpHistPdf2D(pdfName,pdfName,self.w.var(observables[0]),self.w.var(observables[1]),conditional,pdfList,coeffList,1,-1)
            elif len(observables)==3:
		total=ROOT.FastVerticalInterpHistPdf3D(pdfName,pdfName,self.w.var(observables[0]),self.w.var(observables[1]),self.w.var(observables[2]),conditional,pdfList,coeffList,1,-1)
	    
	    if self.w.pdf(pdfName) == None: getattr(self.w,'import')(total,ROOT.RooFit.Rename(pdfName))


    def addShapes(self,name,observables,filename,histoname,systematics=[],conditional = False,order=0,newTag="",filename2="",systematics2=[]):     
        varset=ROOT.RooArgSet()
        varlist=ROOT.RooArgList()
        varPointers=[]
        for var in observables:
            if self.w.var(var) == None: self.w.factory(var+"[0,10000]")
            varPointers.append(self.w.var(var))
            varset.add(self.w.var(var))
            varlist.add(self.w.var(var))

        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag
        print "########################################"
        print tag
        print self.tag
        print "#########################################"
        FR=ROOT.TFile(filename)
        
        #Load PDF
        histo=FR.Get(histoname)
	print histo.GetName()


        if len(systematics)>0:
            histName="_".join([name+"NominalHIST",tag])
            pdfName="_".join([name+"Nominal",self.tag])
        else:
            histName="_".join([name+"HIST",tag])
            pdfName="_".join([name,self.tag])

        roohist = ROOT.RooDataHist(histName,histName,varlist,histo)      
	varset.Print()
	varlist.Print()
        pdf=ROOT.RooHistPdf(pdfName,pdfName,varset,roohist,order)
        getattr(self.w,'import')(roohist,ROOT.RooFit.Rename(histName))
        getattr(self.w,'import')(pdf,ROOT.RooFit.Rename(pdfName))
        #Load SYstematics
        coeffList=ROOT.RooArgList()
        pdfList=ROOT.RooArgList(self.w.pdf(pdfName))

        for systval in systematics:
            splitted=systval.split(':')
            systName=splitted[1]
            syst=splitted[0]
            self.w.factory(systName+"[-1,1]")
            coeffList.add(self.w.var(systName))

            for variation in ["Up","Down"]:
                histo=FR.Get(histoname+"_"+syst+variation)
                print 'loaded',histoname+"_"+syst+variation
                histName="_".join([name+"_"+syst+variation+"HIST",tag])
                roohist = ROOT.RooDataHist(histName,histName,varlist,histo)
       
                pdfName="_".join([name+"_"+syst+variation,self.tag])
                pdf=ROOT.RooHistPdf(pdfName,pdfName,varset,roohist,order)

                getattr(self.w,'import')(roohist,ROOT.RooFit.Rename(histName))
                getattr(self.w,'import')(pdf,ROOT.RooFit.Rename(pdfName))
                pdfList.add(self.w.pdf(pdfName))
         
        
        if len(systematics2)>0:
            FR2 = ROOT.TFile(filename2) 

            for systval in systematics2:
                splitted=systval.split(':')
                systName=splitted[1]
                syst=splitted[0]
                self.w.factory(systName+"[-1,1]")
                coeffList.add(self.w.var(systName))

                for variation in ["Up","Down"]:
                    histo=FR2.Get(histoname+"_"+syst+variation)
                    print 'loaded',histoname+"_"+syst+variation
                    histName="_".join([name+"_"+syst+variation+"HIST",tag])
                    roohist = ROOT.RooDataHist(histName,histName,varlist,histo)
        
                    pdfName="_".join([name+"_"+syst+variation,self.tag])
                    pdf=ROOT.RooHistPdf(pdfName,pdfName,varset,roohist,order)

                    getattr(self.w,'import')(roohist,ROOT.RooFit.Rename(histName))
                    getattr(self.w,'import')(pdf,ROOT.RooFit.Rename(pdfName))
                    pdfList.add(self.w.pdf(pdfName))
                
        

        pdfName="_".join([name,self.tag])
        if len(systematics)>0:
            total=ROOT.FastVerticalInterpHistPdf3D(pdfName,pdfName,self.w.var(observables[0]),self.w.var(observables[1]),self.w.var(observables[2]),conditional,pdfList,coeffList,1,-1)
	    getattr(self.w,'import')(total,ROOT.RooFit.Rename(pdfName))




    def addMJJParametricBackgroundShapeErfExp(self,name,variable,jsonFile,systP0={},systP1={},systP2={}):

        MJJ=variable
        self.w.factory(MJJ+"[0,10000]")
        f=open(jsonFile)
        info=json.load(f)


        p0Systs=[]
        p1Systs=[]
        p2Systs=[]

        p0SystStr='0'
        p1SystStr='0'
        p2SystStr='0'

        for syst,factor in systP0.iteritems():
            self.w.factory(syst+"[0,-0.5,0.5]")
            p0SystStr+="+{factor}*{syst}".format(factor=factor,syst=syst)
            p0Systs.append(syst)
        for syst,factor in systP1.iteritems():
            self.w.factory(syst+"[0,-0.5,0.5]")
            p1SystStr+="+{factor}*{syst}".format(factor=factor,syst=syst)
            p1Systs.append(syst)
        for syst,factor in systP2.iteritems():
            self.w.factory(syst+"[0,-0.5,0.5]")
            p2SystStr+="+{factor}*{syst}".format(factor=factor,syst=syst)
            p2Systs.append(syst)



        p0="_".join(["p0",name,self.tag])
        self.w.factory("expr::{name}('({param})*(1+{vv_syst})',{vv_systs})".format(name=p0,param=info['c_0'],vv_syst=p0SystStr,vv_systs=','.join(p0Systs)))

        p1="_".join(["p1",name,self.tag])
        self.w.factory("expr::{name}('({param})*(1+{vv_syst})',{vv_systs})".format(name=p1,param=info['c_1'],vv_syst=p1SystStr,vv_systs=','.join(p1Systs)))

        p2="_".join(["p2",name,self.tag])
        self.w.factory("expr::{name}('({param})*(1+{vv_syst})',{vv_systs})".format(name=p2,param=info['c_2'],vv_syst=p2SystStr,vv_systs=','.join(p2Systs)))

        pdfName="_".join([name,self.tag])
        erfexp = ROOT.RooErfExpPdf(pdfName,pdfName,self.w.var(MJJ),self.w.function(p0),self.w.function(p1),self.w.function(p2))
        getattr(self.w,'import')(erfexp,ROOT.RooFit.Rename(pdfName))
        f.close()






    def addMJJParametricBackgroundShapeExpo(self,name,variable,jsonFile,systP0={}):

        MJJ=variable
        self.w.factory(MJJ+"[0,10000]")

        f=open(jsonFile)
        info=json.load(f)

        p0Systs=[]
        p0SystStr='0'
        for syst,factor in systP0.iteritems():
            self.w.factory(syst+"[0,-0.5,0.5]")
            p0SystStr+="+{factor}*{syst}".format(factor=factor,syst=syst)
            p0Systs.append(syst)

        p0="_".join(["p0",name,self.tag])
        self.w.factory("expr::{name}('({param})*(1+{vv_syst})',{vv_systs})".format(name=p0,param=info['c_0'],vv_syst=p0SystStr,vv_systs=','.join(p0Systs)))


        pdfName="_".join([name,self.tag])
        self.w.factory("RooExponential::{name}({x},{slope})".format(name=pdfName,x=MJJ,slope=p0))
        f.close()


    def addMJJFloatingBackgroundShapeErfExp(self,name,variable,newTag=""):
        MJJ=variable
        self.w.factory(MJJ+"[0,1000]")

        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag


        p0="_".join(["c_0",tag])
        self.w.factory(p0+"[-0.03,-10,10]")
        p1="_".join(["c_1",tag])
        self.w.factory(p1+"[50,40,160]")
        p2="_".join(["c_2",tag])
        self.w.factory(p2+"[20,1,200]")

        pdfName="_".join([name,self.tag])
        bernsteinPDF = ROOT.RooErfExpPdf(pdfName,pdfName,self.w.var(MJJ),self.w.var(p0),self.w.var(p1),self.w.var(p2))
        getattr(self.w,'import')(bernsteinPDF,ROOT.RooFit.Rename(pdfName))

    def addMJJFloatingBackgroundShapeBifur(self,name,variable,newTag=""):
        MJJ=variable
        self.w.factory(MJJ+"[0,1000]")

        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag


        p0="_".join(["c_0",tag])
        self.w.factory(p0+"[60,40,160]")
        p1="_".join(["c_1",tag])
        self.w.factory(p1+"[20,1,200]")
        p2="_".join(["c_2",tag])
        self.w.factory(p2+"[20,1,200]")


        pdfName="_".join([name,self.tag])
        self.w.factory("RooBifurGauss::{name}({var},{p0},{p1},{p2})".format(name=pdfName,var=MJJ,p0=p0,p1=p1,p2=p2))



    def addMJJFloatingBackgroundShapeExpo(self,name,variable,newTag=""):
        MJJ=variable
        self.w.factory(MJJ+"[0,10000]")

        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag


        p0="_".join(["c_0",tag])

        self.w.factory(p0+"[-10,0]")

        pdfName="_".join([name,self.tag])
        self.w.factory("RooExponential::{name}({var},{p0})".format(name=pdfName,var=MJJ,p0=p0))





    def addMVVBackgroundShapeQCD(self,name,variable,logTerm=False,newTag="",preconstrains={}):

        MVV=variable
        if self.w.var(MVV) == None: self.w.factory(MVV+"[0,10000]")


        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag

        p0="_".join(["CMS_VV_JJ_p0",tag])
	self.w.factory("{name}[{val},10,60]".format(name=p0,val=preconstrains['CMS_p0']['val']))
       
        p1="_".join(["CMS_VV_JJ_p1",tag])
	self.w.factory("{name}[{val},0,5]".format(name=p1,val=preconstrains['CMS_p1']['val']))

        p2="_".join(["CMS_VV_JJ_p2",tag])
        
        if logTerm: self.w.factory("{name}[{val},0,10]".format(name=p2,val=preconstrains['CMS_p2']['val']))
        
        pdfName="_".join([name,self.tag])
        qcd = ROOT.RooQCDPdf(pdfName,pdfName,self.w.var(MVV),self.w.var(p0),self.w.var(p1),self.w.var(p2))
        getattr(self.w,'import')(qcd,ROOT.RooFit.Rename(pdfName))
        
        
    def addMjetBackgroundShapeVJets(self,name,variable,newTag="",preconstrains={}):
        Mjet=variable
        self.w.factory(Mjet+"[0,10000]")
        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag

        p0="_".join(["c_0",tag])
        if "c_0" in preconstrains.keys():
            val = preconstrains['c_0']['val']
            err = preconstrains['c_0']['err']
            #self.addSystematic(p0,"param",[val,err])
        else:
            val = 15.0
            print "attention set value to default in addMjetBackgroundShapeVJets"
        self.w.factory("{name}[{val},-{err},{err}]".format(name=p0,val=val,err=err))
        p1="_".join(["c_1",tag])
        if "c_1" in preconstrains.keys():
            val = preconstrains['c_1']['val']
            err = preconstrains['c_1']['err']
            #self.addSystematic(p1,"param",[val,err])
        else:
            val = 0.001
            print "attention set value to default in addMjetBackgroundShapeVJets"
        self.w.factory("{name}[{val},-{err},{err}]".format(name=p1,val=val,err=err))

        p2="_".join(["c_2",tag])
        if "c_2" in preconstrains.keys():
            val = preconstrains['c_2']['val']
            err = preconstrains['c_2']['err']
            #self.addSystematic(p2,"param",[val,err])
        else:
            val = 0.001
            print "attention set value to default in addMjetBackgroundShapeVJets"
        self.w.factory("{name}[{val},-{err},{err}]".format(name=p2,val=val,err=err))
        pdfName="_".join([name,self.tag])
        wjet = ROOT.RooErfPowPdf(pdfName,pdfName,self.w.var(Mjet),self.w.var(p0),self.w.var(p1),self.w.var(p2))
        getattr(self.w,'import')(wjet,ROOT.RooFit.Rename(pdfName))


    def addMjetBackgroundShapeVJetsGaus(self,name,variable,newTag="",preconstrains={},scale ={},resolution={}):
        Mjet=variable
        if self.w.var(Mjet)==None:
            self.w.factory(Mjet+"[0,10000]")
        
        scaleStr='0'
        resolutionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
        for syst,factor in scale.iteritems():
            #if self.w.var(syst) == None:
	    self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.iteritems():
            #if self.w.var(syst) == None:
	    self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            resolutionSysts.append(syst)

        if newTag !="": tag=newTag
        else: tag=name+"_"+self.tag

        mean="_".join(["mean",tag])
        self.w.factory("expr::{name}('{param}*(1+{vv_syst})',{vv_systs})".format(name=mean,param=preconstrains['mean']['val'],vv_syst=scaleStr,vv_systs=','.join(scaleSysts)))
                
        sigma="_".join(["sigma",tag])
        self.w.factory("expr::{name}('({param})*(1+{vv_syst})',{vv_systs})".format(name=sigma,param=preconstrains['sigma']['val'],vv_syst=resolutionStr,vv_systs=','.join(resolutionSysts)))
        
        pdfName="_".join([name,self.tag])
        wjet = ROOT.RooGaussian(pdfName,pdfName,self.w.var(Mjet),self.w.function(mean),self.w.function(sigma))
        getattr(self.w,'import')(wjet,ROOT.RooFit.Rename(pdfName))

    def addMjetBackgroundShapeVJetsRes(self,name,variable,newTag="",preconstrains={},scale ={},resolution={}):
        print "start importing shapes for resonant Vjets backgorund "+str(name)
        
        Mjet=variable
        self.w.var(Mjet)
        if self.w.var(Mjet)==None:
            self.w.factory(Mjet+"[0,10000]")
        
        scaleStr='0'
        resolutionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
        for syst,factor in scale.iteritems():
            self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.iteritems():
            self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            resolutionSysts.append(syst)

        if newTag !="": tag=newTag
        else: tag=name+"_"+self.tag

        mean="_".join(["mean",tag])
        if len(scaleSysts)>=1:
            self.w.factory("expr::{name}('{param}*(1+{vv_syst})',{vv_systs})".format(name=mean,param=preconstrains['mean']['val'],vv_syst=scaleStr,vv_systs=','.join(scaleSysts)))
        else:
            self.w.factory("{name}[{val},{err},{err}]".format(name=mean,val=preconstrains['mean']['val'],err=preconstrains['mean']['val']))
            self.w.var(mean).setConstant(1)
                
        sigma="_".join(["sigma",tag])
        if len(resolutionSysts)>=1:
            self.w.factory("expr::{name}('({param})*(1+{vv_syst})',{vv_systs})".format(name=sigma,param=preconstrains['sigma']['val'],vv_syst=resolutionStr,vv_systs=','.join(resolutionSysts)))
        else:
            self.w.factory("{name}[{val},{err},{err}]".format(name=sigma,val=preconstrains['sigma']['val'],err=preconstrains['sigma']['val']))
            self.w.var(sigma).setConstant(1)
        
        alpha="_".join(["alpha",tag])

        if "alpha" in preconstrains.keys():
            val = preconstrains['alpha']['val']
            err = 0
        else:
            val = 15.0
            print "attention set value to default in addMjetBackgroundShapeVJets"

        self.w.factory("{name}[{val},{err},{err}]".format(name=alpha,val=val,err=val))
        self.w.var(alpha).setConstant(1)
        
        n="_".join(["n",tag])

        if "n" in preconstrains.keys():
            val = preconstrains['n']['val']
            err = 0
        else:
            val = 15.0
            print "attention set value to default in addMjetBackgroundShapeVJets"
        self.w.factory("{name}[{val},{err},{err}]".format(name=n,val=val,err=val))
	self.w.var(n).setConstant(1)

        alpha2="_".join(["alpha2",tag])
        if "alpha2" in preconstrains.keys():
            val = preconstrains['alpha2']['val']
            err = 0
        else:
            val = 15.0
            print "attention set value to default in addMjetBackgroundShapeVJets"
        self.w.factory("{name}[{val},{err},{err}]".format(name=alpha2,val=val,err=val))
	self.w.var(alpha2).setConstant(1)

        n2="_".join(["n2",tag])
        if "n2" in preconstrains.keys():
            val = preconstrains['n2']['val']
            err = 0
        else:
            val = 15.0
            print "attention set value to default in addMjetBackgroundShapeVJets"
        self.w.factory("{name}[{val},{err},{err}]".format(name=n2,val=val,err=val))
	self.w.var(n2).setConstant(1)

        pdfName="_".join([name,self.tag])
        if len(resolutionSysts)>=1:
            wjet = ROOT.RooDoubleCB(pdfName,pdfName,self.w.var(Mjet),self.w.function(mean),self.w.function(sigma),self.w.var(alpha),self.w.var(n),self.w.var(alpha2),self.w.var(n2))
        #getattr(self.w,'import')(vvMass,ROOT.RooFit.Rename(pdfName))
        else:
            wjet = ROOT.RooDoubleCB(pdfName,pdfName,self.w.var(Mjet),self.w.var(mean),self.w.var(sigma),self.w.var(alpha),self.w.var(n),self.w.var(alpha2),self.w.var(n2))
        getattr(self.w,'import')(wjet,ROOT.RooFit.Rename(pdfName))
        
        
    

    def addMVVBackgroundShapePow(self,name,variable,newTag="",preconstrains={}):
        
        MVV=variable
        self.w.factory(MVV+"[0,13000]")

        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag

        p0="_".join(["p0",tag])
        if "p0" in preconstrains.keys():
            val = preconstrains['p0']['val']
            err = preconstrains['p0']['err']
            self.addSystematic(p0,"param",[val,err])
        else:
            val = -4
        self.w.factory("{name}[{val},-100,0]".format(name=p0,val=val))

        pdfName="_".join([name,self.tag])
        qcd = ROOT.RooPower(pdfName,pdfName,self.w.var(MVV),self.w.var(p0))
        getattr(self.w,'import')(qcd,ROOT.RooFit.Rename(pdfName))

    def addMjetBackgroundShapeVJetsRes2(self,name,variable,newTag="",preconstrains={},preconstrains2={},scale ={},resolution={}):
        print "start importing shapes for resonant Vjets backgorund "+str(name)
        
        Mjet=variable
        self.w.var(Mjet)
        if self.w.var(Mjet)==None:
            self.w.factory(Mjet+"[0,10000]")
        
        scaleStr='0'
        resolutionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
        for syst,factor in scale.iteritems():
            self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.iteritems():
            self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            resolutionSysts.append(syst)

        if newTag !="": tag=newTag
        else: tag=name+"_"+self.tag

        mean="_".join(["mean",tag])
        if len(scaleSysts)>=1:
            self.w.factory("expr::{name}('{param}*(1+{vv_syst})',{vv_systs})".format(name=mean,param=preconstrains['mean']['val'],vv_syst=scaleStr,vv_systs=','.join(scaleSysts)))
        else:
            self.w.factory("{name}[{val},-{err},{err}]".format(name=mean,val=preconstrains['mean']['val'],err=0))
            self.w.var(mean).setConstant(1)
                
        sigma="_".join(["sigma",tag])
        if len(resolutionSysts)>=1:
            self.w.factory("expr::{name}('({param})*(1+{vv_syst})',{vv_systs})".format(name=sigma,param=preconstrains['sigma']['val'],vv_syst=resolutionStr,vv_systs=','.join(resolutionSysts)))
        else:
            self.w.factory("{name}[{val},-{err},{err}]".format(name=sigma,val=preconstrains['sigma']['val'],err=0))
            self.w.var(sigma).setConstant(1)
        
        alpha="_".join(["alpha",tag])
        self.w.factory("{name}[{val},-{err},{err}]".format(name=alpha,val=preconstrains['alpha']['val'],err=0))
        self.w.var(alpha).setConstant(1)
        
        n="_".join(["n",tag])
        self.w.factory("{name}[{val},-{err},{err}]".format(name=n,val=preconstrains['n']['val'],err=0))
        self.w.var(n).setConstant(1)
        
        alpha2="_".join(["alpha2",tag])
        self.w.factory("{name}[{val},-{err},{err}]".format(name=alpha2,val=preconstrains['alpha2']['val'],err=0))
        self.w.var(alpha2).setConstant(1)
        
        n2="_".join(["n2",tag])
        self.w.factory("{name}[{val},-{err},{err}]".format(name=n2,val=preconstrains['n2']['val'],err=0))
	self.w.var(n2).setConstant(1)
	
	
	mean_g="_".join(["mean_Vjets_nonRes",tag])
        self.w.factory("{name}[{val},-{err},{err}]".format(name=mean_g,val=preconstrains2['c0']['val'],err=0))
	self.w.var(mean_g).setConstant(1)
	
	sigma_g="_".join(["sigma_Vjets_nonRes",tag])
        self.w.factory("{name}[{val},-{err},{err}]".format(name=sigma_g,val=preconstrains2['c1']['val'],err=0))
	self.w.var(sigma_g).setConstant(1)
	
	
	fVar="CMS_VV_JJ_VJets_f"#"_".join(["f",tag])
        self.w.factory("{name}[{val},-{err},{err}]".format(name=fVar,val=preconstrains2['cf']['val'],err=0))
	self.w.var(fVar).setConstant(1)
	print "******************************************************************"
	print fVar

        pdfName="_".join([name,self.tag])
        if len(resolutionSysts)>=1:
            wjet_res = ROOT.RooDoubleCB(pdfName+'_Res',pdfName+'_Res',self.w.var(Mjet),self.w.function(mean),self.w.function(sigma),self.w.var(alpha),self.w.var(n),self.w.var(alpha2),self.w.var(n2))
            gaus = ROOT.RooGaussian(pdfName+'_Vjets_nonRes',pdfName+'_Vjets_nonRes',self.w.var(Mjet),self.w.var(mean_g),self.w.var(sigma_g))
            wjet = ROOT.RooAddPdf(pdfName,pdfName,wjet_res,gaus,self.w.var(fVar))
        else:
            wjet = ROOT.RooDoubleCB(pdfName,pdfName,self.w.var(Mjet),self.w.var(mean),self.w.var(sigma),self.w.var(alpha),self.w.var(n),self.w.var(alpha2),self.w.var(n2))
        getattr(self.w,'import')(wjet,ROOT.RooFit.Rename(pdfName))
    
    def addMVVBackgroundShapeEXPN(self,name,variable,newTag="",preconstrains={}):
        
        MVV=variable
        self.w.factory(MVV+"[0,13000]")

        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag

        p0="_".join(["p0",tag])
        if "p0" in preconstrains.keys():
            val = preconstrains['p0']['val']
            err = preconstrains['p0']['err']
            self.addSystematic(p0,"param",[val,err])
        else:
            val = -0.001
        self.w.factory("{name}[{val},-5,0]".format(name=p0,val=val))

        p1="_".join(["p1",tag])
        if "p0" in preconstrains.keys():
            val = preconstrains['p1']['val']
            err = preconstrains['p1']['err']
            self.addSystematic(p1,"param",[val,err])
        else:
            val = 0
        self.w.factory("{name}[{val},-10,10]".format(name=p1,val=val))


        pdfName="_".join([name,self.tag])
        qcd = ROOT.RooExpNPdf(pdfName,pdfName,self.w.var(MVV),self.w.var(p0),self.w.var(p1))
        getattr(self.w,'import')(qcd,ROOT.RooFit.Rename(pdfName))


    def addMVVBackgroundShapeEXPNFromMC(self,name,variable,newTag,filename):
        f=ROOT.TFile(filename)
        hist=f.Get(name)
        func = ROOT.TF1("func","[0]*exp([1]*x+[2]/x)",500,10000)
        hist.Fit(func)
        MVV=variable
        self.w.factory(MVV+"[0,13000]")

        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag

        p0="_".join(["p0",tag])
        val = -0.001
        self.w.factory("{name}[{val}]".format(name=p0,val=func.GetParameter(1)))

        p1="_".join(["p1",tag])
        val = 0
        self.w.factory("{name}[{val}]".format(name=p1,val=func.GetParameter(2)))


        pdfName="_".join([name,self.tag])
        qcd = ROOT.RooExpNPdf(pdfName,pdfName,self.w.var(MVV),self.w.var(p0),self.w.var(p1))
        getattr(self.w,'import')(qcd,ROOT.RooFit.Rename(pdfName))


    def addMVVBackgroundShapeErfPow(self,name,variable,newTag="",preconstrains={}):
        
        MVV=variable
        self.w.factory(MVV+"[0,13000]")

        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag

        p0="_".join(["p0",tag])
        if "p0" in preconstrains.keys():
            val = preconstrains['p0']['val']
            err = preconstrains['p0']['err']
            self.addSystematic(p0,"param",[val,err])
        else:
            val = -0.1
        self.w.factory("{name}[{val},-20,0]".format(name=p0,val=val))


        p1="_".join(["p1",tag])
        if "p1" in preconstrains.keys():
            val = preconstrains['p1']['val']
            err = preconstrains['p1']['err']
            self.addSystematic(p1,"param",[val,err])
        else:
            val = 700
        self.w.factory("{name}[{val},0,2000]".format(name=p1,val=val))


        p2="_".join(["p2",tag])
        if "p2" in preconstrains.keys():
            val = preconstrains['p2']['val']
            err = preconstrains['p2']['err']
            self.addSystematic(p2,"param",[val,err])
        else:
            val = 1000
        self.w.factory("{name}[{val},0,5000]".format(name=p2,val=val))


        pdfName="_".join([name,self.tag])
        qcd = ROOT.RooErfPowPdf(pdfName,pdfName,self.w.var(MVV),self.w.function(p0),self.w.function(p1),self.w.function(p2))

        getattr(self.w,'import')(qcd,ROOT.RooFit.Rename(pdfName))



    def addParametricMVVBKGShapeErfPow(self,name,MVV,MJJ,jsonFile,newTag="",systs0={},systs1={},systs2={},pdfTag=""):

        syst0Str='0'
        syst1Str='0'
        syst2Str='0'


        systsV0=[]
        systsV1=[]
        systsV2=[]

        for syst,factor in systs0.iteritems():
            self.w.factory(syst+"[0,-0.1,0.1]")
            syst0Str+="+{factor}*{syst}".format(factor=factor,syst=syst)
            systsV0.append(syst)

        for syst,factor in systs1.iteritems():
            self.w.factory(syst+"[0,-0.1,0.1]")
            syst1Str+="+{factor}*{syst}".format(factor=factor,syst=syst)
            systsV1.append(syst)

        for syst,factor in systs2.iteritems():
            self.w.factory(syst+"[0,-0.1,0.1]")
            syst2Str+="+{factor}*{syst}".format(factor=factor,syst=syst)
            systsV2.append(syst)


       
        self.w.factory(MVV+"[0,13000]")
        self.w.factory(MJJ+"[0,1000]")

        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag


        f=open(jsonFile)
        info=json.load(f)

        p0="_".join(["p0",name,tag])
        if syst0Str.find(MJJ)!=-1 or str(info['p0']).find('mjj')!=-1:
            self.w.factory("expr::{name}('({param})*(1+0*{MJJ}+{syst})',{MJJ},{systs})".format(name=p0,param=str(info['p0']).replace("mjj",MJJ),MJJ=MJJ,syst=syst0Str,systs=','.join(systsV0)))
        else:
            self.w.factory("expr::{name}('({param})*(1+{syst})',{systs})".format(name=p0,param=str(info['p0']).replace("mjj",MJJ),syst=syst0Str,systs=','.join(systsV0)))


        p1="_".join(["p1",name,tag])
        if syst1Str.find(MJJ)!=-1 or str(info['p1']).find('mjj')!=-1:
            self.w.factory("expr::{name}('({param})*(1+0*{MJJ}+{syst})',{MJJ},{systs})".format(name=p1,param=str(info['p1']).replace("mjj",MJJ),MJJ=MJJ,syst=syst1Str,systs=','.join(systsV1)))
        else:
            self.w.factory("expr::{name}('({param})*(1+{syst})',{systs})".format(name=p1,param=str(info['p1']).replace("mjj",MJJ),syst=syst1Str,systs=','.join(systsV1)))



        p2="_".join(["p2",name,tag])
        if syst2Str.find(MJJ)!=-1 or str(info['p2']).find('mjj')!=-1:
            self.w.factory("expr::{name}('({param})*(1+0*{MJJ}+{syst})',{MJJ},{systs})".format(name=p2,param=str(info['p2']).replace("mjj",MJJ),MJJ=MJJ,syst=syst2Str,systs=','.join(systsV2)))
        else:
            self.w.factory("expr::{name}('({param})*(1+{syst})',{systs})".format(name=p2,param=str(info['p2']).replace("mjj",MJJ),syst=syst2Str,systs=','.join(systsV2)))
        
        if pdfTag=="":
            pdfName="_".join([name,self.tag])
        else:
            pdfName="_".join([name,pdfTag])

        erfexp = ROOT.RooErfPowPdf(pdfName,pdfName,self.w.var(MVV),self.w.function(p0),self.w.function(p1),self.w.function(p2))
        getattr(self.w,'import')(erfexp,ROOT.RooFit.Rename(name))



    def addParametricMVVBKGShapePow(self,name,MVV,MJJ,jsonFile,newTag="",systs0={}):
        syst0Str='0'
        systsV0=[]
        for syst,factor in systs0.iteritems():
            self.w.factory(syst+"[0,-0.1,0.1]")
            syst0Str+="+{factor}*{syst}".format(factor=factor,syst=syst)
            systsV0.append(syst)
       
        self.w.factory(MVV+"[0,13000]")
        self.w.factory(MJJ+"[0,1000]")

        if newTag !="":
            tag=newTag
        else:
            tag=name+"_"+self.tag


        f=open(jsonFile)
        info=json.load(f)



        p0="_".join(["p0",name,tag])
        self.w.factory("expr::{name}('({param})*(1+0*{MJJ}+{syst})',{MJJ},{systs})".format(name=p0,param=info['p0'].replace("mjj",MJJ),MJJ=MJJ,syst=syst0Str,systs=','.join(systsV0)))

        pdfName="_".join([name,self.tag])
        qcd = ROOT.RooPower(pdfName,pdfName,self.w.var(MVV),self.w.function(p0))
        getattr(self.w,'import')(erfexp,ROOT.RooFit.Rename(name))




    def sum(self,name,pdf1,pdf2,sumVar,sumVarExpr=""):
        pdfName="_".join([name,self.tag])
        pdfName1="_".join([pdf1,self.tag])
        pdfName2="_".join([pdf2,self.tag])
        #if sumVarExpr=='':
        #    self.w.factory(sumVar+"[0,1]")
        #else:    
        #    self.w.factory("expr::"+sumVar+"("+sumVarExpr+")")
        self.w.factory("SUM::{name}({f}*{name1},{name2})".format(name=pdfName,name1=pdfName1,f=sumVar,name2=pdfName2))
        
    def sumPdf(self,name,pdf1,pdf2,variable):
        #self.w.factory(variable+"[0,0.,1.]")
        print "variable "+str(variable)
        self.w.factory(variable+"[0.5,0.5,0.5]")
        self.w.var(variable).setConstant(1)
        pdfName="_".join([name,self.tag])
        pdfName1="_".join([pdf1,self.tag])
        pdfName2="_".join([pdf2,self.tag])
        
           
        self.w.factory("SUM::{name}(expr::{name2}_ratio1('(1-{f})',{f})*{name1},{f}*{name2})".format(name=pdfName,name1=pdfName1,f=variable,name2=pdfName2))
       
            
      

    def conditionalProduct(self,name,pdf1,varName,pdf2,pdf3,tag1="",tag2="",tag3=""):
        pdfName="_".join([name,self.tag])

        if tag1=="":
            pdfName1="_".join([pdf1,self.tag])
        else:
            pdfName1="_".join([pdf1,tag1])
        if tag2=="":    
            pdfName2="_".join([pdf2,self.tag])
        else:
            pdfName2="_".join([pdf2,tag2])
        if tag3=="":    
            pdfName3="_".join([pdf3,self.tag])
        else:
            pdfName3="_".join([pdf3,tag3])
	    
        self.w.factory("PROD::{name}({name1}|{x},{name2}|{x},{name3})".format(name=pdfName,name1=pdfName1,x=varName,name2=pdfName2,name3=pdfName3))
        
    def conditionalProduct2(self,name,pdf1,pdf2,pdf3,varName,tag1="",tag2="",tag3=""):
        pdfName="_".join([name,self.tag])

        if tag1=="":
            pdfName1="_".join([pdf1,self.tag])
        else:
            pdfName1="_".join([pdf1,tag1])
        if tag2=="":    
            pdfName2="_".join([pdf2,self.tag])
        else:
            pdfName2="_".join([pdf2,tag2])
        if tag3=="":    
            pdfName3="_".join([pdf3,self.tag])
        else:
            pdfName3="_".join([pdf3,tag3])
	    
        self.w.factory("PROD::{name}({name1},{name2},{name3}| {x})".format(name=pdfName,name1=pdfName1,x=varName,name2=pdfName2,name3=pdfName3))    

    def product(self,name,pdf1,pdf2):
        pdfName="_".join([name,self.tag])
        pdfName1="_".join([pdf1,self.tag])
        pdfName2="_".join([pdf2,self.tag])
        self.w.factory("PROD::{name}({name1},{name2})".format(name=pdfName,name1=pdfName1,name2=pdfName2))

    def product3D(self,name,pdf1,pdf2,pdf3):
        print self.tag
        pdfName="_".join([name,self.tag])
        pdfName1="_".join([pdf1,self.tag])
        pdfName2="_".join([pdf2,self.tag])
        pdfName3="_".join([pdf3,self.tag])
        print pdfName
        print pdfName1
        print pdfName2
        print pdfName3
        self.w.factory("PROD::{name}({name1},{name2},{name3})".format(name=pdfName,name1=pdfName1,name2=pdfName2,name3=pdfName3))
    
    
    def envelope(self,name,pdfs):
        catName = "envelope_"+name+"_"+self.tag
        pdfName="_".join([name,self.tag])
        
        pdfList=[]
        pdfArgList = ROOT.RooArgList()
        for p in pdfs:
            pdfList.append(p+"_"+self.tag)
            pdfArgList.add(self.w.pdf(p+"_"+self.tag))

        pdfStr = ','.join(pdfList)
        self.w.factory("{cat}[{list}]".format(cat=catName,list=pdfStr))
        multiPDF = ROOT.RooMultiPdf(pdfName,pdfName,self.w.cat(catName),pdfArgList)
        getattr(self.w,'import')(multiPDF)
        self.addSystematic(catName,"discrete","")

    def conditionalDoubleProduct(self,name,pdf1,pdf2,varName,pdf3):
        pdfName="_".join([name,self.tag])
        pdfName1="_".join([pdf1,self.tag])
        pdfName2="_".join([pdf2,self.tag])
        pdfName3="_".join([pdf3,self.tag])
        self.w.factory("PROD::{name}({name1}|{x},{name2}|{x},{name3})".format(name=pdfName,name1=pdfName1,x=varName,name2=pdfName2,name3=pdfName3))


    def product(self,name,pdf1,pdf2):
        pdfName="_".join([name,self.tag])
        pdfName1="_".join([pdf1,self.tag])
        pdfName2="_".join([pdf2,self.tag])
        self.w.factory("PROD::{name}({name1},{name2})".format(name=pdfName,name1=pdfName1,name2=pdfName2))

    def addParametricYield(self,name,ID,jsonFile,constant=1.0):
        f=open(jsonFile)
        info=json.load(f)

        pdfName="_".join([name,self.tag])
        pdfNorm="_".join([name,self.tag,"norm"])
        self.w.factory("expr::{name}('({param})*{lumi}*{constant}',MH,{lumi})".format(name=pdfNorm,param=info['yield'],lumi=self.physics+"_"+self.period+"_lumi",constant=constant))       
        f.close()
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':1.0})



    def addParametricYieldWithUncertainty(self,name,ID,jsonFile,constant,uncertaintyName,uncertaintyFormula,uncertaintyValue):
        f=open(jsonFile)
        info=json.load(f)
        pdfName="_".join([name,self.tag])
        pdfNorm="_".join([name,self.tag,"norm"])
        self.w.factory(uncertaintyName+'[0,-1,1]')
        #self.w.factory("expr::{name}('({param})*{lumi}*({constant}+{unc}*{form})',MH,{lumi},{unc})".format(name=pdfNorm,param=info['yield'],lumi=self.physics+"_"+self.period+"_lumi",constant=constant,unc=uncertaintyName,form=uncertaintyFormula))
        #this formula is better: y*exp(log(unc)*theta)
        self.w.factory("expr::{name}('({param})*{lumi}*(exp(log({form})*{unc}))',MH,{lumi},{unc})".format(name=pdfNorm,param=info['yield'],lumi=self.physics+"_"+self.period+"_lumi",constant=constant,unc=uncertaintyName,form=uncertaintyFormula))
        self.addSystematic(uncertaintyName,"param",[0,uncertaintyValue])
        f.close()
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':1.0})


    def addParametricYieldWithCrossSection(self,name,ID,jsonFile,jsonFileCS,sigmaStr,BRStr):
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        from array import array
        #first load cross section
        fCS=open(jsonFileCS)
        info=json.load(fCS)
        xArr=[]
        yArr=[]

        for m in sorted(map(float,info.keys())):
            xArr.append(float(m))
            #I know this is stupid 
            yArr.append(float(info[str(int(m))][sigmaStr])*float(info[str(int(m))][BRStr]))


        pdfSigma="_".join([name,self.tag,"sigma"])
        spline=ROOT.RooSpline1D(pdfSigma,pdfSigma,self.w.var("MH"),len(xArr),array('d',xArr),array('d',yArr))    
        getattr(self.w,'import')(spline,ROOT.RooFit.Rename(pdfSigma))
        fCS.close()

        f=open(jsonFile)
        info=json.load(f)      
        pdfName="_".join([name,self.tag])
        pdfNorm="_".join([name,self.tag,"norm"])
        self.w.factory("expr::{name}('({param})*{lumi}*({sigma})',MH,{lumi},{sigma})".format(name=pdfNorm,param=info['yield'],lumi=self.physics+"_"+self.period+"_lumi",sigma=pdfSigma))       
        f.close()
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':1.0})

    def addParametricYieldWithCrossSectionAndUncertainties(self,name,ID,jsonFile,jsonFileCS,sigmaStr,BRStr,sigmaStrP,sigmaStrM,constant,uncertaintyName,uncertaintyFormula,uncertaintyValue):
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        from array import array
        #first load cross section
        fCS=open(jsonFileCS)
        info=json.load(fCS)
        xArr=[]
        yArr=[]
        yErrArr=[]

        for m in sorted(map(float,info.keys())):
            xArr.append(float(m))
            #I know this is stupid 
            yArr.append(float(info[str(int(m))][sigmaStr])*float(info[str(int(m))][BRStr]))
            yErrArr.append(0.5*(float(info[str(int(m))][sigmaStrP])-float(info[str(int(m))][sigmaStrM]))/float(info[str(int(m))][sigmaStr]))

        pdfSigma="_".join([name,self.tag,"sigma"])
        spline=ROOT.RooSpline1D(pdfSigma,pdfSigma,self.w.var("MH"),len(xArr),array('d',xArr),array('d',yArr))    
        getattr(self.w,'import')(spline,ROOT.RooFit.Rename(pdfSigma))

        pdfErrSigma="_".join([name,self.tag,"sigmaErr"])
        splineErr=ROOT.RooSpline1D(pdfErrSigma,pdfErrSigma,self.w.var("MH"),len(xArr),array('d',xArr),array('d',yErrArr))    
        getattr(self.w,'import')(splineErr,ROOT.RooFit.Rename(pdfErrSigma))
        fCS.close()

        f=open(jsonFile)
        info=json.load(f)      
        pdfName="_".join([name,self.tag])
        pdfNorm="_".join([name,self.tag,"norm"])

        self.w.factory("sigmaUnc[-3,3]")
        self.w.factory(uncertaintyName+'[0,-1,1]')

        self.w.factory("expr::{name}('({param})*{lumi}*({sigma})*(1+sigmaUnc*{sigmaErr})*({constant}+{unc}*{form})',MH,{lumi},{sigma},{sigmaErr},sigmaUnc,{unc})".format(name=pdfNorm,param=info['yield'],lumi=self.physics+"_"+self.period+"_lumi",sigma=pdfSigma,sigmaErr=pdfErrSigma,constant=constant,unc=uncertaintyName,form=uncertaintyFormula))       
        f.close()
        self.addSystematic('sigmaUnc',"param",[0,1])
        self.addSystematic(uncertaintyName,"param",[0,uncertaintyValue])
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':1.0})



    def addParametricYieldHVTBR(self,name,ID,jsonFile,jsonFileCS,sigmaStr,BRStr,constant,uncertaintyName,uncertaintyFormula,uncertaintyValue):
        print 'I will only assume the BRs from HVT and float the cross section'
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        from array import array
        #first load cross section
        fCS=open(jsonFileCS)
        info=json.load(fCS)
        xArr=[]
        yArr=[]
        yErrArr=[]

        sigmaStr = sigmaStr.split(",")
        for m in sorted(map(float,info.keys())):
            xArr.append(float(m))
            #I know this is stupid 
	    value = 0
	    for i,s in enumerate(sigmaStr):
	     value+=float(info[str(int(m))][sigmaStr[i]])*float(info[str(int(m))][BRStr])
            yArr.append(value)
	    print m,value

        pdfSigma="_".join([name,self.tag,"sigma"])
        spline=ROOT.RooSpline1D(pdfSigma,pdfSigma,self.w.var("MH"),len(xArr),array('d',xArr),array('d',yArr))    
        getattr(self.w,'import')(spline,ROOT.RooFit.Rename(pdfSigma))

        f=open(jsonFile)
        info=json.load(f)      
        pdfName="_".join([name,self.tag])
        pdfNorm="_".join([name,self.tag,"norm"])

        self.w.factory(uncertaintyName+'[0,-1,1]')

        self.w.factory("expr::{name}('({param})*{lumi}*({sigma})*({constant}+{unc}*{form})',MH,{lumi},{sigma},{unc})".format(name=pdfNorm,param=info['yield'],lumi=self.physics+"_"+self.period+"_lumi",sigma=pdfSigma,constant=constant,unc=uncertaintyName,form=uncertaintyFormula))       
        f.close()
        self.addSystematic(uncertaintyName,"param",[0,uncertaintyValue])
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':1.0})



    def addFloatingYield(self,name,ID,events,mini=0,maxi=1e+9,constant=False):
        pdfName="_".join([name,self.tag])
        pdfNorm="_".join([name,self.tag,"norm"])
        self.w.factory("{name}[{val},{mini},{maxi}]".format(name=pdfNorm,val=events,mini=mini,maxi=maxi))       
        if constant:
            self.w.var(pdfNorm).setConstant(1)
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':1.0})


    def addConstrainedYield(self,name,ID,events,nuisance,uncertainty):
        pdfName="_".join([name,self.tag])
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':events})
        self.addSystematic(nuisance,"lnN",{name:1+uncertainty})

    def addConstrainedYieldFromFile(self,name,ID,filename,histoName,nuisance,uncertainty):
        pdfName="_".join([name,self.tag])

        f=ROOT.TFile(filename)
        histogram=f.Get(histoName)
        events=histogram.Integral()*self.luminosity
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':events})
        self.addSystematic(nuisance,"lnN",{name:1+uncertainty})

    def addFixedYieldFromFile(self,name,ID,filename,histoName,constant=1.0):
        pdfName="_".join([name,self.tag])
        f=ROOT.TFile(filename)
        histogram=f.Get(histoName)
        events=histogram.Integral()*self.luminosity*constant
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':events})

    def addYieldWithRateParameter(self,name,ID,paramName,formula,values):#jen        
        pdfName="_".join([name,self.tag])
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':1.0})            
        kind = 'rateParam\t%s\t%s\t%s'%(self.tag,name,formula)
        self.systematics.append({'name':paramName,'kind':kind,'values':values })
        
    def addYieldWithRateParameterFromFile(self,name,ID,paramName,filename,histoName,values=[],scaleFactor=1):#jen        
        pdfName="_".join([name,self.tag])
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':1.0})    

        f=ROOT.TFile(filename)
        histogram=f.Get(histoName)
        events=histogram.Integral()*self.luminosity*scaleFactor
        
        kind = 'rateParam\t%s\t%s\t%f'%(self.tag,name,events)
        self.systematics.append({'name':paramName,'kind':kind,'values':values})

    def addYieldWithRateParameterFromFileWithUncertainty(self,name,ID,paramName,filename,histoName,uncertaintyName,uncertaintyFormula,uncertaintyValue,values=[],scaleFactor=1):#jen        
        pdfName="_".join([name,self.tag])
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':1.0})    

        f=ROOT.TFile(filename)
        histogram=f.Get(histoName)
        events=histogram.Integral()*self.luminosity*scaleFactor
        
        kind = 'rateParam\t%s\t%s\t%f'%(self.tag,name,events)
        self.systematics.append({'name':paramName,'kind':kind,'values':values})

        pdfNorm="_".join([name,self.tag,"norm"])
        if not self.w.var(uncertaintyName):
         self.w.factory(uncertaintyName+'[0,-1,1]')
         self.addSystematic(uncertaintyName,"param",[0,uncertaintyValue])
        self.w.factory("expr::{name}('exp(log({form})*{unc})',MH,{unc})".format(name=pdfNorm,unc=uncertaintyName,form=uncertaintyFormula))

    def addYieldWithRateParameterWithUncertainty(self,name,ID,paramName,formula,values,uncertaintyName,uncertaintyFormula,uncertaintyValue):#jen        
        pdfName="_".join([name,self.tag])
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':1.0})            
        kind = 'rateParam\t%s\t%s\t%s'%(self.tag,name,formula)
        self.systematics.append({'name':paramName,'kind':kind,'values':values})
        
        pdfNorm="_".join([name,self.tag,"norm"])
        if not self.w.var(uncertaintyName):
         self.w.factory(uncertaintyName+'[0,-1,1]')
         self.addSystematic(uncertaintyName,"param",[0,uncertaintyValue])
        self.w.factory("expr::{name}('exp(log({form})*{unc})',MH,{unc})".format(name=pdfNorm,unc=uncertaintyName,form=uncertaintyFormula))

       
    def makeCard(self):

        if self.cat!="":f = open("datacard_"+self.cat+'.txt','w')
        else: f = open("datacard_"+self.tag+'.txt','w')
        f.write('imax 1\n')
        f.write('jmax {n}\n'.format(n=len(self.contributions)-1))
        f.write('kmax *\n')
        f.write('-------------------------\n')
        for c in self.contributions:
            if self.cat!="": f.write('shapes {name} {channel} {file}.root w:{pdf}\n'.format(name=c['name'],channel=self.tag,file="datacardInputs_"+self.cat,pdf=c['pdf']))
            else: f.write('shapes {name} {channel} {file}.root w:{pdf}\n'.format(name=c['name'],channel=self.tag,file="datacardInputs_"+self.tag,pdf=c['pdf']))
        if self.cat!="": f.write('shapes {name} {channel} {file}.root w:{name}\n'.format(name="data_obs",channel=self.tag,file="datacardInputs_"+self.cat))  
        else: f.write('shapes {name} {channel} {file}.root w:{name}\n'.format(name="data_obs",channel=self.tag,file="datacardInputs_"+self.tag))
        f.write('-------------------------\n')
        f.write('bin '+self.tag+'\n')
        f.write('observation  -1\n')
        f.write('-------------------------\n')
        f.write('bin\t') 

        for shape in self.contributions:
            f.write(self.tag+'\t')
        f.write('\n')

        #Sort the shapes by ID 
 
        shapes = sorted(self.contributions,key=lambda x: x['ID'])
        #print names
        f.write('process\t')
        for shape in shapes:
            f.write(shape['name']+'\t')
        f.write('\n')

        #Print ID
        f.write('process\t')
        for shape in shapes:
            f.write(str(shape['ID'])+'\t')
        f.write('\n')

        #print rates
        f.write('rate\t')
        for shape in shapes:
            f.write(str(shape['yield'])+'\t')
        f.write('\n')


        #Now systematics
        for syst in self.systematics:
            if syst['kind'] == 'param':
                f.write(syst['name']+'\t'+'param\t' +str(syst['values'][0])+'\t'+str(syst['values'][1])+'\n')

            elif 'rateParam' in syst['kind']:
                line = syst['name']+'\t'+str(syst['kind'])+"\t"
                for v in syst['values']:
                 line+=str(v)+","
                line=line[:-1]
                line+='\n' 
                f.write(line)
                
            elif syst['kind'] == 'discrete':
                f.write(syst['name']+'\t'+'discrete\n')

            elif syst['kind'] == 'lnN' :
                f.write(syst['name']+'\t'+ 'lnN\t' )
                for shape in shapes:
                    has=False
                    for name,v in syst['values'].iteritems():
                        if shape['name']==name:
                            f.write(str(v)+'\t' )
                            has=True
                            break;
                    if not has:
                            f.write('-\t' )
                f.write('\n' )
            elif syst['kind'] == 'lnU': 
                f.write(syst['name']+'\t'+ 'lnU\t' )
                for shape in shapes:
                    has=False
                    for name,v in syst['values'].iteritems():
                        if shape['name']==name:
                            f.write(str(v)+'\t' )
                            has=True
                            break;
                    if not has:
                            f.write('-\t' )
                f.write('\n' )
                            
                        
        f.close()


        self.rootFile.cd()
        self.w.Write("w",0,2000000000)
        self.rootFile.Close()
            
    
        

    def importBinnedData(self,filename,histoname,poi,name = "data_obs",scale=1):
        f=ROOT.TFile(filename)
        histogram=f.Get(histoname)
        #histogram.Scale(scale*self.luminosity)	#jen
        histogram.Scale(scale)	#jen
        cList = ROOT.RooArgList()
        for i,p in enumerate(poi):
            cList.add(self.w.var(p))
            if i==0:
                axis=histogram.GetXaxis()
            elif i==1:
                axis=histogram.GetYaxis()
            elif i==2:
                axis=histogram.GetZaxis()
            else:
                print 'Asking for more than 3 D . ROOT doesnt support that, use unbinned data instead'
                return
            mini=axis.GetXmin()
            maxi=axis.GetXmax()
            bins=axis.GetNbins()
            self.w.var(p).setMin(mini)
            self.w.var(p).setMax(maxi)
            self.w.var(p).setBins(bins)
        dataHist=ROOT.RooDataHist(name,name,cList,histogram)

        getattr(self.w,'import')(dataHist,ROOT.RooFit.Rename(name))
        
