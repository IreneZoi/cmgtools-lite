# VV statistical analysis in 10X

Prepare the working directory with Higgs Combine Tools. Use the 10X release compatible with the [UHH framework](https://github.com/UHH2/UHH2). If you have that already installed you do
not need to check out the CMSSW release again. Be carefull this may impact the crab job submission because of a too big sandbox. Also, installation on lxplus is suggested as it will work much faster without many modifications.

```
mkdir VVAnalysisWith2DFit
mkdir CMGToolsForStat10X
cd CMGToolsForStat10X
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_2_10
cd CMSSW_10_2_10/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.0.0
scramv1 b clean && scramv1 b -j 8
```

Fork cmgtools from https://github.com/Diboson3D/cmgtools-lite and checkout the VV statistical tools

```
cd ../..
export GITUSER=`git config user.github`
git clone https://github.com/${GITUSER}/cmgtools-lite CMGTools
cd CMGTools
git remote add Diboson3D https://github.com/Diboson3D/cmgtools-lite -b VV_VH
git fetch Diboson3D
git checkout -b VV_VH Diboson3D/VV_VH
scram b -j 8
cd VVResonances/interactive
ln -s samples_location sample
```

Current sample location with random sorting of jet1 and jet2

```
/eos/cms/store/cmst3/group/exovv/VVtuple/FullRun2VVVHNtuple/2016/
/eos/cms/store/cmst3/group/exovv/VVtuple/FullRun2VVVHNtuple/2018/
```

Make the 3D templates. Several options can be specified to produce QCD, V+Jets or signal templates, normalization etc.
 
```
python makeInputs.py
```

Run closure test of signal fits:

```
python plotSignalShapesFromJSON.py -f JJ_BulkGravWW_2016_MJl1_VV_HPLP.json -v mJ
python plotSignalShapesFromJSON.py -f JJ_BulkGravWW_2016_MJl2_VV_HPLP.json -v mJ -l "l2"
python plotSignalShapesFromJSON.py -f JJ_BulkGravWW_2016_MVV.json -v mVV
```

After producing the QCD background templetes, to improve the agreement between MC and template, it is necessary to fit the  HPLP MC with HPLP kernel ( with the option -p it is possible to select different projections: x for mjet1, y for mjet2 and z for mjj - use just one option at the time to avoid crashes! - e.g. -p z -x 65,105 -y 65,105 gives mjj projection in the mjet1&2 range 65,105) 

The script expects to find the files in a directory called results_year.

```
mkdir postfit_qcd #if you do not have it already
python transferKernel.py -i results_2016/JJ_2016_nonRes_VV_HPLP.root --sample pythia --year 2016 -p z --pdfIn results_2016/JJ_2016_nonRes_3D_VV_HPLP.root
python transferKernel.py -i results_2016/JJ_2016_nonRes_VV_HPLP.root --sample herwig --year 2016 -p z --pdfIn results_2016/JJ_2016_nonRes_3D_VV_HPLP.root 
python transferKernel.py -i results_2016/JJ_2016_nonRes_VV_HPLP.root --sample madgraph --year 2016 -p z --pdfIn results_2016/JJ_2016_nonRes_3D_VV_HPLP.root
```

merge 1Dx2Dx2D HPLP kernels

```
python transferKernel.py -i results_2016/JJ_2016_nonRes_VV_HPLP.root --sample pythia --year 2016 -p z --pdfIn results_2016/JJ_2016_nonRes_3D_VV_HPLP.root --merge
```

Results validation:

```
python Projections3DHisto.py --mc results_2016/JJ_2016_nonRes_VV_HPLP.root,nonRes -k save_new_shapes_2016_pythia_VV_HPLP_3D.root,histo -o control-plots-HPLP-pythia
python Projections3DHisto.py --mc results_2016/JJ_2016_nonRes_VV_HPLP_altshapeUp.root,nonRes -k save_new_shapes_2016_pythia_VV_HPLP_3D.root,histo -o control-plots-HPLP-herwig
python Projections3DHisto.py --mc results_2016/JJ_2016_nonRes_VV_HPLP_altshape2.root,nonRes -k save_new_shapes_2016_pythia_VV_HPLP_3D.root,histo -o control-plots-HPLP-madgraph
```

fit HPHP MC with post-fit HPLP kernel

```
python transferKernel.py -i results_2016/JJ_2016_nonRes_VV_HPHP.root --sample pythia --year 2016 -p z --pdfIn results_2016/JJ_2016_nonRes_3D_VV_HPLP.root 
python transferKernel.py -i results_2016/JJ_2016_nonRes_VV_HPHP.root --sample herwig --year 2016 -p z --pdfIn results_2016/JJ_2016_nonRes_3D_VV_HPLP.root 
python transferKernel.py -i results_2016/JJ_2016_nonRes_VV_HPHP.root --sample madgraph --year 2016 -p z --pdfIn results_2016/JJ_2016_nonRes_3D_VV_HPLP.root
```

merge 1Dx2Dx2D HPHP kernels

```
python transferKernel.py -i results_2016/JJ_2016_nonRes_VV_HPHP.root --sample madgraph --year 2016 -p z --pdfIn results_2016/JJ_2016_nonRes_3D_VV_HPLP.root --merge
```

Results validation:

```
python Projections3DHisto_HPHP.py --mc results_2016/JJ_2016_nonRes_VV_HPHP.root,nonRes -k save_new_shapes_2016_pythia_VV_HPHP_3D.root,histo -o control-plots-HPHP-pythia
python Projections3DHisto_HPHP.py --mc results_2016/JJ_2016_nonRes_VV_HPHP_altshapeUp.root,nonRes -k save_new_shapes_2016_pythia_VV_HPHP_3D.root,histo -o control-plots-HPHP-herwig
python Projections3DHisto_HPHP.py --mc results_2016/JJ_2016_nonRes_VV_HPHP_altshape2.root,nonRes -k save_new_shapes_2016_pythia_VV_HPHP_3D.root,histo -o control-plots-HPHP-madgraph
```


Create datacard and workspace and run post fit (Some major modification are needed to include the new categories and handle the different years + many background are still missing)

```
python makeCard.py
text2workspace.py datacard_JJ_HPHP_13TeV.txt -o JJ_BulkGWW_HPHP_13TeV_workspace.root
python runPostFit.py
```

Run the limits with combine and make final plot (but not tested because datacards missing, will need modification for desy condor)

```
vvSubmitLimits.py JJ_BulkGWW_HPHP_13TeV_workspace.root -s 100 -q "workday" -m 1200 -M 4200 -C 1
find higgsCombineTest.Asymptotic.* -size +1500c | xargs hadd Limits_BulkGWW_HPHP_13TeV.root
vvMakeLimitPlot.py Limits_BulkGWW_HPHP_13TeV.root -x 1200 -X 4200 #(expected limits)
vvMakeLimitPlot.py Limits_BulkGWW_HPHP_13TeV.root -x 1200 -X 4200 -b 0 #(expected+observed limits)
```

Further instructions on how to run the code can be found at:
https://docs.google.com/document/d/1hU84u27mY85UaAK5R11OHYctBMckDU6kX7IcorboZf8/edit?usp=sharing

