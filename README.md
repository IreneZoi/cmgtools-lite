# VV statistical analysis in 10X

Prepare the working directory with Higgs Combine Tools. Use the 10X release compatible with the [UHH framework](https://github.com/UHH2/UHH2). If you have that already installed you do
not need to check out the CMSSW release again.

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
/eos/cms/store/cmst3/group/exovv/VVtuple/FullRun2VVVHNtuple/2016_new/
/eos/cms/store/cmst3/group/exovv/VVtuple/FullRun2VVVHNtuple/2018/
```

Make the 3D templates
 
```
python makeInputs.py
```

Run closure test of signal fits:

```
python plotSignalShapesFromJSON.py -f JJ_BulkGravWW_2016_MJl1_VV_HPLP.json -v mJ
python plotSignalShapesFromJSON.py -f JJ_BulkGravWW_2016_MJl2_VV_HPLP.json -v mJ -l "l2"
python plotSignalShapesFromJSON.py -f JJ_BulkGravWW_2016_MVV.json -v mVV
```

Run closure test of 3D templates versus simulation with Projections3DHisto.py script:

```
python Projections3DHisto.py --mc JJ_nonRes_VH_HPLP.root,nonRes -k JJ_2016_nonRes_3D_VH_HPLP.root,histo -o control-plots-HPLP-pythia

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

