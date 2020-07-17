#!bin/bash
indir=$1
name=$2_OPT3_0and2_OPTXY0and2_OPTZ1and2_PTZ2and2 #OPT31and1_OPTXY2_OPTZ2
postfitdir=postfit_qcd/${name}/
mkdir $postfitdir
categories=('VH_LPHP' 'VH_HPHP' 'VH_HPLP' 'VV_HPLP' 'VV_HPHP')


for cat in ${categories[*]}; do
    #               pythia
    python transferKernelUnique.py -i ${indir}JJ_2016_nonRes_${cat}.root --sample pythia --year 2016 -p x --pdfIn ${indir}JJ_2016_nonRes_3D_NP.root | tee KernelTransf_2016_${cat}_pythia_x_${name}.out
    #python transferKernel.py -i ${indir}JJ_2016_nonRes_${cat}.root --sample pythia --year 2016 -p y --pdfIn ${indir}JJ_2016_nonRes_3D_NP.root | tee KernelTransf_2016_${cat}_pythia_y_${name}.out
    #python transferKernel.py -i ${indir}JJ_2016_nonRes_${cat}.root --sample pythia --year 2016 -p z --pdfIn ${indir}JJ_2016_nonRes_3D_NP.root | tee KernelTransf_2016_${cat}_pythia_z_${name}.out

    #                madgraph
    python transferKernelUnique.py -i ${indir}JJ_2016_nonRes_${cat}.root --sample madgraph --year 2016 -p x --pdfIn ${indir}JJ_2016_nonRes_3D_NP.root | tee KernelTransf_2016_${cat}_madgraph_x_${name}.out                   
    #python transferKernel.py -i ${indir}JJ_2016_nonRes_${cat}.root --sample madgraph --year 2016 -p y --pdfIn ${indir}JJ_2016_nonRes_3D_NP.root | tee KernelTransf_2016_${cat}_madgraph_y_${name}.out
    #python transferKernel.py -i ${indir}JJ_2016_nonRes_${cat}.root --sample madgraph --year 2016 -p z --pdfIn ${indir}JJ_2016_nonRes_3D_NP.root | tee KernelTransf_2016_${cat}_madgraph_z_${name}.out
    
    #         herwig

    python transferKernelUnique.py -i ${indir}JJ_2016_nonRes_${cat}.root --sample herwig --year 2016 -p x --pdfIn ${indir}JJ_2016_nonRes_3D_NP.root | tee KernelTransf_2016_${cat}_herwig_x_${name}.out
    #python transferKernel.py -i ${indir}JJ_2016_nonRes_${cat}.root --sample herwig --year 2016 -p y --pdfIn ${indir}JJ_2016_nonRes_3D_NP.root | tee KernelTransf_2016_${cat}_herwig_y_${name}.out
    #python transferKernel.py -i ${indir}JJ_2016_nonRes_${cat}.root --sample herwig --year 2016 -p z --pdfIn ${indir}JJ_2016_nonRes_3D_NP.root | tee KernelTransf_2016_${cat}_herwig_z_${name}.out

    #        merge 
    python transferKernelUnique.py -i ${indir}JJ_2016_nonRes_${cat}.root --sample pythia --year 2016 -p z --pdfIn ${indir}JJ_2016_nonRes_3D_NP.root --merge
    
    mv postfit_qcd/PostFit_*.* $postfitdir

    #     control plots

    #python Projections3DHisto.py --mc ${indir}JJ_2016_nonRes_${cat}.root,nonRes -k save_new_shapes_2016_pythia_${cat}_3D.root,histo -o control-plots-QCD_pythia_signals_${cat}_${name}/       
    #python Projections3DHisto.py --mc ${indir}JJ_2016_nonRes_${cat}_altshapeUp.root,nonRes -k save_new_shapes_2016_herwig_${cat}_3D.root,histo -o control-plots-QCD_herwig_signals_${cat}_${name}/
    #python Projections3DHisto.py --mc ${indir}JJ_2016_nonRes_${cat}_altshape2.root,nonRes -k save_new_shapes_2016_madgraph_${cat}_3D.root,histo -o control-plots-QCD_madgraph_signals_${cat}_${name}/
    python Projections3DHisto_HPHP.py --mc ${indir}JJ_2016_nonRes_${cat}.root,nonRes -k save_new_shapes_2016_pythia_${cat}_3D.root,histo -o control-plots-QCD_pythia_signals_${cat}_${name}/       
    python Projections3DHisto_HPHP.py --mc ${indir}JJ_2016_nonRes_${cat}_altshapeUp.root,nonRes -k save_new_shapes_2016_herwig_${cat}_3D.root,histo -o control-plots-QCD_herwig_signals_${cat}_${name}/
    python Projections3DHisto_HPHP.py --mc ${indir}JJ_2016_nonRes_${cat}_altshape2.root,nonRes -k save_new_shapes_2016_madgraph_${cat}_3D.root,histo -o control-plots-QCD_madgraph_signals_${cat}_${name}/
    
done    
