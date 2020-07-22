#!/bin/csh -f

if ( $#argv < 4 ) then
   echo Need to give ntoys per dir, ndirs, signal strength and the masspoint mp for the corresponding datacard
   exit
endif
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0- PREPARE SCRIPT TO SUBMIT 
echo "-------------------------------------------------------------"
#echo "0) WARNING: killing all OLD jobs by " $user 
#condor_rm -all

set ntoys = $1
set ndirs = $2
set mu = $3
set mp = $4
set index_toy = $5
set index_fit = $6


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1- SUBMIT JOBS
echo "-------------------------------------------------------------"

set datacard_dir = "/afs/desy.de/user/v/vagneria/cms/analysis2017/2017data/TTBar/Combine/CMSSW_8_1_0/src/Analysis/Models/test/bias_study/datacards"
set i=0
set seed=112456  # initialise seed
while ( $i < $ndirs )

   set exedir = "HTC_$i"
   if ( -d $exedir ) then
      echo "Similar jobs were already submitted. Move ot remove directories "
      exit
   endif

    mkdir -p $exedir
    cd $exedir
    cp ${datacard_dir}/datacard_bias_m${mp}_pdf_index_*.txt . #copy all datacard for some pdf

    echo "Write executable for each dir"
    echo $mp
    cat > exe.sh <<EOF

cd ${datacard_dir}
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
cd -

#Toy Generation
combine datacard_bias_m${mp}_pdf_index_${index_toy}.txt -M GenerateOnly -s $seed --setParameters pdf_index=0 --toysFrequentist -t $ntoys --expectSignal $mu --saveToys -m $mp --freezeParameters pdf_index #-S 0 

#Estimate bias with alternative function
combine datacard_bias_m${mp}_pdf_index_${index_fit}.txt -M FitDiagnostics -s $seed --setParameters pdf_index=0 --toysFile higgsCombineTest.GenerateOnly.mH${mp}.${seed}.root -m $mp -t $ntoys --rMin -10 --rMax 10 --freezeParameters pdf_index --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --X-rtd ADDNLL_RECURSIVE=0 --X-rtd FITTER_NEW_CROSSING_ALGO --cminDefaultMinimizerTolerance 0.01 --robustFit 1  # -v2

EOF
    chmod u+x exe.sh
    ../../htc_sub.sh exe.sh
#    sleep 0.01
    cd -
    @ i++
    @ seed++  # upgrade seed to avoid identical toys 
end

 


