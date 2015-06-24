#! /bin/bash
dir=/afs/cern.ch/user/h/hroskes/work/public/forMeng/ttH/
hepdir=/scratch0/hep/hroskes/pythiatest/     #comment out to store on lxplus, in CMSSW_7_1_14/src/(filename without .lhe)/
hepusername=hroskes

if ! [ $3 ]; then
    echo "    ./hadronize.sh file.lhe events events_per_job"
    exit 1
fi &&

ls $1 > /dev/null &&

CWD=$(pwd -P) &&
cd $dir &&

if ! [ -d CMSSW_7_1_14 ]; then
    scram p CMSSW CMSSW_7_1_14
fi &&

cd CMSSW_7_1_14/src &&
eval $(scram ru -sh) &&
mkdir -p Configuration/GenProduction/python/ThirteenTeV/ &&
if ! [ -f Configuration/GenProduction/python/ThirteenTeV/Hadronizer_TuneCUETP8M1_13TeV_generic_LHE_pythia8_Tauola_cff.py ]; then
    cp ../../script/Hadronizer_TuneCUETP8M1_13TeV_generic_LHE_pythia8_Tauola_cff.py Configuration/GenProduction/python/ThirteenTeV/Hadronizer_TuneCUETP8M1_13TeV_generic_LHE_pythia8_Tauola_cff.py
fi &&
cd Configuration &&
scram b &&
cd .. &&

cd $CWD &&
ln -fs $(readlink -f "$1") "$(cd -)/" &&

cd - > /dev/null &&

echo "Hadronizing $1..." &&

a=$(basename $1) &&
lhefile=$a &&
GENfile=${a/.lhe/-GEN.root} &&
GENcfg=${a/.lhe/-GEN_cfg.py} &&
GENSIMfile=${a/.lhe/-GEN-SIM_py8.root} &&
GENSIMcfg=${a/.lhe/-GEN-SIM_py8_cfg.py} &&
GENSIMcfgtemplate=${GENSIMcfg/_cfg.py/_cfg_template.py} &&
jobdir=${a/.lhe/}
cmsDriver.py step1 --filein file:$lhefile --fileout file:$GENfile --mc --eventcontent LHE --datatier GEN --conditions MCRUN2_71_V1::All --step NONE --python_filename $GENcfg --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n $2 &&
cmsRun $GENcfg &&
cmsDriver.py Configuration/GenProduction/python/ThirteenTeV/Hadronizer_TuneCUETP8M1_13TeV_generic_LHE_pythia8_Tauola_cff.py --filein file:$GENfile --fileout file:$GENSIMfile --mc --eventcontent RAWSIM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --conditions MCRUN2_71_V1::All --step GEN,SIM --magField 38T_PostLS1 --python_filename $GENSIMcfg --no_exec -n $3 &&
if [ $2 -ne $3 ]; then
    ../../script/converttotemplate.sh $GENSIMcfg &&
    echo '
        #!/bin/bash
        tmpdir=$(mktemp -d)
        cd $tmpdir
        scram p CMSSW CMSSW_7_1_14
        cd CMSSW_7_1_14/src
        eval $(scram ru -sh)
        mkdir '"$jobdir
        cd $jobdir
        ln -fs $dir/CMSSW_7_1_14/src/$GENfile .
        ln -fs $dir/CMSSW_7_1_14/src/$jobdir/CFGFILE .
        cmsRun CFGFILE
        if [ $hepdir ]; then
            rsync -az ${GENSIMfile/.root/_JOBNUMBER.root} $hepusername""@hep.pha.jhu.edu:$hepdir
        else
            rsync -az ${GENSIMfile/.root/_JOBNUMBER.root} $dir/CMSSW_7_1_14/src/$jobdir
        fi
    " > template.sh &&
    python ../../script/submitJobs.py $GENSIMcfgtemplate $jobdir $2 $3
else
    cmsRun $GENSIMcfg
fi
