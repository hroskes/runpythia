#!/bin/bash
cd /afs/cern.ch/user/h/hroskes/work/public/ttHdecay/pythiatest/CMSSW_7_1_14/src/OUTDIR
#cd /afs/cern.ch/work/d/dsperka/Run2MC/CMSSW_7_1_14/src/OUTDIR/
eval `scram runtime -sh`
cmsRun CFGFILE
