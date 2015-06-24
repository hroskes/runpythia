To run: just change the top few lines of hadronize.sh

   ./hadronize.sh file.lhe events events_per_job

if events=events_per_job, then it is assumed to be a small number, and pythia is run interactively.  Otherwise, it will submit jobs to the batch queue specified at the top of hadronize.sh

If running in parallel, the output files can be copied to a specified directory on hep, or left on lxplus in CMSSW_7_1_14/src/(file minus .lhe)/ by commenting out the assignment of hepdir.
