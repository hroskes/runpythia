infile=$1
outfile=${infile/_cfg.py/_cfg_template.py}
cat $infile | sed "s/\(fileNames = .*root')\)/\1,/" |
              sed "/fileNames = .*root'),/a eventsToProcess = cms.untracked.VEventRange('1:FIRSTEVENT-1:LASTEVENT')" |
              sed "s/^events/    events/" |
              sed "s/\(fileName = .*\)\.root'),/\1_JOBNUMBER.root'),/" > $outfile
