import os, commands, time, sys

try:
    cfg_template = sys.argv[1]
    outdir = sys.argv[2]
    nevents = int(sys.argv[3])
    eventsperjob = int(sys.argv[4])
    njobs = (nevents + eventsperjob - 1) / eventsperjob
    queue = sys.argv[5]
except IndexError:
    print "Not enough arguments!"
    print "Usage:"
    print "    python submitJobs.py mytemplate_cfg.py outdir nevents eventsperjob queue"
    sys.exit(1)

resubmit = False
testing = False

csh_template = 'template.sh'

if (not os.path.exists(outdir)):
    os.system('mkdir '+outdir)

startDir = os.getcwd()
os.chdir(os.path.join(startDir,outdir))

def processCmd(cmd, quite = 0):
    if testing and "bsub" in cmd:
        print cmd
        return ""
    status, output = commands.getstatusoutput(cmd)
    if (status !=0 and not quite):
        print 'Error in processing command:\n   ['+cmd+']'
        print 'Output:\n   ['+output+'] \n'
    return output


njobsresubmitted = 0
for job in range(1,njobs+1):

    str_job = str(job)
    cfg_job = cfg_template.replace('_template','_'+str(job))
    csh_job = csh_template.replace('template','job_'+str(job))

    if (resubmit):
            
        jobexists=False
        with open('../existing_jobs.txt', 'r') as f:
            for line in f:
                if line.rstrip()==str_job: jobexists=True

        jobfailed=False
        with open('../failed_jobs.txt', 'r') as f:
            for line in f:
                if line.rstrip()==str_job: jobfailed=True

        if (jobexists and (not jobfailed)): continue
        
        print str_job,'exists?',jobexists,'failed?',jobfailed
        print 'resubmitting job',str(job)
        njobsresubmitted +=1
        cmd = 'bsub -q 1nh -J '+cfg_job+' '+csh_job
        output = processCmd(cmd)
        while ('error' in output):
            time.sleep(1.0);
            output = processCmd(cmd)
            if ('error' not in output):
                print output
                print 'Submitting after retry - job '+str(job)
            print output
        print output
        continue


    output = processCmd('cp ../'+cfg_template+' '+cfg_job)
    output = processCmd("sed -i 's~JOBNUMBER~"+str(job)+"~g' "+cfg_job)

    firstevent = nevents/njobs*(job-1)+1
    lastevent = nevents/njobs*(job)

    output = processCmd("sed -i 's~FIRSTEVENT~"+str(firstevent)+"~g' "+cfg_job)
    output = processCmd("sed -i 's~LASTEVENT~"+str(lastevent)+"~g' "+cfg_job)

    output = processCmd('cp ../'+csh_template+' '+csh_job)
    output = processCmd("sed -i 's~OUTDIR~"+outdir+"~g' "+csh_job)
    output = processCmd("sed -i 's~CFGFILE~"+cfg_job+"~g' "+csh_job)
    output = processCmd("sed -i 's~JOBNUMBER~"+str(job)+"~g' "+csh_job)

    #if (job>1): continue

    if cfg_job in processCmd("bjobs -w -q " + queue):
        continue

    #output = processCmd('bsub -q 1nh -J '+cfg_job+' '+csh_job)
    print 'submitting job',str(job)
    cmd = 'bsub -q ' + queue + ' -J '+cfg_job+' '+csh_job
    output = processCmd(cmd)
    while ('error' in output):
        time.sleep(1.0);
        output = processCmd(cmd)
        if ('error' not in output):
            print output
            print 'Submitting after retry - job '+str(job)
    print output

if (resubmit): print njobsresubmitted,'jobs resubmitted'
os.chdir(startDir)

