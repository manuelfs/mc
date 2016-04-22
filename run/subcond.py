#!/usr/bin/env python

#### Script to send jobs to condor at UCSD

import os,sys
import time
import glob

models = [
  "TChiHH_mChi-400_onlyGEN",
]
nevents = "5000"
total_jobs = 10


# Folder for executables
if not os.path.exists("run"):
  os.mkdir("run")

for model in models:
  #gp_path = "/nfs-7/userdata/manuelf/gridpacks/"+model+"_tarball.tar.xz"
  gp_path = "/nfs-7/userdata/manuelf/gridpacks/TChiHH_mChi-400_tarball.tar.xz"
  outdir = "/hadoop/cms/store/user/manuelf/private_74X/"+model+"/"
  if not os.path.exists(outdir):
    os.mkdir(outdir)

  exefile = "/home/users/manuelf/code/mc/run/run_"+model+".sh"
  
  for i in range(total_jobs):
      fname = "log_"+model+"_batch_"+str(i)
      subfile = "run/cond_"+model+"_batch_"+str(i)+".cmd"
      f = open(subfile,"w")
      f.write("Universe = grid\n")
      f.write("Grid_Resource = condor cmssubmit-r1.t2.ucsd.edu glidein-collector.t2.ucsd.edu\n")
      f.write("x509userproxy=/tmp/x509up_u31569\n")
      f.write("+DESIRED_Sites=\"T2_US_UCSD\"\n")
      f.write("Executable = "+exefile+"\n")
      f.write("arguments =  "+' '.join([nevents,str(i),outdir])+"\n")
      f.write("Transfer_Executable = True\n")
      f.write("should_transfer_files = YES\n")
      f.write("transfer_input_files = "+gp_path+"\n")
      f.write("Notification = Never\n")
      f.write("Log=run/"+fname+".log.$(Cluster).$(Process)\n")
      f.write("output=run/"+fname+".out.$(Cluster).$(Process)\n")
      f.write("error=run/"+fname+".err.$(Cluster).$(Process)\n")
      f.write("queue 1\n")
      f.close()
  
      cmd = "condor_submit " + subfile
      #print cmd
      os.system(cmd)

