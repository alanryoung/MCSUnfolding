#!/usr/bin/env python
import os, subprocess

def main():
    number_of_mcfiles=972
    file=open("reduce_MC_empty_FG.sh","w")
    file.write("#!/bin/bash\n\n")
    file.write("/data/mice/bin/global data_zero_3200_FG")
    for i in range(number_of_mcfiles):
#        print(i)
        numstr=str(i+1).zfill(5)
        linewrite=" \\\n                       ./Empty_absorber/200MeV/MCproduction_run_123/"+numstr+"_mc/"+numstr+"_sim.root"        
        file.write(linewrite)
    file.close()
    execdir = "."
    cmd = [os.path.join(execdir,"reduce_MC_empty_FG.sh")]
    print(cmd) 
    proc = subprocess.Popen(cmd)
    proc.wait()
   

if __name__ == "__main__":
    main()

