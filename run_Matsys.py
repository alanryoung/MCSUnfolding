import libxml2
import os, subprocess
from ROOT import TRandom3, TMath


def updatesysval(infile, outfile, name, value):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/sys"):
        if node.prop("name").find(name) >= 0:
            node.setProp("value", str(value))
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()

def updatecutval(infile, outfile, name, value):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/cuts"):
        if node.prop("name").find(name) >= 0:
            node.setProp("value", str(value))
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()

def updatefilename(infile, outfile, newname):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/file"):
        if node.prop("id").find("outfile") >= 0:
            node.setProp("name", newname)
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()

def submitUnfolding(xmlfile):
    execdir = "/data/neutrino02/rbayes/MICE/MCS_selection/Unfolding/v2.5.0"
    cmd = [os.path.join(execdir,"MCSUnfolding"), xmlfile] 
    proc = subprocess.Popen(cmd)
    proc.wait()

def print_batch_submission(configfile, settings):
    batch_data = '''
#PBS -N %(name)s
#PBS -q medium6
#PBS -l walltime=5:59:00,mem=3000Mb
#PBS -e %(working_dir)s/logs/%(name)s.err
#PBS -o %(working_dir)s/logs/%(name)s.out

cd %(working_dir)s
. %(maus_root_dir)s/local_env.sh
%(execcmd)s %(xmlfile)s
'''% settings
    
    outfile = open(configfile, 'w+')
    outfile.write(batch_data)
    outfile.close()

def submit_to_batch(xmlfile):
    execdir = "/data/neutrino02/rbayes/MICE/MCS_selection/Unfolding/v2.5.0"
    execcmd = os.path.join(execdir,"MCSUnfolding")
    maus_root_dir = "/data/neutrino02/rbayes/MICE/MCS_selection/Unfolding/v2.5.0"
    working_dir = os.getcwd()
    name = xmlfile[:-4]
    settings = {"xmlfile":xmlfile, "name":name, "working_dir":working_dir, "maus_root_dir":maus_root_dir, "execcmd":execcmd}
    batch_file = os.path.join(working_dir,  name + ".job")
    print_batch_submission(batch_file, settings)
    cmd = ['qsub', batch_file]
    q = subprocess.Popen(cmd)
    q.wait()


llim_172 = 1.
ulim_172 = 28.9

llim_200 = 1.
ulim_200 = 27.9

llim_240 = 1.
ulim_240 = 27.3

for i in [-10,-3,-2,-1,1,2,3,10]:

    updatefilename("material_0/LiHMu_3172_0.xml",
                   "material_0/LiHMu_3172_res_lim"+str(llim_172 + i*0.0109)+".xml",
                   "material_0/LiHMu_3172_res_lim"+str(llim_172 + i*0.0109)+".root")
    updatesysval("material_0/LiHMu_3172_res_lim"+str(llim_172 + i*0.0109)+".xml",
                 "material_0/LiHMu_3172_res_lim"+str(llim_172 + i*0.0109)+".xml",
                 "resX", llim_172 + i*0.0109)
    updatesysval("material_0/LiHMu_3172_res_lim"+str(llim_172 + i*0.0109)+".xml",
                 "material_0/LiHMu_3172_res_lim"+str(llim_172 + i*0.0109)+".xml",
                 "resY", llim_172 + i*0.0109)
    submit_to_batch("material_0/LiHMu_3172_res_lim"+str(llim_172 + i*0.0109)+".xml")

    updatefilename("material_0/LiHMu_3200_0.xml",
                   "material_0/LiHMu_3200_res_lim"+str(llim_200 + i*0.0109)+".xml",
                   "material_0/LiHMu_3200_res_lim"+str(llim_200 + i*0.0109)+".root")
    updatesysval("material_0/LiHMu_3200_res_lim"+str(llim_200 + i*0.0109)+".xml",
                 "material_0/LiHMu_3200_res_lim"+str(llim_200 + i*0.0109)+".xml",
                 "resX", llim_200 + i*0.0109)
    updatesysval("material_0/LiHMu_3200_res_lim"+str(llim_200 + i*0.0109)+".xml",
                 "material_0/LiHMu_3200_res_lim"+str(llim_200 + i*0.0109)+".xml",
                 "resY", llim_200 + i*0.0109)
    submit_to_batch("material_0/LiHMu_3200_res_lim"+str(llim_200 + i*0.0109)+".xml")

    
    updatefilename("material_0/LiHMu_3240_0.xml",
                   "material_0/LiHMu_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                   "material_0/LiHMu_3240_res_lim"+str(llim_240 + i*0.0109)+".root")
    updatesysval("material_0/LiHMu_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                 "material_0/LiHMu_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                 "resX", llim_240 + i*0.0109)
    updatesysval("material_0/LiHMu_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                 "material_0/LiHMu_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                 "resY", llim_240 + i*0.0109)
    submit_to_batch("material_0/LiHMu_3240_res_lim"+str(llim_240 + i*0.0109)+".xml")

    '''
    updatefilename("material_0/XePion_3240.xml",
                   "material_0/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                   "material_0/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".root")
    updatesysval("material_0/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                 "material_0/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                 "resX", llim_240 + i*0.0109)
    updatesysval("material_0/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                 "material_0/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                 "resY", llim_240 + i*0.0109)
    submit_to_batch("material_0/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".xml")

'''
