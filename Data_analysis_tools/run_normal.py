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

def updatedatafilename(infile, outfile, newname):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/file"):
        if node.prop("id").find("datafile") >= 0:
            node.setProp("name", newname)
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()

def updatetrainfilename(infile, outfile, newname):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/file"):
        if node.prop("id").find("trainfile") >= 0:
            node.setProp("name", newname)
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()

def updatemodelfilename(infile, outfile, newname):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/file"):
        if node.prop("id").find("modelfile") >= 0:
            node.setProp("name", newname)
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()

def submitUnfolding(xmlfile):
    execdir = "/data/mice/bin"
    cmd = [os.path.join(execdir,"MCSUnfolding"), xmlfile]
    print(cmd) 
    proc = subprocess.Popen(cmd)
    proc.wait()

def print_batch_submission(configfile, settings):
    testscript = '''#!/bin/bash
cd %(working_dir)s
. %(maus_root_dir)s/local_env.sh
%(execcmd)s %(xmlfile)s
'''% settings
    
    outfilename = '''%(shellscript)s'''% settings
    outfile0 = open(outfilename, 'w+')
    outfile0.write(testscript)
    outfile0.close()

    batch_data = '''
universe       = vanilla
executable     = %(shellscript)s
output         = %(working_dir)s/%(name)s.out
error          = %(working_dir)s/%(name)s.err             
log            = test.log
requirements   = OpSysAndVer == "CentOS7"

queue
'''% settings

    outfile = open(configfile, 'w+')
    outfile.write(batch_data)
    outfile.close()

def submit_to_batch(xmlfile):
    execdir = "/home/ppe/j/jnugent/workarea/Unfolding/fiducial"
    execcmd = os.path.join(execdir,"../MCSUnfolding")
    maus_root_dir = "/home/ppe/j/jnugent/workarea/Unfolding/fiducial"
    working_dir = os.getcwd()
    name = xmlfile[:-4]
    settings = {"shellscript":xmlfile[:-4]+".sh", "xmlfile":xmlfile, "name":name, "working_dir":working_dir, "maus_root_dir":maus_root_dir, "execcmd":execcmd}
    batch_file = os.path.join(working_dir,  name + ".job")
    print_batch_submission(batch_file, settings)
    #cmd = ['qsub', batch_file]
    cmd = ['condor_submit', batch_file]
    q = subprocess.Popen(cmd)
    q.wait()

submitUnfolding("LiHMu_3172.xml")
submitUnfolding("LiHMu_3200.xml")
submitUnfolding("LiHMu_3240.xml")

''' 
     updatefilename("XePion_3240.xml", "XePion_3240_R_"+str(i)+"_G_"+str(j)+".xml", "XePion_3240_R_"+str(i)+"_G_"+str(j)+".root")
     updatecutval("XePion_3240_R_"+str(i)+"_G_"+str(j)+".xml", "XePion_3240_R_"+str(i)+"_G_"+str(j)+".xml", "fid_rad", 100 + i*10)
     updatecutval("XePion_3240_R_"+str(i)+"_G_"+str(j)+".xml", "XePion_3240_R_"+str(i)+"_G_"+str(j)+".xml", "fid_grad", 0.005*j)
     submit_to_batch("XePion_3240_R_"+str(i)+"_G_"+str(j)+".xml")
'''
# submitUnfolding("LiHMu_3172.xml")
# submitUnfolding("LiHMu_3200.xml")
# submitUnfolding("LiHMu_3240.xml")

'''
rndm = TRandom3()
for i in range(30):
    updatesysval("LiHMu_3172.xml", "LiHMu_3172_"+str(i)+".xml", "alXUS", -0.356 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "alYUS",  1.13 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "alXDS", -0.126 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "alYDS",  0.118 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "thXUS", -0.0031 +  rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "thYUS", -0.005 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "thXDS", -0.003 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "thYDS", -0.0151 + rndm.Gaus(0, 0.004))
    updatefilename("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".root")
    submitUnfolding("LiHMu_3172_"+str(i)+".xml")


for i in range(30):
    updatesysval("LiHMu_3200.xml", "LiHMu_3200_"+str(i)+".xml", "alXUS", -0.356 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "alYUS",  1.13 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "alXDS", -0.126 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "alYDS",  0.118 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "thXUS", -0.0031 +  rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "thYUS", -0.005 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "thXDS", -0.003 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "thYDS", -0.0151 + rndm.Gaus(0, 0.004))
    updatefilename("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".root")
    submitUnfolding("LiHMu_3200_"+str(i)+".xml")


for i in range(30):
    updatesysval("LiHMu_3240.xml", "LiHMu_3240_"+str(i)+".xml", "alXUS", -0.356 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "alYUS",  1.13 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "alXDS", -0.126 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "alYDS",  0.118 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "thXUS", -0.0031 +  rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "thYUS", -0.005 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "thXDS", -0.003 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "thYDS", -0.0151 + rndm.Gaus(0, 0.004))
    updatefilename("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".root")
    submitUnfolding("LiHMu_3240_"+str(i)+".xml")
'''
