#!/bin/bash

cd materials
python run_Matsys.py
cd .. 

cd fiducial
python run_MCSbyxml.py
cd ..

#cd alignment
#python run_alignmentsys.py
#cd .. 

#cd TOFScan
#python run_TOFAllscan.py
#cd .. 

#cd TOFsys
#python run_TOFsys.py
#cd .. 

#cd angdef
#python run_angdefsys.py
#cd .. 





