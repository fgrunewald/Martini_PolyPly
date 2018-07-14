#!/bin/bash

rm -r test
mkdir test
cd test 

cp ../martini.itp ./
cp ../*.itp ./
cp ../system.top ./

polyply -polymer PEL_lipid.martini.2 PEO.martini.2 -n_mon 1 15 -o 5_PEL.itp -name PEL -links ../links.itp -endgroup ../dum.itp ../SP2.itp
insane -o bilayer.gro -p topol.top -l DPPC:20 -l DOPE:1 -salt 0 -sol W -x 4 -y 4 -z 15

gmx_mpi grompp -f ../min.mdp -p topol.top -c bilayer.gro
gmx_mpi mdrun

polyply -env bilayer -lipid DOPE -sys confout.gro -p system.top -name PEL -sol W -o PEGylated_bilayer.gro
