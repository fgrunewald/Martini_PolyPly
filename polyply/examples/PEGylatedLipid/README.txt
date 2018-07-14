++++++++++++++++++++++++++++++++++++++++++++++++ GROW PEGylated LIPIDS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

1. Generate a itp file for PEG attached to a DOPE lipid
   polyply -polymer PEL_lipid PEO -n_mon 1 5 -o 5_PEL.itp -name PEL -links links.itp -endgroup

2. Generate the lipid bilayer by insane and minimize
   insane -o bilayer.gro -p topol.top -l DPPC:20 -l DOPE:1 -salt 0 -sol W -x 4 -y 4 -z 15

3. Minimize the isane bilayer
   gmx grompp -f min.mdp -p topol.top -c bilayer.gro
   gmx mdrun

4. Grow the PEGylated lipid ontop
   polyply -env bilayer -lipid DOPE -sys confout.gro -p system.top -name PEL -sol W

5. Visulaize the output with vmd
   a) open out.gro and insert one space in the heading line and close.
   b) run "vmd out.gro" 
