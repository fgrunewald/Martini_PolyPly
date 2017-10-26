# Martini_PolyPly

## Functionality 
PolyPly can be used to **generate GROMACS itp files** of polymers  and **starting confomrations** for polymer systems in principle of any type of force-field (FF). It has mainly been developed and tested for [MARTINI polymers](http://www.cgmartini.nl/index.php/force-field-parameters/polymers). At the moment it only supports a limited number of GROMACS bonded interactions, which limits it's applicability to polymers, which use no other bonded interactions than those included. As time progresses they will be updated and completed. The tool can also be used to generate intial structure files for more or less complex systems. We aim to offer the functionality for creating melts, single chains in solvent and multiple chains in solvent. The tool is in principle also applicable for any other GROMACS based force-fields, but again it has mainly been tested and optimized for MARTINI polymers. Besides offering a partical tool for easy simulation of polymer systems within GROMACS, the aim is also to make MD simulations more reproducable by offering an easy to use robust tool.

## Usage
### Itp file generation for simple polymer
```
polyply -itp [monomer.itp] -n_mon [number_of_repeat_units] -o [outname] 
```
#### Format monomer.itp
The itp generation tool essentially takes a monomer itp file, which uses the standard GROMACS itp file format. The itp-file has to contain all bonded paramters of one repeat unit including those with the following repeat units, if there are any. Always make sure for a short chain that the output contains everything you would expect. Some example monomer.itp files can be found in the monomer_itps directory. **Feel free to propose new monomer itp files to be included in the monomer directory.**  

#### Block-Copolymers and non-Block-Copolymers (in the beta stage! not everything works yet!)
The tool can also be used to generate block-copolymer itp files from two or more seperate monomer itp files. In this case the sequence of blocks and their length has to match the one supplied in the first two arguments above. For example to generate a block-copolymer of 100 repeat units PS and 100 repeat units PEO the input would be:
```
polyply -itp PS_monomer.itp PEO_monomer.itp -n_mon 100 100 -o PS-b-PEO.itp 
```
In the above case the programm will put a bond of force constant 3000 kj mol^-1 nm^-1 between the atoms of consecutive blocks. This force constant works for MARTINI block-copolymers. Other force-fields might require a different specification or one might want to include a special linker group. In that case please prepare them as seperate monomer.itp and put them in the sequence as shown below:   
```
polyply -itp PS_monomer.itp Linker_group.itp PEO_monomer.itp -n_mon 100 1 100 -o PS-b-PEO.itp 
```
In future we might supply more generic options. The same principle as above has to be applied for other types of non-block copolymers.
#### 

### Initial System generation
```
polyply -sys [sol, melt] -s system.top -o system.gro -flory [good, theta, bad] -n_chains 10<int<100
```
This function is currently in the development stage and not yet available.

## MARTINI Polymers 
In the Monomer-Itps directory monomer ITP files can be found. For the MARTINI itps make sure to use the correct force-field version. For most cases the force-field deviates slightly from the standard MARTINI version 2.x. The correct versions for some MARTINI polymers are listed below: 
* [PS](http://www.cgmartini.nl/images/applications/polymers/martini_v2.1_PS.itp)
* [PEO_Monticelli_et_al](http://perso.ibcp.fr/luca.monticelli/MARTINI/index.html)
* [PEO_Lee_et_al](http://www.cgmartini.nl/images/parameters/ITP/martini_v2.2.itp) make sure to change the SN0-SN0 interactions and the SN0-P4 interactions as detailed [here](http://www.cgmartini.nl/index.php/force-field-parameters/polymers)
* [P3HT_Alessandri_et_al](http://www.cgmartini.nl/images/parameters/ITP/martini_v2.2.itp)

## Authors

## License

This project is licensed under the GNU general public license - see the [LICENSE](LICENSE) file for details.
