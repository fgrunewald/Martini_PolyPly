# Martini_PolyPly

## CURRENTLY UNDER CONSTRUCTION 

## Functionality 
PolyPly can be used to **generate GROMACS itp files** of polymers  and **starting confomrations** for polymer systems in principle of any type of force-field (FF). It has mainly been developed and tested for [MARTINI polymers](http://www.cgmartini.nl/index.php/force-field-parameters/polymers). At the moment it only supports a limited number of GROMACS bonded interactions, which limits it's applicability to polymers, which use no other bonded interactions than those included. As time progresses they will be updated and completed. The tool can also be used to generate intial structure files for more or less complex systems. We aim to offer the functionality for creating melts, single chains in solvent and multiple chains in solvent. The tool is in principle also applicable for any other GROMACS based force-fields, but again it has mainly been tested and optimized for MARTINI polymers. Besides offering a partical tool for easy simulation of polymer systems within GROMACS, the aim is also to make MD simulations more reproducable by offering a 

## Usage
### Itp file generation
```
polyply -itp monomer_A.itp monomer_B.itp ... -n_mon int_A int_B ... -o polymer.itp
```
#### Format monomer.itp
The itp generation tool essentially takes a monomer itp file, which uses the standard GROMACS itp file format. There are only two things to consider. On the one hand the file has to contain all bonded paramter of one repeat unit including those with the following repeat units, if there are any. On the other hand comments on the same line as any section indicator or term is not permitted. Any line with the ';' character will be treated as comment. In any case always make sure for a short chain that the output contains everything you would expect. Some example monomer.itp files can be found in the monomers directory. **Feel free to commit new monomer itp files to the monomer directory**  

#### 

### Initial System generation
```
polyply -sys [sol, melt] -s polymer.top -o system.gro -flory [good, theta, bad] -n_chains 10<int<100
```

## MARTINI Polymers 
In the Monomer-Itps directory monomer ITP files can be found. For the MARTINI itps make sure to use the correct force-field version. For most cases the force-field deviates slightly from the standard MARTINI version 2.x. The correct versions for some MARTINI polymers are listed below: 
* [PS](http://www.cgmartini.nl/images/applications/polymers/martini_v2.1_PS.itp)
* [PEO_Monticelli_et_al](http://perso.ibcp.fr/luca.monticelli/MARTINI/index.html)
* [PEO_Lee_et_al](http://www.cgmartini.nl/images/parameters/ITP/martini_v2.2.itp) make sure to change the SN0-SN0 interactions and the SN0-P4 interactions as detailed [here](http://www.cgmartini.nl/index.php/force-field-parameters/polymers)
* [P3HT_Alessandri_et_al](http://www.cgmartini.nl/images/parameters/ITP/martini_v2.2.itp)

## Authors

## License

This project is licensed under the GNU general public license - see the [LICENSE](LICENSE) file for details.
