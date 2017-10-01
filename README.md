# Martini_PolyPly

## Functionality 
PolyPly can be used to **generate GROMACS itp files** of polymers  and **starting confomrations** for polymer systems in principle of any type of force-field (FF). It has mainly been developed and tested for [MARTINI polymers](http://www.cgmartini.nl/index.php/force-field-parameters/polymers). At the moment it only supports a limited number of GROMACS bonded interactions, which limits it's applicability to polymers, which use no other bonded interactions than those included. As time progresses they will be updated and completed. The tool can also be used to generate intial structure files for more or less complex systems. We aim to offer the functionality for creating melts, single chains in solvent and multiple chains in solvent. The tool is in principle also applicable for any other GROMACS based force-fields, but again it has mainly been tested and optimized for MARTINI polymers. Besides offering a partical tool for easy simulation of polymer systems within GROMACS, the aim is also to make MD simulations more reproducable by offering a 

## Usage
### Itp file generation
Make sure to read the usage section before using  It reads an monomer ITP file, which also has to contain ALL bonded interactions to other monomers. It then multiplies these terms up to x repeats and removes the all terms beyond the last monomer.
### Initial System generation

## MARTINI Polymers 
In the Monomer-Itps directory monomer ITP files can be found. For the MARTINI itps make sure to use the correct force-field version. For most cases the force-field deviates slightly from the standard MARTINI version 2.x. The correct versions for some MARTINI polymers are listed below: 
* [PS](http://www.cgmartini.nl/images/applications/polymers/martini_v2.1_PS.itp)
* [PEO_Monticelli_et_al](http://perso.ibcp.fr/luca.monticelli/MARTINI/index.html)
* [PEO_Lee_et_al](http://www.cgmartini.nl/images/parameters/ITP/martini_v2.2.itp) make sure to change the S0-S0 interactions and the S0-P4 interactions as detailed [here](http://www.cgmartini.nl/index.php/force-field-parameters/polymers)
* [P3HT_Alessandri_et_al](http://www.cgmartini.nl/images/parameters/ITP/martini_v2.2.itp)

## Authors

## License

This project is licensed under the GNU License - see the [LICENSE.md](LICENSE.md) file for details
