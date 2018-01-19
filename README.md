# Martini_PolyPly
# NOT UP TO DATE !!!
# Please wait for release!

## Functionality 
PolyPly can be used to **generate GROMACS itp files** of polymers from a monomer itp file and **produce polymer molecules** from topology files. In principle the program can be used with any type of force-field (FF) as long as the files are in GROMACS format. It has mainly been developed and tested for [MARTINI polymers](http://www.cgmartini.nl/index.php/force-field-parameters/polymers). At the moment it only supports a limited number of GROMACS bonded interactions, which limits it's applicability to polymer with these bonded interactions. As time progresses they will be updated and completed. The tool can also be used to **generate intial structure files for more or less complex systems**. At the moment it supports growing polymers on DOPE lipids to form for example PEGylated bilayers and into solvated systems of any kind. Read carefully the usage section, since the tool uses strict file formats. 

## 1. Itp file generation
To generate a polymer itp file the tool offers two options. Either polymers from the library (mostly MARTINI) can be generated or one can generate polymers from custom monomer-itp files. 

#### a - from library
To generate a polymer, which is part of the library simply execute the following command:
```
polyply -polymer [name of polymer] -n_mon [number of repeat units] -o [name of outfile] -name [name of polymer]
```
The 'n_mon' option sets the number of repeat units (i.e. n monomers) and the name option gives the possibility to give another than the default name to the polymer. The names of the polymers in the library can be found by executing the command:
```
polyply -lib
```
**Feel free to propose new monomer itp files to be included in the monomer directory.**  
#### b - generating block-copolymers
The tool can also be used to generate a block-copolymers of a custom number of blocks. This can be done by supplying multiple names via the polymer directive and multiple numbers via the n_repeat diretive. For example to generate a block-copolymer of 100 repeat units of PS and 100 repeat units of PEO the input would be:
```
polyply -polymer PS PEO -n_mon 100 100 -o PS-b-PEO.itp -linker links.itp -name PS-b-PEO
```
The above command generates a block-copolymer itp for PS and PEO each block with a length of 100 monomers. The linkage between atom 100 and 101, has to be specified in an external file here called 'links.itp'. The format is the same as for a normal itp. The only difference is that there is no '[ moleculetype ]' and '[ atoms ]' directive. It only contains the bonded parameters corresponding to the link. Of course one needs to make sure the atom numbering is correct. In the case where we just want a single bond to link the PS and PEO block the linker.itp could look like this:
```
[ bonds ]
; atom1  atom2    ref     fc
   100   101     0.33    7000
```
#### c - custom itp files monomer.itp
PolyPly also takes custom monomer itp files, as long as they use the standard GROMACS itp file format. The itp-file has to contain all bonded parameters of one repeat unit including those with the following repeat units, if there are any. Always make sure that the output contains everything you would expect. Some example monomer.itp files can be found in the monomer_itps directory. The only difference to the library command is, that the itp file name is supplied via the -itp option as shown in the following scheme:
```
polyply -itp [name of polymer itp file] -n_mon [number of repeat units] -o [name of outfile] -name [name of polymer]
```
Note by supplying multiple files also custom block-copolymer can be generated. 
### 2. Initial Strucutre generation
PolyPly also offers the possibility to grow polymers into existing systems. To indicate that a structure needs to be generated supply the -grow flag with one of the system classifications defined as follows:

option  | system
----------------------------
vac     | vacuum
sol     | solution
bilayer | PEGylated lipid

```
polyply -grow [ vac, bilayer, sol] -s system.top -o final.gro -name  -n_chains 10<int<100
```
## Authors

## License

This project is licensed under the GNU general public license - see the [LICENSE](LICENSE) file for details.
