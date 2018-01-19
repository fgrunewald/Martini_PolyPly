# Martini_PolyPly
# NOT UP TO DATE !!!
# Please wait for release!

## Functionality 
PolyPly can be used to **generate GROMACS itp files** of polymers from a monomer itp file and **produce polymer molecules** from topology files. In principle the program can be used with any type of force-field (FF) as long as the files are in GROMACS format. It has mainly been developed and tested for [MARTINI polymers](http://www.cgmartini.nl/index.php/force-field-parameters/polymers). At the moment it only supports a limited number of GROMACS bonded interactions, which limits it's applicability to polymer with these bonded interactions. As time progresses they will be updated and completed. The tool can also be used to **generate intial structure files for more or less complex systems**. At the moment it supports growing polymers on DOPE lipids to form for example PEGylated bilayers and into solvated systems of any kind. Read carefully the usage section, since the tool uses strict file formats. 

## Itp file generation
To generate a polymer itp file the tool offers two options. Either polymers from the library (mostly MARTINI) can be generated or one can generate polymers from custom monomer-itp files. 

#### from library
To generate a polymer, which is part of the library simply execute the following command:
```
polyply -polymer [name of polymer] -n_mon [number of repeat units] -o [name of outfile] -name [name of polymer]
```
The name option gives the possibility to give another than the default name to the polymer. The names of the polymers in the library can be found by executing the command:
```
polyply -lib
```
#### the environment option

#### the bilayer option

#### generating block-copolymers
If you want to generate a block-copolymer, you can do this by supplying multiple names and monomer repeat units. For example to generate a block-copolymer of 100 repeat units PS and 100 repeat units PEO the input would be:
```
polyply -itp PS_monomer.itp PEO_monomer.itp -n_mon 100 100 -o PS-b-PEO.itp -linker links.itp
```
The above command generates a block-copolymer itp for PS and PEO each block with a length of 100 monomers. The linkage between atom 100 and 101, has to be specified in the 'links.itp' file. The format is the same as for a normal itp. The only difference is that there is no '[ moleculetype ]' and '[ atoms ]' directive. It only contains the bonded parameters corresponding to the link. Of course one needs to make sure the atom numbering is correct. In the case where we just want a single bond to link the PS and PEO block the linker.itp could look like this:
```
[ bonds ]
; atom1  atom2    ref     fc
   100   101     0.33    7000
```
#### custom itp files monomer.itp
The itp generation tool essentially takes a monomer itp file, which uses the standard GROMACS itp file format. The itp-file has to contain all bonded parameters of one repeat unit including those with the following repeat units, if there are any. Always make sure for a short chain that the output contains everything you would expect. Some example monomer.itp files can be found in the monomer_itps directory. **Feel free to propose new monomer itp files to be included in the monomer directory.**  

#### Block-Copolymers and non-Block-Copolymers (in the beta stage! not everything works yet!)
The tool can also be used to generate block-copolymer itp files from two or more separate monomer itp files. In this case the sequence of blocks and their length has to match the one supplied in the first two arguments above. For example to generate a block-copolymer of 100 repeat units PS and 100 repeat units PEO the input would be:
```
polyply -itp PS_monomer.itp PEO_monomer.itp -n_mon 100 100 -o PS-b-PEO.itp -linker links.itp
```
The above command generates a block-copolymer itp for PS and PEO each block with a length of 100 monomers. The linkage between atom 100 and 101, has to be specified in the 'links.itp' file. The format is the same as for a normal itp. The only difference is that there is no '[ moleculetype ]' and '[ atoms ]' directive. It only contains the bonded parameters corresponding to the link. Of course one needs to make sure the atom numbering is correct. In the case where we just want a single bond to link the PS and PEO block the linker.itp could look like this:
```
[ bonds ]
; atom1  atom2    ref     fc
   100   101     0.33    7000
```
#### 

### Initial System generation
```
polyply -sys [sol, melt] -s system.top -o system.gro -flory [good, theta, bad] -n_chains 10<int<100
```
This function is currently in the development stage and not yet available.

## Authors

## License

This project is licensed under the GNU general public license - see the [LICENSE](LICENSE) file for details.
