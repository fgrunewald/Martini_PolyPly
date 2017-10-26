from Martini_PolyPly.itp_tool.itp_I import *

ref=read_itp("PS100_Amap.itp")
new=read_itp("test_PS.itp")
ref_PEO=read_itp("PEO_Lee_37.itp")
new_PEO=read_itp("test_PEO.itp")

#ref_P3HT=read_itp("P3HT_Alessandri.itp")
#new_P3HT=read_itp("test_P3HT.itp")

bonds, angles, constraints, dihedrals, exclusions, virtual = [], [], [], [], [], []

for i, j in zip(ref["bonds"], new["bonds"]):
    for a, b in zip(i,j):
        bonds += [a==b]

for i, j in zip(sorted(ref["angles"]), sorted(new["angles"])):
     for a, b in zip(i,j):
         angles += [a==b]

for i, j in zip(sorted(ref["constraints"]), sorted(new["constraints"])):
     for a, b in zip(i,j):
         angles += [a==b]

for i, j in zip(ref_PEO["bonds"], new_PEO["bonds"]):
    for a, b in zip(i,j):
        bonds += [a==b]

for i, j in zip(sorted(ref_PEO["angles"]), sorted(new_PEO["angles"])):
     for a, b in zip(i,j):
         angles += [a==b]

for i, j in zip(sorted(ref_PEO["dihedrals"]), sorted(new_PEO["dihedrals"])):
     for a, b in zip(i,j):
         dihedrals += [a==b]

for i, j in zip(ref_P3HT["bonds"], new_P3HT["bonds"]):
    for a, b in zip(i,j):
        bonds += [a==b]

for i, j in zip(sorted(ref_P3HT["angles"]), sorted(new_P3HT["angles"])):
     for a, b in zip(i,j):
         angles += [a==b]

for i, j in zip(sorted(ref_P3HT["dihedrals"]), sorted(new_P3HT["dihedrals"])):
     for a, b in zip(i,j):
         dihedrals += [a==b]

for i, j in zip(sorted(ref_P3HT["virtual_sitesn"]), sorted(new_P3HT["virtual_sitesn"])):
     for a, b in zip(i,j):
         virtual += [a==b]

for i, j in zip(sorted(ref_P3HT["exclusions"]), sorted(new_P3HT["exclusions"])):
     for a, b in zip(i,j):
         exclusions += [a==b]

for i, j in zip(sorted(ref_P3HT["constraints"]), sorted(new_P3HT["constraints"])):
     for a, b in zip(i,j):
         constraints += [a==b]

print("Bonds:", all(bonds))
print("Angles:", all(angles))
print("Constraints:", all(constraints))
print("dihedrals:", all(dihedrals))
print("virtuals:", all(dihedrals))
print("exclusions:", all(dihedrals))

