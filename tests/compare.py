from Martini_PolyPly.itp_tool.itp_I import *

ref=read_itp("PS100_Amap.itp")
new=read_itp("test.itp")
bonds, angles, constraints = [], [], []

for i, j in zip(ref["bonds"], new["bonds"]):
    for a, b in zip(i,j):
        bonds += [a==b]

for i, j in zip(sorted(ref["angles"]), sorted(new["angles"])):
     for a, b in zip(i,j):
         angles += [a==b]

for i, j in zip(sorted(ref["constraints"]), sorted(new["constraints"])):
     for a, b in zip(i,j):
         angles += [a==b]

print("Bonds:", all(bonds))
print("Angles:", all(angles))
print("Constraints:", all(constraints))
