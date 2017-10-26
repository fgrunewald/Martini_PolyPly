rm -f test_PS.itp
rm -f test_PEO.itp
python3.5 ../polyply.py -n_mon 100 -itp ../monomer_itps/PS_Monticelli_mon.itp -o "test_PS.itp" -name PS100
python3.5 ../polyply.py -n_mon 37 -itp ../monomer_itps/PEO_Lee_mon.itp -o "test_PEO.itp" -name PEO37
python3.5 compare.py
