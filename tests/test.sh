rm -f test.itp
python3.5 ../polyply.py -n_mon 100 -itp ../monomer_itps/PS_Monticelli_mon.itp -o "test.itp" -name PS100
python3.5 compare.py
