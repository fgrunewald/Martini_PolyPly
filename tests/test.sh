rm -f test_PS.itp
rm -f test_PEO.itp
rm -r test_P3HT.itp
python3.5 ../polyply.py -n_mon 48 -polymer P3HT -o "test_P3HT.itp" -name P3HT48
python3.5 ../polyply.py -n_mon 100 -polymer PS -o "test_PS.itp" -name PS100
#python3.5 ../polyply.py -n_mon 37 -polymer old/PEO_Lee_mon.itp -o "test_PEO.itp" -name PEO37
python3.5 compare.py
error_A=$(python3.5 ../polyply.py -n_mon 1 -itp false_format_dope_I.itp -o "test_PEO.itp" -name DOPE | grep -c Error)
error_B=$(python3.5 ../polyply.py -n_mon 1 -itp false_format_dope_II.itp -o "test_PEO.itp" -name DOPE | grep -c Error)
let "error = error_A + error_B"

if [ ${error} == 2 ]
then
   echo "Error handeling: True"
else
   echo "Error handeling: False"
fi
