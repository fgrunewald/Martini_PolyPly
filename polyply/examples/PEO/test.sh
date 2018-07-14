rm -r test
mkdir test
cd test
cp ../PEO.top ./
cp ../PEO_OH.top ./

polyply -polymer PEO.martini.2 -name PEO -n_mon 10 -o PEO.itp
polyply -p PEO.top -o PEO.gro -env vac -n_mon 10 -name PEO


polyply -polymer PEO.martini.2 -name PEO -n_mon 10 -o PEO_OH.itp -endgroup ../SP2.itp ../SP2.itp
polyply -p PEO_OH.top -o PEO_OH.gro -env vac -n_mon 10 -name PEO
