rm -r test
mkdir test
cd test
polyply -p ../system.top -sys ../BENZ.gro -o out.gro -env sol -n_mon 5 -sol BENZ -name PS
