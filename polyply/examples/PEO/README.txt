Produce an itp of PEO grown a single chain with and without SP2 terminus.

# Unterminated PEO

polyply -polymer PEO.martini.2 -name PEO -n_mon 10 -o PEO.itp
polyply -p ../system.top -o out.gro -env vac -n_mon 10 -name PEO


# Terminated PEO

polyply -polymer PEO.martini.2 -name PEO -n_mon 10 -o PEO.itp -endgroup SP2.itp SP2.itp
polyply -p ../system.top -o out.gro -env vac -n_mon 10 -name PEO
