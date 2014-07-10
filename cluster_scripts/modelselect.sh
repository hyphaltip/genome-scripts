#PBS -l nodes=1:ppn=8 -q js -N modelselection 
module load RAxML

ProteinModelSelection.pl in.phy >& modelselect.log
