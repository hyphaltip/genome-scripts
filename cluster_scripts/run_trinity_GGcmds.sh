#PBS -l nodes=1:ppn=4,mem=2gb -N TrinityJob -j oe

module load trinity-rnaseq/r2013-05-23
module load bowtie2/2.1.0
module load bowtie
module load gmap
module load samtools

INFILE=trinity_GG.cmds
N=$PBS_ARRAYID

if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then
 echo "Need a PBS_ARRAYID or cmdline number"
 exit;
fi

line=`head -n $N $INFILE  | tail -n 1`

$line
