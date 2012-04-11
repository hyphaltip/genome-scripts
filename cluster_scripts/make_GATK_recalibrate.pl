#!/usr/bin/perl -w
use strict;
use Env qw(HOME);
use File::Spec;
use Getopt::Long;

my $jobdir = File::Spec->catdir($HOME,'jobs');
my $reference_genome;
my $gatk_jarfile = '/opt/stajichlab/GATK/latest/GenomeAnalysisTK.jar';
my $picard = '/usr/local/java/common/lib/picard-tools';
my $sort = '/usr/local/java/common/lib/picard-tools';
my $javamem ="-Xmx6g";
my $threads = '4';
mkdir("$HOME/bigdata/tmp");
my $targetsp;
my $tmpdir  = "$HOME/bigdata/tmp/";
my $force = 0;
my $ctr = 'Cornell';
my $platform = 'illumina';
GetOptions(
    'r|ref:s'    => \$reference_genome,
    't|target:s' => \$targetsp,
    'jar|gatk:s' => \$gatk_jarfile,
    'job:s'      => \$jobdir,
    'threads:i'  => \$threads,
    'javamem:s'  => \$javamem,
    'f|force!'  => \$force,
    'tmpdir:s'  => \$tmpdir,
    'c|ctr|center:s' => \$ctr,
    'pl|platform:s' => \$platform,
);

if( ! defined $reference_genome )  {
    die "must provide reference genome with -r or --ref\n";
} 

if( ! defined $targetsp )  {
    die "must provide target species/strain abbreviation of the reference genome with -t or --target\n";
} 

$reference_genome = File::Spec->rel2abs($reference_genome);
my @bamfiles = @ARGV;
for my $bam ( @bamfiles ) {  
    $bam = File::Spec->rel2abs($bam);
    my (undef,$bdir,$fname) = File::Spec->splitpath($bam);
	warn("$bam, $fname\n");
    next unless ( $fname =~ /(\S+)\.$targetsp\.sam$/);
    my $base = $1;
	warn("base is $base\n");
	my ($strain) = split(/\./,$base);
#    next if $fname =~ /(sort|sort_ordered|dedup|realign)/;

#    $bdir =~ s/\/$//;
    open(my $jobfh => ">$jobdir/$base.gatk_realign.sh") || die $!;
    print $jobfh "#PBS -N $base.realn -l nodes=1:ppn=$threads,walltime=99:99:99,mem=8gb -q js -j oe\n";
    print $jobfh "hostname\n";
    print $jobfh "module load stajichlab\nmodule load java\n";
    print $jobfh "cd $bdir\n";
    if( $force ) {
	unlink("$bdir/$base.intervals");
    }
    print $jobfh "if [ ! -f $base.$targetsp.bam ]; then\n" unless $force;
    print $jobfh "  java -Djava.io.tmpdir=$tmpdir $javamem -jar $picard/AddOrReplaceReadGroups.jar INPUT=$base.$targetsp.sam OUTPUT=$base.$targetsp.bam RGID=$base RGLB=$base RGPL=$platform RGPU=Genomic RGSM=$strain RGCN=$ctr SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=1000000\n";
    print $jobfh "fi\n" unless $force;

    print $jobfh "if [ ! -f $base.sort_ordered.$targetsp.bam ]; then\n" unless $force;
    print $jobfh "  java -Djava.io.tmpdir=$tmpdir $javamem -jar $picard/ReorderSam.jar INPUT=$base.$targetsp.bam OUTPUT=$base.sort_ordered.$targetsp.bam REFERENCE=$reference_genome MAX_RECORDS_IN_RAM=1000000 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT;\n";
    print $jobfh "fi\n" unless $force;

    print $jobfh "if [ ! -f $base.dedup.$targetsp.bam ]; then\n" unless $force;
    print $jobfh "  java -Djava.io.tmpdir=$tmpdir $javamem -jar $picard/MarkDuplicates.jar INPUT=$base.sort_ordered.$targetsp.bam OUTPUT=$base.dedup.$targetsp.bam METRICS_FILE=$base.dedup.$targetsp.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT;\n";
    print $jobfh "fi\n" unless $force;

    print $jobfh "if [ ! -f $base.$targetsp.intervals ]; then\n" unless $force;
    print $jobfh "  java -Djava.io.tmpdir=$tmpdir $javamem -jar $gatk_jarfile -T RealignerTargetCreator -R $reference_genome -o $base.$targetsp.intervals -I $base.dedup.$targetsp.bam;\n";
    print $jobfh "fi\n" unless $force;

    print $jobfh "if [ ! -f $base.realign.$targetsp.bam ]; then\n"unless $force;
    print $jobfh "  java -Djava.io.tmpdir=$tmpdir $javamem -jar $gatk_jarfile -T IndelRealigner -R $reference_genome -targetIntervals $base.$targetsp.intervals -I $base.dedup.$targetsp.bam -o $base.realign.$targetsp.bam;\n";
    print $jobfh "fi\n" unless $force;
    
    #print $jobfh "if [ ! -f $base.snps.raw.vcf ]; then\n" unless $force;
    #print $jobfh "  java $tmpdir $javamem -jar $gatk_jarfile -nt $threads -T UnifiedGenotyper -R $reference_genome -I $base.realign.bam -o $base.snps.raw.vcf -glm BOTH -rf BadCigar ;\n";
    #print $jobfh "fi\n" unless $force;
}
