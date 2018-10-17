#########################################################
# This script does the first 2 steps of QC for PE runs:
#  
# 1. Takes R1 and R2 namesd in a file,
# 2. Runs FASTQC on each
# 3. Trims adaptors using trim_galore (PE mode)
# 4. Runs FASTQC on the trimmed files
#
#########################################################

# Open infile and read the lines:
if (scalar(@ARGV) < 1) {die "USAGE: class_ex.pl INPUT_FILE\n";}
my $inFileName = $ARGV[0];
open IN, "<$inFileName" or die "can't open '$inFileName'";
# infile shpuld include the fastq files in pairs, separated by " "
# For example: SC931_24a_HeLa_input_S8_R1_001 SC931_24a_HeLa_input_S8_R2_001

my $line;

foreach $line (<IN>){	
 	chomp $line;
	@split = split(" ", $line);
	$r1 = $split[0];
	$r2 = $split[1];
	system("/nadata/software/FastQC_0.11/fastqc $r1".".fastq.gz");
	system("/nadata/software/FastQC_0.11/fastqc $r2".".fastq.gz");
	system("/nadata/software/trim_galore --paired -q 20 --length 30 --trim-n  $r1 $r2");	
	system("/nadata/software/FastQC_0.11/fastqc $r1"."_val_1.fq.gz");
        system("/nadata/software/FastQC_0.11/fastqc $r2"."_val_2.fq.gz");
	}
