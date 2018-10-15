use strict;

#####################################################################################################
# This script takes a SAM line, and creates a hash where the chr+start+CIGAR+UMI is saved as key
# the value being the entire line itself from the SAM file.
# It compares each new key to existing keys;
# If a key is new, print this line out to a new SAM file
# if a key exists, go to next line (skipping this one)
# The new file is in fact the collapsed SAM, as it only contains uniq keys, therefor unique reads
# that have CIGAR and UMIs that do not repeat themselves for the same genomic location.
####################################################################################################

# Open infile and read the lines:
if (scalar(@ARGV) < 1) {die "USAGE: class_ex.pl INPUT_FILE\n";}
my $inFileName = $ARGV[0];
my $outFileName = $ARGV[1];
open IN, "<$inFileName" or die "can't open '$inFileName'";
open OUT, ">$outFileName" or die "can't open '$outFileName'";

my $line;
my (%hash, $hash);
my ($key);
my @key_list;
my $counter = 0;
foreach $line (<IN>){	
	chomp $line;
	if ($line =~ m/^\@SQ/) {
		print OUT $line."\n";
# A typical line looks like this: HISEQ:188:CC1HKANXX:7:1101:1455:1998 #
	} else {if ($line =~ m/\w+\:\d+\:\w+\d+\w+\:\d+\:\d+\d+\:\d+\:\d+/){
		#print OUT $line."\n";
		$counter++;
		my @split = split("\t",$line);
		#my $split;
		# Now I create the key's name, which contains the chr+start+CIGAR+UMI 
		my $last = pop @split;
		$key = $split[2]."_".$split[3]."_".$split[5]."_".$last."\n";
} else {print "no hash exists".$line."\n"}
		
# Now check if there is already an identical key in the hash! if there is, I would like to skip this line.
# if this is a new key, I want to put it into a hash and print it out to a new SAM file
			if (!exists $hash{$key}){
			print OUT $line."\n";
			# create hash - the key is the chr+start+CIGAR+UMI, the value is the SAM line 		
			$hash{$key} = $line;
			} else {
				next
				}
		
	}
}

print "This is the number of reads checked $counter \n\n";







close IN;
close OUT;
