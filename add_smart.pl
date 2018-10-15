use strict;

#####################################################################################################
# This script takes a line with a read header and a UMI comprised of bith R1 and R2 UMIs
# and creates a hash where the header is the key and the UMIs are the value.
# Then, it opens the SAM file, and looks for a the same header as the keys.
# If a keu exists that is in the hash, it prints out the sam line + the value of this hash
####################################################################################################

# Open infile and read the lines:
if (scalar(@ARGV) < 1) {die "USAGE: NO_INPUT_FILE\n";}
my $inFileName1 = $ARGV[0];
my $inFileName2 = $ARGV[1];
my $outFileName = $ARGV[2];
open IN1, "<$inFileName1" or die "can't open '$inFileName1'"; # UMIs file
open IN2, "<$inFileName2" or die "can't open '$inFileName2'"; # SAM file
open OUT, ">$outFileName" or die "can't open '$outFileName'"; # out SAM file

my (%hash, $hash, $umis);
my ($key,$value);

# create a hash that contains the header of the read as the key, and the UMI in it's SAM-friendly format as a value:
foreach $umis (<IN1>){
	chomp $umis;
	my @split = split("\t",$umis);
	my $key = $split[0];
	$value = $split[1];
	$hash{$key} = $value;
}
#while( my( $key, $value ) = each %hash ){
#    print "$key \t $value\n";
#}
print "Done reading UMIs\n";

#Now, open the SAM file:
my $sam_header;
my $line;

#foreach SAM line, search for the header in SAM, and check to see if the header appears as a key (the header from the UMIs file) in the hash I just created.
#If so, print the line (SAM line) and the value (UMI) of this hash
foreach $line (<IN2>){
	chomp $line;
	if ($line =~ m/(HISEQ\:\d+\:\w+\d+\w+\:\d+\:\d+\d+\:\d+\:\d+)/){
		$sam_header = $1;
#		print $sam_header;
		if (exists($hash{$sam_header})){
			print OUT $line."\t"."Xm:Z:".$hash{$sam_header}."\n";
			} else {
				print "There are headers that have not been added UMIs"."\n";
			}
		}
}
close IN1;
close IN2;
close OUT;
