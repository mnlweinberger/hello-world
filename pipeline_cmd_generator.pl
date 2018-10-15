use strict;

if (scalar(@ARGV) < 1) {die "USAGE: pipeline_cmd_generator INPUT_FILE\n";}
my $inFileName = $ARGV[0];
my $outFileName = $ARGV[1];
open IN, "<$inFileName" or die "can't open '$inFileName'";
open OUT, ">$outFileName" or die "can't open '$outFileName'";
my ($line, $iodine, $untreated);


foreach $line (<IN>){	
	chomp $line;
	my (@split, $split);
	@split = split(" ", $line);
	$iodine = $split[0];
	$untreated = $split[1];
#	print $iodine."\t";
#	print $untreated."\n";
	print OUT "set -x\n";

# First step - check the alignment rates to make sure everything is going OK at this point
	print OUT "echo check alignment rate\n";
	print OUT  "echo $iodine.sam\n";
	print OUT  "samtools view -S -c -F 4 $iodine.sam\n";
	print OUT  "samtools view -S -c -f 4 $iodine.sam\n";
	print OUT  "echo $untreated.sam\n";
	print OUT  "samtools view -S -c -F 4 $untreated.sam\n";
	print OUT  "samtools view -S -c -f 4 $untreated.sam\n";

	print OUT "\n\n";

# Now create a sepeate header and sam body files for each file, so there will be no problems in the next steps of adding UMIs and collapsing PCR duplicates
	print OUT  "echo prepare sam and header file\n";
	print OUT  "samtools view -S $iodine.sam > noheader_$iodine.sam\n";
	print OUT  "samtools view -S -H $iodine.sam > header_$iodine.txt\n";
	print OUT  "samtools view -S $untreated.sam > noheader_$untreated.sam\n";
	print OUT  "samtools view -S -H $untreated.sam > header_$untreated.txt\n";

	print OUT "\n\n";

# Next, run add_smart.pl that is the UMI adding script, adding a field to the sam file containing the UMI
	print OUT "echo The next step is adding UMIs to the sam file, with ulimit -d 33,000,000 perl collapse_good.pl\n";
	print OUT "ulimit -d 33000000\n";
	print OUT "perl add_smart.pl F_UMIs_$iodine.txt noheader_$iodine.sam UMIs_noheader_$iodine.sam\n";
	print OUT "ulimit -d 33000000\n";
	print OUT "perl add_smart.pl F_UMIs_$untreated.txt noheader_$untreated.sam UMIs_noheader_$untreated.sam\n";
	print OUT "echo finished adding UMIs\n";
	
	print OUT "echo Now adding headers\n";
	print OUT "cat header_$iodine.txt UMIs_noheader_$iodine.sam > F_UMIs_$iodine.sam\n";
	print OUT "cat header_$untreated.txt UMIs_noheader_$untreated.sam > F_UMIs_$untreated.sam\n";

	print OUT "\n\n";

# Now that the sam file is ready, I want to collapse it
	print OUT "echo Started collapsing for IP files\n";
	print OUT "ulimit -d 33000000\n";
	print OUT "perl -w collapse_smart.pl F_UMIs_$iodine.sam Collapsed_$iodine.sam\n";
	print OUT "ulimit -d 33000000\n";
	print OUT "perl -w collapse_smart.pl F_UMIs_$untreated.sam Collapsed_$untreated.sam\n";

# Samtools is now used to turn SAM into BAM, then sort and index files
	print OUT "samtools view -bS Collapsed_$iodine.sam 1>Collapsed_$iodine.bam 2>Collapsed_$iodine.view\n";
	print OUT "samtools view -bS Collapsed_$untreated.sam 1>Collapsed_$untreated.bam 2>Collapsed_$untreated.view\n";

	print OUT "samtools sort Collapsed_$iodine.bam Collapsed_$iodine\_sorted\n";
	print OUT "samtools sort Collapsed_$untreated.bam Collapsed_$untreated\_sorted\n";

	print OUT "samtools index Collapsed_$iodine\_sorted.bam\n";
	print OUT "samtools index Collapsed_$untreated\_sorted.bam\n";

	print OUT "\n\n";

# This step is time consuming as the pileup needs to calculate depth of reads and mutations in each position of each chromosome
	print OUT "echo Now, continue to clalculate mutations\n";
	print OUT "samtools  mpileup -d 2500000 -f /nadata/data/genomes/human/hg19/bowtie2/hg19.fa -A Collapsed_$iodine\_sorted.bam > Collapsed_$iodine.pileup\n";
	print OUT "samtools  mpileup -d 2500000 -f /nadata/data/genomes/human/hg19/bowtie2/hg19.fa -A Collapsed_$untreated\_sorted.bam > Collapsed_$untreated.pileup\n";

	print OUT "\n\n";

	print OUT "cat Collapsed_$iodine.pileup";
	print OUT q{| awk -F "\t" '{if ($4>5) {split($5,a,"");i=1;while(i<=length(a)){start=i;if(a[i]=="." || a[i]==",") {C["Ref"]++;i++} if(a[i]~/[acgtACGT]/) {C[toupper(a[i])]++;i++}; if(a[i]=="-" || a[i]=="+") {c=a[i];num=a[i+1];j=i+2;while(a[j]~/[[:digit:]]/) {num=num*10+a[j];j++} c=c num;for(k=0;k<num;k++) c=c a[j+k];C[c]++;i=j+k} if(a[i]=="^") i+=2;if(a[i]=="<" || a[i]==">") C["skip"]++;if(start==i) i++} num=0;for(val in C) if(val!~/skip/) num+=C[val]; if(num-C["Ref"]>0) {printf "%s\t%i\t%s\t%i\t",$1,$2,$3,num;for(val in C) {printf "%s:%i;",val,C[val]; printf "\t";}printf "\n";fflush(); } delete C;} }' > };
	print OUT "parsed_$iodine\n";
	
	print OUT "cat Collapsed_$untreated.pileup";
	print OUT q{| awk -F "\t" '{if ($4>5) {split($5,a,"");i=1;while(i<=length(a)){start=i;if(a[i]=="." || a[i]==",") {C["Ref"]++;i++} if(a[i]~/[acgtACGT]/) {C[toupper(a[i])]++;i++}; if(a[i]=="-" || a[i]=="+") {c=a[i];num=a[i+1];j=i+2;while(a[j]~/[[:digit:]]/) {num=num*10+a[j];j++} c=c num;for(k=0;k<num;k++) c=c a[j+k];C[c]++;i=j+k} if(a[i]=="^") i+=2;if(a[i]=="<" || a[i]==">") C["skip"]++;if(start==i) i++} num=0;for(val in C) if(val!~/skip/) num+=C[val]; if(num-C["Ref"]>=0) {printf "%s\t%i\t%s\t%i\t",$1,$2,$3,num;for(val in C) {printf "%s:%i;",val,C[val]; printf "\t";}printf "\n";fflush(); } delete C;} }' > };
	print OUT "parsed_$untreated\n";

	print OUT "\n\n";

# Fixing indels to reflect correct position
	print OUT "cat parsed_$iodine";
	print OUT q{| awk '{indel=0;for(i=5;i<=NF;i++) {if($i~/Ref/) {split($i,a,";");split(a[1],b,":");} if($i~/[-+]/) {split($i,a,";");split(a[1],c,":");indel+=c[2]}} print $1"_"$2+1,indel}' > };
	print OUT "parsed_$iodine"."1\n";
	print OUT "cat parsed_$iodine";
	print OUT q{|awk '{not_ref=0;for(i=5;i<=NF;i++) {if($i~/Ref/) {split($i,a,";");split(a[1],b,":");} if($i!~/[-+]/) {split($i,a,";");split(a[1],c,":");not_ref+=c[2]}} print $1"_"$2,$3,b[2],not_ref}' > };
	print OUT "parsed_$iodine"."2\n";

	print OUT "cat parsed_$untreated";
	print OUT q{| awk '{indel=0;for(i=5;i<=NF;i++) {if($i~/Ref/) {split($i,a,";");split(a[1],b,":");} if($i~/[-+]/) {split($i,a,";");split(a[1],c,":");indel+=c[2]}} print $1"_"$2+1,indel}' > };
	print OUT "parsed_$untreated"."1\n";
	print OUT "cat parsed_$untreated";
	print OUT q{|awk '{not_ref=0;for(i=5;i<=NF;i++) {if($i~/Ref/) {split($i,a,";");split(a[1],b,":");} if($i!~/[-+]/) {split($i,a,";");split(a[1],c,":");not_ref+=c[2]}} print $1"_"$2,$3,b[2],not_ref}' > };
	print OUT "parsed_$untreated"."2\n";

	print OUT "\n\n";
	
	print OUT "echo finished calculating mutations and indels\n";
	print OUT "echo putting the two files together for accurte calculations\n";

# Iodine treated	
	print OUT "echo putting the two files together for accurte calculations\n";
	print OUT q{awk '{if(NR==FNR) {C[$1]=$2} else {OFS="\t";if($1 in C) {print $1,$2,$3,C[$1]+$NF} else {$(NF+1)=1;NF--;print}}}'};
	print OUT " parsed_$iodine"."1 parsed_$iodine"."2 ";
	print OUT q{| awk '{if($4>0){print $1"_"$2,1-($3/$4)} else {print $1"_"$2,0}}' > };
	print OUT "$iodine"."_mut_rate\n";

# Untreated
	print OUT "echo putting the two files together for accurte calculations\n";
	print OUT q{awk '{if(NR==FNR) {C[$1]=$2} else {OFS="\t";if($1 in C) {print $1,$2,$3,C[$1]+$NF} else {$(NF+1)=1;NF--;print}}}'};
	print OUT " parsed_$untreated"."1 parsed_$untreated"."2 ";
	print OUT q{| awk '{if($4>0){print $1"_"$2,1-($3/$4)} else {print $1"_"$2,0}}' > };
	print OUT "$untreated"."_mut_rate\n";

#Now create the iodine-untreated pairs for running cutoff critiria
 	print OUT "echo finished putting the two files together for accurte calculations\n";
	print OUT "echo making comparisons in mutation rates\n";
	print OUT q{awk '{if(NR==FNR) {C[$1]=$2} else {OFS="\t";if($1 in C) print $1,C[$1],$2}}'};
	print OUT " $iodine"."_mut_rate "."$untreated"."_mut_rate > $iodine"."-"."$untreated\n"; 

# Clean out the possible problem of upper-case/lower-case
	print OUT "echo the pairs are ready for intersection\n";
	print OUT "cat $iodine"."-"."$untreated";
	print OUT q{ | awk '{split($1,a,"_"); print a[1],a[2],a[3],$2,$3}' | awk '{if(toupper($3)=="A" || toupper($3)=="T") print $0}' > };
	print OUT "$iodine"."-"."$untreated"."_AT\n";

	print OUT "\n\n";
# Now it's time for running the critiria : FC>2, mut_rate>0.1 and diff>0.5
	print OUT "echo Running Xushen's critiria\n";
	print OUT q{awk 'OFS="\t" {a=10; if($5==0 && $4>0.1 && $4-$5>0.05) print $1,$2,$3,$4/a,$4-$5}{if($5!=0 && $4>0.1 && $4/$5>2 && $4-$5>0.05) print $1,$2,$3,$4/$5,$4-$5}' };
	print OUT "$iodine"."-"."$untreated"."_AT > "."$iodine"."-"."$untreated"."_AT_good_points.txt\n";

#Turn files to bed files to allow intersection
	print OUT q{awk 'OFS="\t" {print $1,$2-1,$2,"*","*",$3}'};
	print OUT "$iodine"."-"."$untreated"."_AT_good_points.txt > "."$iodine"."-"."$untreated"."_AT_good_points.bed\n";

# Running the intersection command
	print OUT "echo Now intersecting with known databases of small RNAs\n";
	print OUT "bedtools intersect \-a $iodine"."-"."$untreated"."_AT_good_points.bed \-b \/nadata\/users\/mor\/DASHR.bed \-wao";
	print OUT q{ | awk '{if($7~/chr/) print $0}' > };
	print OUT "$iodine"."-"."$untreated"."iodine_untreated_DASHR\n";
	print OUT "bedtools intersect \-a $iodine"."-"."$untreated"."_AT_good_points.bed \-b \/nadata\/users\/mor\/hg19-tRNAs.bed \-wao";
	print OUT q{ | awk '{if($7~/chr/) print $0}' > };
	print OUT "$iodine"."-"."$untreated"."iodine_untreated_tRNAs\n";

	print OUT "\n\n";

	print OUT "\# The script for $iodine and $untreated is ready\!\n";
	print "The script for $iodine and $untreated is ready\!\n";

}

