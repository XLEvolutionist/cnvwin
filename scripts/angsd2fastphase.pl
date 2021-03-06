#!usr/bin/perl
use strict;
use warnings;
use List::MoreUtils qw(indexes);

# usage: script.pl <genotype.angsd.output> <ancestral.fa> > output

# if the file is stil gziped then the inout file should be:
# <(zcat genotype.angsd.output | awk '{$1 == 1 print}' )
# to go through chromosome one..

				##################
####################################################
# Convert a ANGSD genotype file to fastPHASE input #
####################################################
				##################

##																						##
# This script works on single chr information, thus data need to be provided on a		 #
# per chromosomes basis. The scrip WILL NOT WORK otherwise.								 #
#																						 #
# using the code:'awk '$1==1 {print}' infile > out.file' to select chr one, for example  #
##																						##

# define some useful variables
my %top;
my %bottom;
my @pos;

# open a file to which bed file info will be exported
open (BED , ">$ARGV[2].temp.bed") || die "could not open file $ARGV[2].temp.bed:$!\n";

# open the input file

open ( IN , "<$ARGV[0]" ) || die "could not open file $ARGV[0]:$!\n";

# cycle through the file
while ( <IN> ) {
	chomp;
	my ( $chr , $pos , $major , $minor , @samples ) = split /\s+/;
	my $bedPos = $pos-1;
	print BED "$chr\t$bedPos\t$pos\n";
	foreach my $sampleNum ( 0 .. $#samples ) {
		my @snps = split( // , $samples[$sampleNum] );
		push(@{$top{$sampleNum}} , $snps[0] );
		push(@{$bottom{$sampleNum}} , $snps[1] );
	}# foreach
	push(@pos,$pos);
}#while

###############################################
# Make a fasta file of all the variable sites #
###############################################

system("bedtools getfasta -fi $ARGV[1] -bed $ARGV[2].temp.bed -fo $ARGV[2].temp.fasta");

########################################################################
# load in that data into an array, this contains the ancestral alleles #
########################################################################

open ( INFASTA , "<$ARGV[2].temp.fasta" ) || die "could not open file $ARGV[2].temp.fasta:$!\n";

# slurp the file in
my @fasta = <INFASTA>;

# remove the useless header lines
chomp @fasta;
@fasta = grep(!/>/ , @fasta );

# now find the indices of the "N" nucleotides, so the can be removed later
my @x = indexes { !/N/i } @fasta;
@pos = @pos[@x];

#open up an output file to record the included positions
open(POS , ">$ARGV[2].pos" ) || die "Could not open file >$ARGV[2].pos:$!\n";
#record the final positions and create a SNP file for `rehh`

for my $i ( 0 .. $#pos ) {
	print POS "$i\t$ARGV[2]\t$pos[$i]\t0\t1\n";
}#foreach

# cleanse  @fasta info of Ns
@fasta = grep(!/N/i , @fasta);

my @keys =  keys ( %top );
@keys= sort { $a<=>$b } @keys; 

###########################
# Print the output header #
###########################

# the number of individuals
print scalar @keys , "\n";
# the number of segregating sites
print scalar @fasta , "\n";

# fastPHASE cannot take more than 1000000 on each line, so we have to break up the lines
my $n=499999;

########################################################
# Loop through each samples, modify and print the data #
########################################################

foreach my $key ( @keys ) {
    #remove nucleotides where TRIP is "N" based on the index of matches
	@{$top{$key}} = (@{$top{$key}})[@x];
	@{$bottom{$key}} = (@{$bottom{$key}})[@x];
    #convert nuc values to 0 and 1 according to if they have a match to TRIP.
    #0 is ancestral 1 is derived, so 0 if match, 1 if no match
	for my $i ( 0 .. $#fasta ) {
		
		#########################################################################
		# Convert the data to 0s and 1s according to matches to ancestral state #
		#########################################################################
		if ( ${$top{$key}}[$i] eq "N"  ) { ${$top{$key}}[$i]= "?"}# if
		# if the snp matches the ancestral (i.e. fasta) then make the value 0
		elsif ( uc($fasta[$i]) eq ${$top{$key}}[$i]  ) { ${$top{$key}}[$i]= 0}# if
		# if it does not macth, it is derived and should have value 1
		else { ${$top{$key}}[$i]= 1 }
		
		# dito here
		if ( ${$bottom{$key}}[$i] eq "N"  ) { ${$bottom{$key}}[$i]= "?"}# if
		elsif ( uc($fasta[$i]) eq ${$bottom{$key}}[$i]  ) { ${$bottom{$key}}[$i]= 0}# if
		else { ${$bottom{$key}}[$i]= 1 }
		 
	}# for
	#break up the line according to the max number of characters allowed per line
	print "# id ", $key , "\n";
	while ( my @x = splice ( @{$top{$key}} , 0 , $n ) )  {
	 	print @x ,"\n";
	}#while
	while ( my @x = splice ( @{$bottom{$key}} , 0 , $n ) )  {
	 	print @x ,"\n";
	}#while
	#print @{$top{$key}} , "\n";
	#print @{$bottom{$key}} , "\n";
}#for

exit;

###############
# Subroutines #
###############

