#!usr/bin/perl
use strict;
use warnings;

# usage: script.pl <mafs.gz>

#####
# Convert a ANGSD mafs.gz file to SweepFinder format.
#####

#declare some useful variables
my $previousChr=0;

# first open the maf output from ANGSD

open( MAFS,"<$ARGV[0]" ) || die "Could not open file:$!\n";

# print a needed header

# scan through the input file
while ( <MAFS> ) {
	next if m/position/;
	chomp;
	# a variable to hold the frequency
	my $freq;
	# assume the SNP is folded
	my $fold=0;
	# split the incoming line;
	my ( $chromo, $position, $major, $minor, $ref, $anc, $knownEM, $nInd ) =
			split /\t/;
	
	# see if we need to open up a new file
	if ( $previousChr ne $chromo ) {
		close OUT;
		# open up the initial output file
		my $file = "$chromo"."_sweepfinder_input.txt";
		open( OUT, ">$file" );
		print OUT "position\tx\tn\tfolded\n";
	}# if
	$previousChr = $chromo;		
	# if we have a ancestral state declare the SNP to be unfolded.
	if ( $anc eq "N" ) { $fold=1; }
	# next if the frequency of the minor allele is zero.
	next if $knownEM == 0;
	# if the minor is equal to the anc, then we need to flip the frequency
	# because the major allele is now the derived.
	if ( $minor eq $anc ) {
		$freq = (2*$nInd)-($knownEM*($nInd*2));
	}# if
	else {
		$freq=$knownEM*(2*$nInd);
	}# else
	next if $freq == 0;
	print OUT $position , "\t"; 
	printf OUT ("%.0f",$freq);
	print OUT "\t";
	print OUT 2*$nInd , "\t" , $fold , "\n";
}#while

exit;