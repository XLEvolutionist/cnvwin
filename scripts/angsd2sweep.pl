#!usr/bin/perl
use strict;
use warnings;

# usage: script.pl <mafs.gz>

#####
# Convert a ANGSD mafs.gz file to SweepFinder format.
#####

# first open the maf output from ANGSD

open( MAFS,"<$ARGV[0]" ) || die "Could not open file:$!\n";

# print a needed header
print "position\tx\tn\tfolded\n";

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
	# if we have a ancestral state declare the SNP to be unfolded.
	if ( $anc eq "N" ) { $fold=1; }
	# next if the frequency of the minor allele is zero.
	next if $knownEM == 0;
	# if the minor is equal to the anc, then we need to flip the frequency
	# because the major allele is now the derived.
	if ( $minor eq $anc ) {
		$freq = $nInd-($knownEM*$nInd);
	}# if
	else {
		$freq=$knownEM*$nInd;
	}# else
	print $position , "\t" , $freq, "\t" , $nInd , "\t" , $fold , "\n";
}#while

exit;