#!usr/bin/perl
use strict;
use warnings;

# usage: script.pl <genotype.angsd.output> >output

#####
# Convert a ANGSD genotype file to fastphase input
#####

#
my %top;
my %bottom;
# open the input file

open ( IN , "<$ARGV[0]" ) || die "could not open file:$!\n";

while ( <IN> ) {
	chomp;
	my ( $chr , $pos , $major , $minor , @samples ) = split /\s+/;
	foreach my $sampleNum ( 0 .. $#samples ) {
		my @snps = split( // , $samples[$sampleNum] );
		push(@{$top{$sampleNum}} , $snps[0] );
		push(@{$bottom{$sampleNum}} , $snps[1] );
	}#
}#while

my @keys =  keys ( %top );
@keys= sort(@keys); 


######
# Print the output data
######

# print the header info
print scalar @keys , "\n";

# now print the genotype info
foreach my $key ( @keys ) {
	print "# ", $key , "\n";
	print @{$top{$key}} , "\n";
	print @{$bottom{$key}} , "\n";
}#for

