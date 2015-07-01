#!/usr/bin/perl
use strict;
use warnings;

# Simon Renny-Byfield, UC Davis, March 17th 2015

#usage script.pl <angsd output> >outfile

# This script will take angsd output (.pestPG file) and convert format suitable 
# for input into GenomicRanges. The purpose is to associate regions in maize genome
# that have estimates of Tajima's D and relate those estimates to particular genes. 
# In short each window from the angsd analysis will be labelled with the corresponding 
# gene name using GenomicRanges. 

# The output from angsd looks like this..

# ## thetaStat VERSION: 0.01 build:(Oct 23 2014,15:39:05)
# #(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)	Chr	WinCenter	tW	tP	tF	tH	tL	Tajima	fuf	fud	fayh	zeng	nSites
# (25,53)(1778322,1778350)(1778250,1778350)	1	1778300	1.661604	0.550991	3.218805	1.968841	1.259916	-2.122087	-1.758895	-1.220118	-0.736941	-0.208055	28
# (25,96)(1778322,1778532)(1778300,1778400)	1	1778350	3.663915	2.025260	7.192596	3.483926	2.754594	-1.622298	-1.918688	-1.580332	-0.370203	-0.240505	71
# (53,96)(1778350,1778532)(1778350,1778450)	1	1778400	2.002310	1.474270	3.973791	1.515086	1.494677	-0.869529	-1.449146	-1.362891	-0.017994	-0.225741	43
# (96,127)(1778532,1778731)(1778500,1778600)	1	1778550	0.235098	0.097436	0.000001	0.005128	0.051282	-0.980203	0.205032	0.575739	0.213366	-0.367884	31
# (114,127)(1778550,1778731)(1778550,1778650)	1	1778600	0.235098	0.097436	0.000001	0.005128	0.051282	-0.980203	0.205032	0.575739	0.213366	-0.367884	13
# (127,146)(1778731,1778750)(1778650,1778750)	1	1778700	0.000504	0.000107	0.002142	0.000003	0.000055	-0.066506	-0.096080	-0.089562	0.006211	-0.021466	19

#but a suitable input into GenomicRanges might look like this

# ## thetaStat VERSION: 0.01 build:(Oct 23 2014,15:39:05)
# #(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)	Chr	WinCenter	tW	tP	tF	tH	tL	Tajima	fuf	fud	fayh	zeng	nSites
# 1 1778250	1778350	1778300	1.661604	0.550991	3.218805	1.968841	1.259916	-2.122087	-1.758895	-1.220118	-0.736941	-0.208055	28

#where the first three columns are chr, start, stop.

####
# Read in the data
####

# open the input file

open ( IN , "<$ARGV[0]" ) || die "could not open file:$!\n";

####
# Cycle through the file modifying each line
####

while ( <IN> )  {
	if ( m/##/ ) {
		print;
		next
	}
	if ( m/#\(/ ) {
		print "chr	start	stop	WinCenter	tW	tP	tF	tH	tL	Tajima	fuf	fud	fayh	zeng	nSites\n";
		next
	}#if
	chomp;
	my @data = split "\t";
	# grab the start and stop of each window
	$data[0] =~ m/\(\d+,\d+\)\(\d+,\d+\)\((\d+,\d+)\)/;
	my $new = $1;
	#replace a comma with tab
	$new =~ s/,/\t/;
	#print the new line as you want it.
	print $data[1] , "\t" , $new , "\t" , join( "\t" , @data[2 .. $#data]) , "\n";
}#while

exit;
