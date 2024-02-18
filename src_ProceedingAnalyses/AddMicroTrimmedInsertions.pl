#!/usr/bin/perl
#author:wangxin
### function: The script is to extract the MAT sequences and add the microhomology and trimmed size

use strict;
use warnings;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","m:s","o:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i} || !defined $opts{m} ||  !defined $opts{o}) {
       	die "************************************************************************
       	Usage: $0
			-i: the requested MATa insertion events
			-m: index of MATa insertions with microhomolgy
			-o: Output string of MAT with microhimology
************************************************************************\n";
}

my $input=$opts{i};
my $index=$opts{m};
my $output=$opts{o};





my %stringA; 
open IND, "$index" or die $!;
while (<IND>){
	chomp;
	s/\r//g;
	my $indexSeq=(split/\t/,$_)[0];
	$stringA{$indexSeq}=$_;
}

close IND;

open IN, "$input" or die $!;
open OUT,">$output" or die $!;
while (<IN>){
	chomp;
	my $indexquery=(split/\t/,$_)[0];
	s/\r//g;
	if ($indexquery eq "Index"){
		print OUT "$_\tIndex\tSequence\tInsertionSize\tMATaDeletion\tMicrohomologySize\tMatL\tInsMatR\n";
	}else{
		
		my $string=(exists $stringA{$indexquery})?$stringA{$indexquery}:"0\t0\t0\t0\t0\t0\t0";
		
		print OUT "$_\t$string\n";
	}
	
}

close IN;

close OUT;


