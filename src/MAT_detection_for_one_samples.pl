#!/usr/bin/perl
#author:wangxin
### function: The script is to extract the MAT sequences and calcuated the reads counts

use strict;
use warnings;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","m:s","o:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i} || !defined $opts{m} || !defined $opts{o}) {
       	die "************************************************************************
       	Usage: $0.pl
			-i: work directory
			-m: index of sample
			-o: Output string of MAT with coverage
************************************************************************\n";
}

my $dir=$opts{i};
my $input1=$opts{m};
#my $input2=$opts{n};
my $output=$opts{o};

my %stringA; my $id1;
open IN1, "$dir/$input1/$input1\_merge.assembled.fasta" or die $!;
while (<IN1>){
	chomp;
	if (/>(\S+)/){
		$id1=$1;
	}else {
		$stringA{$id1}.=substr $_,3, -3;
	}
}

close IN1;



my %hash;
foreach my $i (keys %stringA){
	next unless ($stringA{$i}=~/TACTTCAGTATA/|| $stringA{$i}=~/TACTTCAGCATA/);
	my $inf1=$stringA{$i};
	$hash{$inf1}->{A}++;
	
}


open OUT, ">$output.MAT.txt" or die $!;
print OUT "MATsequence\t$input1\n";
foreach my $n (keys %hash){
	
	my $numA=(exists $hash{$n}->{A})?$hash{$n}->{A}:0;
	#my $numB=(exists $hash{$n}->{B})?$hash{$n}->{B}:0;
	print OUT "$n\t$numA\n";
}

close OUT;


