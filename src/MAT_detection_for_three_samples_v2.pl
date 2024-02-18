#!/usr/bin/perl
#author:wangxin
### function: The script is to extract the MAT sequences and calcuated the reads counts

use strict;
use warnings;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","m:s","n:s","o:s","p:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i} || !defined $opts{m} || !defined $opts{n}|| !defined $opts{o}|| !defined $opts{p} ) {
       	die "************************************************************************
       	Usage: $0
			-i: work directory
			-m: index of first sample
			-n: index of second sample
			-p: index of third sample
			-o: Output string of MAT with coverage
************************************************************************\n";
}

my $dir=$opts{i};
my $input1=$opts{m};
my $input2=$opts{n};
my $input3=$opts{p};
my $output=$opts{o};

my %stringA; my $id1;
open IN1, "$dir/$input1/Detection/$input1\_merged.assembled.fasta" or die $!;
while (<IN1>){
	chomp;
	if (/>(\S+)/){
		$id1=$1;
	}else {
		$stringA{$id1}.=substr $_,3, -3;
	}
}

close IN1;

my %stringB; my $id2;
open IN2, "$dir/$input2/Detection/$input2\_merged.assembled.fasta" or die $!;
while (<IN2>){
	chomp;
	if (/>(\S+)/){
		$id2=$1;
	}else {
		$stringB{$id2}.=substr $_,3, -3;
	}
}
close IN2;

my %stringC; my $id3;
open IN1, "$dir/$input3/Detection/$input3\_merged.assembled.fasta" or die $!;
while (<IN1>){
	chomp;
	if (/>(\S+)/){
		$id3=$1;
	}else {
		$stringC{$id3}.=substr $_,3, -3;
	}
}

close IN1;



my %hash;
foreach my $i (keys %stringA){
	next unless ($stringA{$i}=~/TACTTCAGTATA/|| $stringA{$i}=~/TACTTCAGCATA/);
	my $inf1=$stringA{$i};
	$hash{$inf1}->{A}++;
	
}

foreach my $i (keys %stringB){
	next unless ($stringB{$i}=~/TACTTCAGTATA/|| $stringB{$i}=~/TACTTCAGCATA/);
	my $inf2=$stringB{$i};
	$hash{$inf2}->{B}++;	
}


foreach my $i (keys %stringC){
	next unless ($stringC{$i}=~/TACTTCAGTATA/|| $stringC{$i}=~/TACTTCAGCATA/);
	my $inf3=$stringC{$i};
	$hash{$inf3}->{C}++;
	
}



open OUT, ">$output.MAT.txt" or die $!;
print OUT "MATsequence\tType\tTotalCov\tCov.$input1\tCov.$input2\tCov.$input3\n";
foreach my $n (keys %hash){
	
	my $numA=(exists $hash{$n}->{A})?$hash{$n}->{A}:0;
	my $numB=(exists $hash{$n}->{B})?$hash{$n}->{B}:0;
	
	my $numC=(exists $hash{$n}->{C})?$hash{$n}->{C}:0;

	
	my $num=$numA+$numB+$numC;
	print OUT "$n\t$output\t$num\t$numA\t$numB\t$numC\n";
}

close OUT;


