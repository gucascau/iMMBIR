#!/usr/bin/perl
#author:wangxin
### function: The script is to extract the MAT sequences and calcuated the reads counts

use strict;
use warnings;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","m:s","n:s","o:s","p:s","q:s","r:s","f:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i} || !defined $opts{m} || !defined $opts{n}|| !defined $opts{o}|| !defined $opts{p} || !defined $opts{q} || !defined $opts{r} || !defined $opts{f}) {
       	die "************************************************************************
       	Usage: $0
			-i: work directory
			-m: index of first sample
			-n: index of second sample
			-p: index of third sample
			-q: index of fourth sample			
			-f: index of fifth sample
			
			-r: query MATa search string
			-o: Output string of MAT with coverage
************************************************************************\n";
}

my $dir=$opts{i};
my $input1=$opts{m};
my $input2=$opts{n};
my $input3=$opts{p};
my $input4=$opts{q};
my $input5=$opts{f};
my $output=$opts{o};

my $qstring=$opts{r};

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

my %stringD; my $id4;
open IN4, "$dir/$input4/Detection/$input4\_merged.assembled.fasta" or die $!;
while (<IN4>){
	chomp;
	if (/>(\S+)/){
		$id4=$1;
	}else {
		$stringD{$id4}.=substr $_,3, -3;
	}
}
close IN4;


my %stringE; my $id5;
open IN5, "$dir/$input5/Detection/$input5\_merged.assembled.fasta" or die $!;
while (<IN5>){
	chomp;
	if (/>(\S+)/){
		$id5=$1;
	}else {
		$stringE{$id4}.=substr $_,3, -3;
	}
}
close IN5;


my %hash;
foreach my $i (keys %stringA){
	next unless ($stringA{$i}=~/$qstring/);
	my $inf1=$stringA{$i};
	$hash{$inf1}->{A}++;
	
}

foreach my $i (keys %stringB){
	next unless ($stringB{$i}=~/$qstring/);
	my $inf2=$stringB{$i};
	$hash{$inf2}->{B}++;	
}


foreach my $i (keys %stringC){
	next unless ($stringC{$i}=~/$qstring/);
	my $inf3=$stringC{$i};
	$hash{$inf3}->{C}++;
	
}

foreach my $i (keys %stringD){
	next unless ($stringD{$i}=~/$qstring/);
	my $inf4=$stringD{$i};
	$hash{$inf4}->{D}++;	
}

foreach my $i (keys %stringE){
	next unless ($stringE{$i}=~/$qstring/);
	my $inf5=$stringE{$i};
	$hash{$inf5}->{E}++;	
}


open OUT, ">$output.MAT.txt" or die $!;
print OUT "MATsequence\tType\tTotalCov\tCov.$input1\tCov.$input2\tCov.$input3\tCov.$input4\tCov.$input5\n";
foreach my $n (keys %hash){
	
	my $numA=(exists $hash{$n}->{A})?$hash{$n}->{A}:0;
	my $numB=(exists $hash{$n}->{B})?$hash{$n}->{B}:0;
	
	my $numC=(exists $hash{$n}->{C})?$hash{$n}->{C}:0;
	my $numD=(exists $hash{$n}->{D})?$hash{$n}->{D}:0;
	my $numE=(exists $hash{$n}->{E})?$hash{$n}->{E}:0;
	
	my $num=$numA+$numB+$numC+$numD+$numE;
	print OUT "$n\t$output\t$num\t$numA\t$numB\t$numC\t$numD\t$numE\n";
}

close OUT;


