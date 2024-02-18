#!/usr/bin/perl
#author:wangxin
### function: The script is to extract the MAT sequences and calcuated the reads counts

use strict;
use warnings;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","m:s","n:s","o:s","p:s","q:s","f:s","g:s","a:s","b:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i} || !defined $opts{m} || !defined $opts{n}|| !defined $opts{o}|| !defined $opts{p} || !defined $opts{q}|| !defined $opts{f} || !defined $opts{g}|| !defined $opts{a} || !defined $opts{b}) {
       	die "************************************************************************
       	Usage: $0
			-i: work directory
			-m: index of first sample
			-n: index of second sample
			-p: index of third sample
			-q: index of fourth sample
			-f: index of fifth sample
			-g: index of six sample
			-a: index of seven sample
			-b: index of eight sample
			-o: Output string of MAT with coverage
************************************************************************\n";
}

my $dir=$opts{i};
my $input1=$opts{m};
my $input2=$opts{n};
my $input3=$opts{p};
my $input4=$opts{q};
my $input5=$opts{f};
my $input6=$opts{g};
my $input7=$opts{a};
my $input8=$opts{b};
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
		$stringE{$id5}.=substr $_,3, -3;
	}
}
close IN5;


my %stringF; my $id6;
open IN6, "$dir/$input6/Detection/$input6\_merged.assembled.fasta" or die $!;
while (<IN6>){
	chomp;
	if (/>(\S+)/){
		$id6=$1;
	}else {
		$stringF{$id6}.=substr $_,3, -3;
	}
}
close IN6;



my %stringG; my $id7;
open IN7, "$dir/$input7/Detection/$input7\_merged.assembled.fasta" or die $!;
while (<IN7>){
	chomp;
	if (/>(\S+)/){
		$id7=$1;
	}else {
		$stringG{$id7}.=substr $_,3, -3;
	}
}
close IN7;




my %stringH; my $id8;
open IN8, "$dir/$input8/Detection/$input8\_merged.assembled.fasta" or die $!;
while (<IN8>){
	chomp;
	if (/>(\S+)/){
		$id8=$1;
	}else {
		$stringH{$id8}.=substr $_,3, -3;
	}
}
close IN8;

my %hash;
foreach my $i (keys %stringA){
	next unless ((($stringA{$i}=~/TTCAGTATA/) || ($stringA{$i}=~/TTCAGCATA/)) && ($stringA{$i}!~/CTTCAGTATA/) && ($stringA{$i}!~/CTTCAGCATA/));
	my $inf1=$stringA{$i};
	$hash{$inf1}->{A}++;
	
}

foreach my $i (keys %stringB){
	next unless  ((($stringB{$i}=~/TTCAGTATA/) || ($stringB{$i}=~/TTCAGCATA/)) && ($stringB{$i}!~/CTTCAGTATA/) && ($stringB{$i}!~/CTTCAGCATA/));
	my $inf2=$stringB{$i};
	$hash{$inf2}->{B}++;	
}


foreach my $i (keys %stringC){
	next unless ((($stringC{$i}=~/TTCAGTATA/) || ($stringC{$i}=~/TTCAGCATA/)) && ($stringC{$i}!~/CTTCAGTATA/) && ($stringC{$i}!~/CTTCAGCATA/));
	my $inf3=$stringC{$i};
	$hash{$inf3}->{C}++;
	
}

foreach my $i (keys %stringD){
	next unless ((($stringD{$i}=~/TTCAGTATA/) || ($stringD{$i}=~/TTCAGCATA/)) && ($stringD{$i}!~/CTTCAGTATA/) && ($stringD{$i}!~/CTTCAGCATA/));

	my $inf4=$stringD{$i};
	$hash{$inf4}->{D}++;	
}

foreach my $i (keys %stringE){
	next unless ((($stringE{$i}=~/TTCAGTATA/) || ($stringE{$i}=~/TTCAGCATA/)) && ($stringE{$i}!~/CTTCAGTATA/) && ($stringE{$i}!~/CTTCAGCATA/));
	
	my $inf5=$stringE{$i};
	$hash{$inf5}->{E}++;
	
}

foreach my $i (keys %stringF){
	next unless ((($stringF{$i}=~/TTCAGTATA/) || ($stringF{$i}=~/TTCAGCATA/)) && ($stringF{$i}!~/CTTCAGTATA/) && ($stringF{$i}!~/CTTCAGCATA/));
	my $inf6=$stringF{$i};
	$hash{$inf6}->{F}++;	
}


foreach my $i (keys %stringG){
	next unless ((($stringG{$i}=~/TTCAGTATA/) || ($stringG{$i}=~/TTCAGCATA/)) && ($stringG{$i}!~/CTTCAGTATA/) && ($stringG{$i}!~/CTTCAGCATA/));
	my $inf7=$stringG{$i};
	$hash{$inf7}->{G}++;	
}

foreach my $i (keys %stringH){
	next unless ((($stringH{$i}=~/TTCAGTATA/) || ($stringH{$i}=~/TTCAGCATA/)) && ($stringH{$i}!~/CTTCAGTATA/) && ($stringH{$i}!~/CTTCAGCATA/));
	my $inf8=$stringH{$i};
	$hash{$inf8}->{H}++;	
}


open OUT, ">$output.MAT.txt" or die $!;
print OUT "MATsequence\tType\tTotalCov\tCov.$input1\tCov.$input2\tCov.$input3\tCov.$input4\tCov.$input5\tCov.$input6\tCov.$input7\tCov.$input8\n";
foreach my $n (keys %hash){
	
	my $numA=(exists $hash{$n}->{A})?$hash{$n}->{A}:0;
	my $numB=(exists $hash{$n}->{B})?$hash{$n}->{B}:0;
	
	my $numC=(exists $hash{$n}->{C})?$hash{$n}->{C}:0;
	my $numD=(exists $hash{$n}->{D})?$hash{$n}->{D}:0;
	
	my $numE=(exists $hash{$n}->{E})?$hash{$n}->{E}:0;
	my $numF=(exists $hash{$n}->{F})?$hash{$n}->{F}:0;
	
	my $numG=(exists $hash{$n}->{G})?$hash{$n}->{G}:0;
	my $numH=(exists $hash{$n}->{H})?$hash{$n}->{H}:0;
	
	
	my $num=$numA+$numB+$numC+$numD+$numE+$numF+$numG+$numH;
	print OUT "$n\t$output\t$num\t$numA\t$numB\t$numC\t$numD\t$numE\t$numF\t$numG\t$numH\n";
}

close OUT;


