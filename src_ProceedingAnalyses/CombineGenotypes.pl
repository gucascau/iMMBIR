#!/usr/bin/perl
use strict;
use warnings;
#use Data::Dump qw(dump);
#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

#use Bio::Seq;
#use Bio::SeqIO;


#### the script is to combine MAT insertion with the latest new primer for each genotype
my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","r:s","h:s","m:s","n:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o} || !defined $opts{r}  ) {
	die "************************************************************************
	Usage: $0.pl -i folder contain lists of MATa insertions  -o Output insertions with index and read counts information across genotype
	
	Request Parameters:
	-i the folder that contain MATa insertions (WT-DSB, Sgs1, Sgs1EXO1, Exo1, Rad27, and DNA2)
	-r Genotype
	-o Output MATa insertions with index and read counts information across genotype
	

************************************************************************\n";
}


### Read the folder of indels. we put the files with name_indel.deletion.stat or name_indel.insertion.stat
my $indels=$opts{i};
my $output=$opts{o};
my $genotype=$opts{r};
my @file;

opendir MT,"$indels" or die $!;
 @file=readdir MT;
closedir MT;
my $n=0; my %hash; my %num; my %seq; my %size; my %matchP;
foreach my $i (@file){
	next unless ($i=~/.txt/);
	
	print "$i\n";

	$n++;
	#my $type=join "_",(split/\_/,$i)[0,1];

	#### read clinical file in each folder: and then combine the mutaiton, clinical restults. 	
	open FILE,"$indels/$i" or die $!;	
	while (<FILE>){
		chomp;
		s/\r//;
		s/^\s+//;
		my ($type,$event,$cov,$matchpos,$Isize,$sequence)=(split/\t/,$_)[1,2,3,4,5,6];
		next if ($type eq "Genotype" || $type eq "StrainType" );
		#$hash{$event}->{$type}=$cov/$totalcounts{$type};
		$num{$event}->{$type}=$cov;
		$seq{$event}=$sequence;
		$hash{$type}++;
		$size{$event}=$Isize;
		$matchP{$event}=$matchpos;
	}
	close FILE;

}


open OUT,">$output" or die $!;
print OUT "ID\tGenotype\tInsertionFlank\t";
foreach my $j (sort keys %hash){
	print OUT "\t$j";
}
print OUT "\tTotalCov\tMatchPosition\tInsSize\tSequence\n";

my $IDnum=0;
foreach my $i (keys %num){
	$IDnum++;
	print OUT "$genotype:$IDnum\t$genotype\t$i";
	
	my $Tnum=0;
	foreach my $j (sort keys %hash){
		my $fnum=(exists $num{$i}->{$j})?($num{$i}->{$j}):0;
		
		print OUT "\t$fnum";
		$Tnum += $fnum;
	}
	print OUT "\t$Tnum\t$matchP{$i}\t$size{$i}\t$seq{$i}\n";
}

close OUT;