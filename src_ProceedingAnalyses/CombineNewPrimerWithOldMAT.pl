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
if (!defined $opts{i} || !defined $opts{o}  || !defined $opts{r} || defined $opts{h} || !defined $opts{m} || !defined $opts{n}  ) {
	die "************************************************************************
	Usage: $0.pl -i MATa insertion short primer -r MATa insertion long primer  -o MATa insertion short primer for the genotype
	
	Request Parameters:
	-i the list of MATa insertion from short primer
	-r the list of MATa insertion from longer primer
	-o output of MATa insertions for the genotye
	
	-m the number of samples for short primer (default 2)
	-n the number of samples for longer primer (default 2)
	-h Help
************************************************************************\n";
}


### Read the short primer MATa
my $SMat=$opts{i};

### Read the long primer MATa
my $LMat=$opts{r};

## Output folder for the MATa

my $output=$opts{o};

### the number of samples in each categories
my $NS=(defined $opts{m})?$opts{m}:2;
my $NL=(defined $opts{n})?$opts{n}:2;

 my %hash; my $finaltype;
open S, "$SMat" or die $!;
while (<S>){
	chomp;
	s/\r//g;
	my @array=(split/\t/,$_);
	next if ($array[0] eq "ID");
	my $Sn=4+$NS;
	my ($strain,$index,$Scov)=@array[1,2,3];
	my ($MatchPos,$InsSize,$Seq)=@array[$Sn..($Sn+2)];
	my $string=join "\t",@array[4..($Sn-1)];
	$hash{$index}->{Scov}=$Scov;
	$hash{$index}->{SSeq}=$Seq;
	$hash{$index}->{SSize}=$InsSize;
	$hash{$index}->{SMatch}=$MatchPos;
	$hash{$index}->{SStr}=$string;
	$finaltype=$strain;
}

close S;

open L,"$LMat" or die $!;
while (<L>){
	chomp;
	s/\r//g;
	my @array1=(split/\t/,$_);
	next if ($array1[0] eq "ID");
	my $Ln=4+$NL;
	my ($strain,$index,$Scov)=@array1[1,2,3];
	my ($MatchPos,$InsSize,$Seq0,$Seq)=@array1[$Ln..($Ln+3)];
	my $string=join "\t",@array1[4..($Ln-1)];
	$hash{$index}->{Lcov}=$Scov;
	$hash{$index}->{LSeq}=$Seq;
	$hash{$index}->{LSize}=$InsSize;
	$hash{$index}->{LMatch}=$MatchPos;
	$hash{$index}->{LStr}=$string;
}

close L;

open OUT,">$output" or die $!;

my $FtStr;
for my $num (1..($NS+$NL)){
	$FtStr.="\tCov$num";
}

print OUT "ID\tGenotype\tInsertionFlank\tTotalCov$FtStr\tMatchPosition\tInsSize\tSequence\n";

my $n=0;
foreach my $i (keys %hash){
	$n++;
	my $finalCov; my $finalStr; my $finalMatch; my $finalSize; my $finalSeq;
	if (exists $hash{$i}->{Scov} && exists $hash{$i}->{Lcov}){
		
		$finalCov=$hash{$i}->{Scov} +$hash{$i}->{Lcov};
		$finalStr= join "\t", ($hash{$i}->{SStr},$hash{$i}->{LStr});
		$finalSeq=$hash{$i}->{SSeq};
		$finalMatch=$hash{$i}->{SMatch};
		$finalSize=$hash{$i}->{SSize};
		
	}elsif (!exists $hash{$i}->{Scov} && exists $hash{$i}->{Lcov}){
		$finalCov=$hash{$i}->{Lcov};
		$finalStr= join "\t", (("0")x$NS,$hash{$i}->{LStr});
		$finalSeq=$hash{$i}->{LSeq};
		$finalMatch=$hash{$i}->{LMatch};
		$finalSize=$hash{$i}->{LSize};
		
	}else{
		$finalCov=$hash{$i}->{Scov};
		$finalStr= join "\t", ($hash{$i}->{SStr},("0")x$NL,);
		$finalSeq=$hash{$i}->{SSeq};
		$finalMatch=$hash{$i}->{SMatch};
		$finalSize=$hash{$i}->{SSize};
	
	}
	
	print OUT "$finaltype:$n\t$finaltype\t$i\t$finalCov\t$finalStr\t$finalMatch\t$finalSize\t$finalSeq\n";
	
}

close OUT;