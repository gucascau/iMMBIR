#!/usr/bin/perl
use strict;
use warnings;


#### This scripts is to perform the MAT deduplications after their detection


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","t:s","g:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o} || !defined $opts{t}   || !defined $opts{g}  ) {
	die "************************************************************************
	Usage: $0.pl 
				-i Mapping results
				-g Input fasta file
				-t Strain type
				-o Insertion reads and deletion reads
************************************************************************\n";
}

my $output=$opts{o};
my $input=$opts{i};
my $Stype=$opts{t};
my $fasta=$opts{g};



#### input of sam file ####


my %cov; my %partial; my %seq; my %hash; my %Isize;

my $idS; my %sequence;
open F,"$fasta" or die "cannot open file $fasta";

while (<F>){
	chomp;
	if (/>(\S+)/){
		$idS=$1;
	}else{
		$sequence{$idS}.=$_;
	}
}

close F;


open I,"$input" or die "cannot open file $input";

while (<I>){
	chomp;
	s/\r//g;
	next if (/^\@/);

	my ($i,$qual,$chr,$start,$match,$string)=(split /\t/,$_)[0,1,2,3,5,9];
	my $length=length $string;
	
	if ($qual>0 && $match !~m/^(\d+)S(\d+)M(\d+)S$/){
		
		print "$_\n";
		next;
	}
	next if (exists $hash{$i});
	
	$cov{$i}= (split/:/,$i)[1];
	
	$hash{$i}++;
	#### the majority of the duplications for 2bp is the deletion events: 33M6D45M, or 26M13D45M
	
	### no other indel; right have one deletion; right have two different deletions 
	
	if($match =~m/^(\d+)M(\d+)D(\d+)M$/){
		
		### normally left side are at least 25bp
		if($1>=25 && $1<41){
			my $part1=$1-5; 
			my $part2=$2; my $part3=$1;
	
			my $string1=substr($string,$part1,5);
			my $string2=substr($string,$part3,5);
		
			my $str=join '',($string1,"D",$part2,"D",$string2);
			$partial{$i}=$str;
			$Isize{$i}=0;
			$seq{$i}=$_;
		}else{
			print "$_\n";
		}
		### here are the main insertion events
	}elsif($match =~m/^(\d+)M(\d+)I(\d+)M$/){
		
		my $part1=$1-5; my $part2=$1; my $part3=$1+$2;
		
		my $string1=substr($string,$part1,5);
		my $string3=substr($string,$part3,5);
		my $string2=substr($string, $part2,$2);
			
		my $str=join '',($string1,$2,"I",$string2,"I",$string3);
		
		$Isize{$i}=$2;
		
		$partial{$i}=$str;
		$seq{$i}=$_;
	
		### there are deletion or insertion caused by sequence errors. Here we need to determine where they come from
	}elsif($match =~m/^(\d+)M(\d+)(I|D)(\d+)M(\d+)(I|D)(\d+)M$/){
		
		####certian sequences errors in the right part
		if ($1>=25 && $1<41 ){
			
			# the insertion or deletion events depend on $3
			## insertion event
			if ($3 eq "I"){
				my $part1=$1-5; my $part2=$1; my $part3=$1+$2;
		
				my $string1=substr($string,$part1,5);
				my $string3=substr($string,$part3,5);
				my $string2=substr($string, $part2,$2);
			
				my $str=join '',($string1,$2,"I",$string2,"I",$string3);
		
				$Isize{$i}=$2;
		
				$partial{$i}=$str;
				$seq{$i}=$_;
			## deletion event	
			}else{
				my $part1=$1-5; 
				my $part2=$2; my $part3=$1;
	
				my $string1=substr($string,$part1,5);
				my $string2=substr($string,$part3,5);
		
				my $str=join '',($string1,"D",$part2,"D",$string2);
				$partial{$i}=$str;
				$Isize{$i}=0;
				$seq{$i}=$_;
			}
			
		}elsif(($1+$4)>=25 && ($1+$4)<41 ){
			
			## insertion event
			if($6 eq "I"){
				
				my $part1=$1+$4-5; my $part2=$1+$4; my $part3=$1+$4+$5;
		
				my $string1=substr($string,$part1,5);
				my $string3=substr($string,$part3,5);
				my $string2=substr($string, $part2,$5);
			
				my $str=join '',($string1,$5,"I",$string2,"I",$string3);
		
				$Isize{$i}=$5;
		
				$partial{$i}=$str;
				$seq{$i}=$_;
			
			### deletion event	
			}else{
				my $part1=$1+$4-5; 
				my $part2=$5; my $part3=$1+$4;
	
				my $string1=substr($string,$part1,5);
				my $string2=substr($string,$part3,5);
		
				my $str=join '',($string1,"D",$part2,"D",$string2);
				$partial{$i}=$str;
				$Isize{$i}=0;
				$seq{$i}=$_;
			}
		}
	}else{
		print "$_\n";
	}	
}

close I;



my %Tcov; my %MaxCov; my %MaxId;my %TcovA;my %TcovB;my %TcovC;my %TcovD;

foreach my $m (keys %partial){
	my $barcode=$partial{$m};
	
	my $id=(split/\t/,$seq{$m})[0];
	
	my ($covA,$covB)=(split/:/,$id)[2,3];

	my $coverage = $cov{$id};
	$Tcov{$barcode}+=$coverage;
	
	$TcovA{$barcode}+=$covA;
	$TcovB{$barcode}+=$covB;
	# $TcovC{$barcode}+=$covC;
# 	$TcovD{$barcode}+=$covD;
	
	
	if (not defined $MaxCov{$barcode}){
		$MaxCov{$barcode}=$coverage;
		$MaxId{$barcode}=$id;
	}
	if($coverage >= $MaxCov{$barcode}){
		$MaxCov{$barcode} =$coverage;
		$MaxId{$barcode}=$id;
	}
	
	
}

open OUT,">$output" or die $!;
open OUTF,">$output.final.txt" or die $!;

print OUTF "ID\tStrainType\tInsertionFlank\tTotalCov\tCovA\tCovB\tMatchPosition\tInsSize\tSequence\n";
my $num=0;
foreach my $n (sort {$Tcov{$b} <=>$Tcov{$a}} keys %Tcov){
	
	my $name=$MaxId{$n};
	
	my $IMAT=(split/\t/,$seq{$name})[5];
	
	my $ISeq=$sequence{$name};
	
	
	
	$num++;
	my $Sid=join ":",($Stype,$num);
	my $Ssize=$Isize{$name};
	print OUTF "$Sid\t$Stype\t$n\t$Tcov{$n}\t$TcovA{$n}\t$TcovB{$n}\t$IMAT\t$Ssize\t$ISeq\n";
	print OUT "$name\t$n\t$Tcov{$n}\t$TcovA{$n}\t$TcovB{$n}\t$Ssize\t$ISeq\t$seq{$name}\n";
}

close OUT;
close OUTF;

