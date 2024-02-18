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
	##### Matches to identify MAT insertions, we ignore the indel mismatches in the MAT-L  and MAT-R #####
	## right part have several bp deletions or no indels only one insertions or Right MAT side have two different deletions	
	
	### no other indel; right have one deletion; right have two different deletions 
	if($match =~m/^(\d+)M(\d+)I(\d+)M$/ || $match =~m/^(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M$/ || $match =~m/^(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M$/ ||$match =~m/^(\d+)M(\d+)I(\d+)M(\d+)S$/ || $match =~m/^(\d+)S(\d+)M(\d+)S$/){
		
		my $part1=$1-5; my $part2=$1; my $part3=$1+$2;
		
		my $string1=substr($string,$part1,5);
		my $string3=substr($string,$part3,5);
		my $string2=substr($string, $part2,$2);
			
		my $str=join '',($string1,$2,"I",$string2,"I",$string3);
		
		$Isize{$i}=$2;
		
		$partial{$i}=$str;
		$seq{$i}=$_;
		
	### Left  part have one deletions; Left part have one deletion and right also have one deletion

	}elsif($match =~m/^(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M$/ ||$match =~m/^(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M$/  ){
		
		my $part1=$1+$3-5; my $part2=$1+$3; my $part3=$1+$3+$4;
		
		my $string1=substr($string,$part1,5);
		my $string3=substr($string,$part3,5);
		my $string2=substr($string, $part2,$4);
		
		my $str=join '',($string1,$4,"I",$string2,"I",$string3);
		
		$partial{$i}=$str;
		$Isize{$i}=$4;
		$seq{$i}=$_;
		
	
	### Left or right part have one insertion, or  right have one extral deletion 
	
	}elsif($match =~m/^(\d+)M(\d+)I(\d+)M(\d+)I(\d+)M$/ || $match =~m/^(\d+)M(\d+)I(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M$/ ){
		
		if(($1+$3)>35 && ($1+$3)<45){
			
			my $part1=$1+$2+$3-5;
			my $part2=$1+$2+$3;
			my $part3=$1+$2+$3+$4;
			my $string1=substr($string,$part1,5);
			my $string3=substr($string,$part3,5);
			my $string2=substr($string, $part2,$4);
			
			my $str=join '',($string1,$4,"I",$string2,"I",$string3);
		
			$partial{$i}=$str;
			$Isize{$i}=$4;
			$seq{$i}=$_;
			
		}else {	
			my $part1=$1-5;
			my $part2=$1;
			my $part3=$1+$2;
			
			
			my $string1=substr($string,$part1,5);
			my $string3=substr($string,$part3,5);
			my $string2=substr($string, $part2,$2);
			
			my $str=join '',($string1,$2,"I",$string2,"I",$string3);
		
			$partial{$i}=$str;
			$Isize{$i}=$2;
			$seq{$i}=$_;
		}
	
	
	### Left MAT side have two different deletions
	}elsif ($match =~m/^(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M$/){
		my $part1=$1+$3+$5-5; my $part2=$1+$3+$5;
		my $part3=$1+$3+$5+$6;
		
		my $string1=substr($string,$part1,5);
		my $string3=substr($string,$part3,5);
		my $string2=substr($string, $part2,$6);
		
		my $str=join '',($string1,$6,"I",$string2,"I",$string3);
		
		$partial{$i}=$str;
		$Isize{$i}=$6;
		$seq{$i}=$_;
	
	
	
		### 
	}elsif($match =~m/^(\d+)M(\d+)I(\d+)M(\d+)I(\d+)M(\d+)I(\d+)M$/){
		
		### Right MAT three insertions
		if ($1>35 && $1<45){
			my $part1=$1-5; my $part2=$1;
			my $part3=$1+$2;
			my $string1=substr($string,$part1,5);
			my $string3=substr($string,$part3,5);
			my $string2=substr($string, $part2,$2);
		
			my $str=join '',($string1,$2,"I",$string2,"I",$string3);
		
			$partial{$i}=$str;
			$Isize{$i}=$2;
			$seq{$i}=$_;	
			
		### Left MAT have one insertion , right have two insertions	
		}elsif (($1+$2+$3)>35 && ($1+$2+$3)<45){
			my $part1=$1+$2+$3-5; my $part2=$1+$2+$3;
			my $part3=$1+$2+$3+$4;
		
			my $string1=substr($string,$part1,5);
			my $string3=substr($string,$part3,5);
			my $string2=substr($string, $part2,$4);
		
			my $str=join '',($string1,$4,"I",$string2,"I",$string3);
		
			$partial{$i}=$str;
			$Isize{$i}=$4;
			$seq{$i}=$_;	
			
		### Left have three insertions	
		}else{
			
			my $part1=$1+$2+$3+$4+$5-5; my $part2=$1+$2+$3+$4+$5;
			my $part3=$1+$2+$3+$4+$5+$6;
		
			my $string1=substr($string,$part1,5);
			my $string3=substr($string,$part3,5);
			my $string2=substr($string, $part2,$4);
		
			my $str=join '',($string1,$6,"I",$string2,"I",$string3);
		
			$partial{$i}=$str;
			$Isize{$i}=$6;
			$seq{$i}=$_;		
			
		}
	}elsif($match =~m/^(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M(\d+)I(\d+)M$/){
		
		if($1+$2+$3>35 && ($1+$3+$2)<45){
			my $part1=$1+$3-5; my $part2=$1+$3;
			my $part3=$1+$3+$4;
			
			my $string1=substr($string,$part1,5);
			my $string3=substr($string,$part3,5);
			my $string2=substr($string, $part2,$4);
		
			my $str=join '',($string1,$4,"I",$string2,"I",$string3);
		
			$partial{$i}=$str;
			$Isize{$i}=$4;
			$seq{$i}=$_;
			
		}else{
			my $part1=$1+$3+$4+$5-5; my $part2=$1+$3+$4+$5;
			my $part3=$1+$3+$4+$5+$6;
			
			my $string1=substr($string,$part1,5);
			my $string3=substr($string,$part3,5);
			my $string2=substr($string, $part2,$6);
			
			my $str=join '',($string1,$6,"I",$string2,"I",$string3);
		
			$partial{$i}=$str;
			$Isize{$i}=$6;
			$seq{$i}=$_;	
			
		}
		
	}elsif($match =~m/^(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M$/){
		
		if($1>35 && $1<45){
			my $part1=$1-5; my $part2=$1;
			my $part3=$1+$2;
			
			my $string1=substr($string,$part1,5);
			my $string3=substr($string,$part3,5);
			my $string2=substr($string, $part2,$2);
			
			my $str=join '',($string1,$2,"I",$string2,"I",$string3);
		
			$partial{$i}=$str;
			$Isize{$i}=$2;
			$seq{$i}=$_;	
		}else{
			my $part1=$1+$2+$3+$5-5; my $part2=$1+$2+$3+$5;
			my $part3=$1+$2+$3+$5+$6;
			
			my $string1=substr($string,$part1,5);
			my $string3=substr($string,$part3,5);
			my $string2=substr($string, $part2,$6);
			
			my $str=join '',($string1,$6,"I",$string2,"I",$string3);
		
			$partial{$i}=$str;
			$Isize{$i}=$6;
			$seq{$i}=$_;	
		}
	
	}elsif($match =~m/^(\d+)S(\d+)M(\d+)I(\d+)M$/){
		
		my $part1=$1+$2-5; my $part2=$1+$2;
		my $part3=$1+$2+$3;
		
		my $string1=substr($string,$part1,5);
		my $string3=substr($string,$part3,5);
		my $string2=substr($string, $part2,$3);
		
		my $str=join '',($string1,$3,"I",$string2,"I",$string3);
	
		$partial{$i}=$str;
		$Isize{$i}=$3;
		$seq{$i}=$_;	
		
		
		
	}else{
		print "$_\n";
	}
}

close ;



my %Tcov; my %MaxCov; my %MaxId;
#my %TcovA;my %TcovB;my %TcovC;my %TcovD;

foreach my $m (keys %partial){
	my $barcode=$partial{$m};
	
	my $id=(split/\t/,$seq{$m})[0];
	
	my ($covA,$covB)=(split/:/,$id)[2,3];

	my $coverage = $cov{$id};
	$Tcov{$barcode}+=$coverage;
	
	# $TcovA{$barcode}+=$covA;
# 	$TcovB{$barcode}+=$covB;
# 	# $TcovC{$barcode}+=$covC;
# # 	$TcovD{$barcode}+=$covD;
	
	
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

print OUTF "ID\tStrainType\tInsertionFlank\tTotalCov\tMatchPosition\tInsSize\tSequence\n";
my $num=0;
foreach my $n (sort {$Tcov{$b} <=>$Tcov{$a}} keys %Tcov){
	
	my $name=$MaxId{$n};
	
	my $IMAT=(split/\t/,$seq{$name})[5];
	
	my $ISeq=$sequence{$name};
	
	
	
	$num++;
	my $Sid=join ":",($Stype,$num);
	my $Ssize=$Isize{$name};
	print OUTF "$Sid\t$Stype\t$n\t$Tcov{$n}\t$IMAT\t$Ssize\t$ISeq\n";
	print OUT "$name\t$n\t$Tcov{$n}\t$Ssize\t$ISeq\t$seq{$name}\n";
}

close OUT;
close OUTF;

