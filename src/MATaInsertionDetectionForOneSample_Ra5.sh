#!/bin/bash
# Sample batchscript to run a simple job array to create 5 different files, and modify the files independently with one batchscript
 
#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=1:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=yYY521_indel # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails
#SBATCH --mail-user=xin.wang@childrens.harvard.edu # Email address to send the job status
#SBATCH --output=output_%A_%a.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=20 # Number of cpu cores on one node
#SBATCH --mem=20G



# pipeline to detect the large insertion events of MATa
source /programs/biogrids.shrc

SampleID=yYY517A_S10

Samplename=yYY517A_S10

in=/project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/20220812_ChroAging/Update_V210924


softpath=/scratch/ch220812/Software


Chr3M=${softpath}/iMMBIR/Database/Chr3.masked.fasta

genomeseq=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa

cd ${in}


mkdir ${SampleID}CombinedMATaInsertion
cd ${in}/${SampleID}CombinedMATaInsertion

perl ${softpath}/iMMBIR/src/MAT_detection_for_one_samples_v2.pl -i ${in}/${SampleID}/Detection -m ${SampleID} -o ${SampleID}

perl -ne '{chomp; my @array=split/\t/,$_; next if ($array[0] eq "MATsequence"); $n++;my $id=join ":",@array[1..2],$n; print ">$id\n$array[0]\n"}' ${SampleID}.MAT.txt >${SampleID}.MAT.fasta

#bedtools maskfasta -fi Chr3.fasta -bed Masked.bed -fo Chr3.masked.fasta
#bwa index Chr3.masked.fasta
bwa mem -A2 -E1 ${Chr3M} ${SampleID}.MAT.fasta >${SampleID}.MAT.sam

perl ${softpath}/iMMBIR/src/MAT_detection_for_deduplication_one.pl -i ${SampleID}.MAT.sam -g ${SampleID}.MAT.fasta -o ${Samplename}.uniq -t ${SampleID}



