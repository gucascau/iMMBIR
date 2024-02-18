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
SampleID=DNA2

SampleID1=yWH475-1_S10
SampleID2=yWH475-2_S11
SampleID3=yYY423-1_S12
SampleID4=yYY423-2_S13
SampleID5=yWH475_S1
SampleID6=yWH475-colony1_S2




#in=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/Updated_V0903/UsedMutants/WT/Finished_blast

in=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/Updated_V0903/UsedMutants/DNA2/Finished_blast
softpath=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/Updated_V0903/UsedMutants/TestSoftWare/Software


Chr3M=${softpath}/iMMBIR/Database/Chr3.masked.fasta

genomeseq=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa

cd ${in}


mkdir ${SampleID}CombinedMATaInsertionSix
cd ${in}/${SampleID}CombinedMATaInsertionSix

perl ${softpath}/iMMBIR/src/MAT_detection_for_six_samples_v2.pl -i ${in}/ -m ${SampleID1} -n ${SampleID2}  -p ${SampleID3}  -q ${SampleID4} -f ${SampleID5} -g  ${SampleID6} -o ${SampleID}

perl -ne '{chomp; my @array=split/\t/,$_; next if ($array[0] eq "MATsequence"); $n++;my $id=join ":",@array[1..8],$n; print ">$id\n$array[0]\n"}' ${SampleID}.MAT.txt >${SampleID}.MAT.fasta

#bedtools maskfasta -fi Chr3.fasta -bed Masked.bed -fo Chr3.masked.fasta
#bwa index Chr3.masked.fasta
bwa mem -A2 -E1 ${Chr3M} ${SampleID}.MAT.fasta >${SampleID}.MAT.sam

perl ${softpath}/iMMBIR/src/MAT_detection_for_deduplication_six.pl  -i ${SampleID}.MAT.sam -g ${SampleID}.MAT.fasta -o ${SampleID}.uniq -t ${SampleID}



