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
SampleID=WTD16MATB

SampleID1=yYY517A-M-16d_S19
SampleID2=yYY517B-M-16d_S20
#SampleID3=JKM139-3_S4
#SampleID4=JKM139-4_S5
String=GAAGTTTTATA

in=/project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/20210712_ChroAging/Update_V210924


softpath=/scratch/ch220812/Software



Chr3M=${softpath}/iMMBIR/Database/Chr3.masked.fasta

genomeseq=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa

cd ${in}


mkdir ${SampleID}CombinedMATaInsertion
cd ${in}/${SampleID}CombinedMATaInsertion

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/20210712_ChroAging/Update_V210924/OtherMATaIns/scripts/MAT_detection_for_two_samples_SearchIndex.pl -i ${in}/ -m ${SampleID1} -n ${SampleID2}   -o ${SampleID} -q ${String}

#perl ${softpath}/iMMBIR/src/MAT_detection_for_four_samples_v2.pl -i ${in}/ -m ${SampleID1} -n ${SampleID2}  -p ${SampleID3}  -q ${SampleID4}  -o ${SampleID}

perl -ne '{chomp; my @array=split/\t/,$_; next if ($array[0] eq "MATsequence"); $n++;my $id=join ":",@array[1..4],$n; print ">$id\n$array[0]\n"}' ${SampleID}.MAT.txt >${SampleID}.MAT.fasta

#bedtools maskfasta -fi Chr3.fasta -bed Masked.bed -fo Chr3.masked.fasta
#bwa index Chr3.masked.fasta
bwa mem -A2 -E1 ${Chr3M} ${SampleID}.MAT.fasta >${SampleID}.MAT.sam

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/20210712_ChroAging/Update_V210924/OtherMATaIns/scripts/MAT_detection_for_deduplication_two_SearchIndex.pl -i ${SampleID}.MAT.sam -g ${SampleID}.MAT.fasta -o ${SampleID}.uniq -t ${SampleID}


### 
SampleID2B=Pol32_2BP


### The following is for 2bp

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/Pol32/MAT_detection_for_two_samples_2bp.pl -i ${in}/ -m ${SampleID1} -n ${SampleID2}   -o ${SampleID2B}


perl -ne '{chomp; my @array=split/\t/,$_; next if ($array[0] eq "MATsequence"); $n++;my $id=join ":",@array[1..4],$n; print ">$id\n$array[0]\n"}' ${SampleID2B}.MAT.txt >${SampleID2B}.MAT.fasta

#bedtools maskfasta -fi Chr3.fasta -bed Masked.bed -fo Chr3.masked.fasta
#bwa index Chr3.masked.fasta
bwa mem -A2 -E1 ${Chr3M} ${SampleID2B}.MAT.fasta >${SampleID2B}.MAT.sam

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/Pol32/MAT_detection_for_deduplication_two_2bp_v3.pl -i ${SampleID2B}.MAT.sam -g ${SampleID2B}.MAT.fasta -o ${SampleID2B}.uniq -t ${SampleID2B}


cat ${SampleID}.uniq.final.txt ${SampleID2B}.uniq.final.txt >Pol32_2_4BP.uniq.final.txt


