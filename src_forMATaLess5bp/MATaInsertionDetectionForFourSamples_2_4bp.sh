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
SampleID=Sgs1Exo1_3_4BP

SampleID1=yGI199-A1_S6
SampleID2=yGI199-A2_S7
SampleID3=yGI199-B1_S8
SampleID4=yGI199-B2_S9


in=/project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/Sgs1exo1


softpath=/scratch/ch220812/Software



Chr3M=${softpath}/iMMBIR/Database/Chr3.masked.fasta

genomeseq=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa

cd ${in}


mkdir ${SampleID}CombinedMATaInsertion
cd ${in}/${SampleID}CombinedMATaInsertion

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/Sgs1exo1/MAT_detection_for_four_samples_3_4bp.pl -i ${in}/ -m ${SampleID1} -n ${SampleID2}  -p ${SampleID3}  -q ${SampleID4}  -o ${SampleID}

perl -ne '{chomp; my @array=split/\t/,$_; next if ($array[0] eq "MATsequence"); $n++;my $id=join ":",@array[1..6],$n; print ">$id\n$array[0]\n"}' ${SampleID}.MAT.txt >${SampleID}.MAT.fasta

#bedtools maskfasta -fi Chr3.fasta -bed Masked.bed -fo Chr3.masked.fasta
#bwa index Chr3.masked.fasta
bwa mem -A2 -E1 ${Chr3M} ${SampleID}.MAT.fasta >${SampleID}.MAT.sam

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/Sgs1exo1/MAT_detection_for_deduplication_four_2bp_v3.pl  -i ${SampleID}.MAT.sam -g ${SampleID}.MAT.fasta -o ${SampleID}.uniq -t ${SampleID}


SampleID2B=Sgs1Exo1_2BP

### The following is for 2bp

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/Sgs1exo1/MAT_detection_for_four_samples_2bp.pl -i ${in}/ -m ${SampleID1} -n ${SampleID2}  -p ${SampleID3}  -q ${SampleID4}  -o ${SampleID2B}

perl -ne '{chomp; my @array=split/\t/,$_; next if ($array[0] eq "MATsequence"); $n++;my $id=join ":",@array[1..6],$n; print ">$id\n$array[0]\n"}' ${SampleID2B}.MAT.txt >${SampleID2B}.MAT.fasta

#bedtools maskfasta -fi Chr3.fasta -bed Masked.bed -fo Chr3.masked.fasta
#bwa index Chr3.masked.fasta
bwa mem -A2 -E1 ${Chr3M} ${SampleID2B}.MAT.fasta >${SampleID2B}.MAT.sam

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/Sgs1exo1/MAT_detection_for_deduplication_four_2bp_v3.pl -i ${SampleID2B}.MAT.sam -g ${SampleID2B}.MAT.fasta -o ${SampleID2B}.uniq -t ${SampleID2B}


