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
SampleID=WTDSB_M3_4BP

SampleID1=JKM139-1_S2
SampleID2=JKM139-2_S3
SampleID3=JKM139-3_S4
SampleID4=JKM139-4_S5
SampleID5=JKM139_S1
SampleID6=WT-notreatment-2h_S7




#in=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/Updated_V0903/UsedMutants/WT/Finished_blast

in=/project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/WT

softpath=/scratch/ch220812/Software


Chr3M=${softpath}/iMMBIR/Database/Chr3.masked.fasta

genomeseq=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa

cd ${in}


mkdir ${SampleID}CombinedMATaInsertionSix
cd ${in}/${SampleID}CombinedMATaInsertionSix

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/WT/MAT_detection_for_six_samples_3_4bp.pl -i ${in}/ -m ${SampleID1} -n ${SampleID2}  -p ${SampleID3}  -q ${SampleID4} -f ${SampleID5} -g  ${SampleID6} -o ${SampleID}

perl -ne '{chomp; my @array=split/\t/,$_; next if ($array[0] eq "MATsequence"); $n++;my $id=join ":",@array[1..8],$n; print ">$id\n$array[0]\n"}' ${SampleID}.MAT.txt >${SampleID}.MAT.fasta

#bedtools maskfasta -fi Chr3.fasta -bed Masked.bed -fo Chr3.masked.fasta
#bwa index Chr3.masked.fasta
bwa mem -A2 -E1 ${Chr3M} ${SampleID}.MAT.fasta >${SampleID}.MAT.sam

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/WT/MAT_detection_for_deduplication_six_2bp_v3.pl   -i ${SampleID}.MAT.sam -g ${SampleID}.MAT.fasta -o ${SampleID}.uniq -t ${SampleID}



### For the 2bp potential MATa insertions
SampleID2BP=WTDSB_2BP

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/WT/MAT_detection_for_six_samples_2bp.pl -i ${in}/ -m ${SampleID1} -n ${SampleID2}  -p ${SampleID3}  -q ${SampleID4} -f ${SampleID5} -g  ${SampleID6} -o ${SampleID2BP}

perl -ne '{chomp; my @array=split/\t/,$_; next if ($array[0] eq "MATsequence"); $n++;my $id=join ":",@array[1..8],$n; print ">$id\n$array[0]\n"}'  ${SampleID2BP}.MAT.txt > ${SampleID2BP}.MAT.fasta

#bedtools maskfasta -fi Chr3.fasta -bed Masked.bed -fo Chr3.masked.fasta
#bwa index Chr3.masked.fasta
bwa mem -A2 -E1 ${Chr3M} ${SampleID2BP}.MAT.fasta > ${SampleID2BP}.MAT.sam

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/WT/MAT_detection_for_deduplication_six_2bp_v3.pl  -i  ${SampleID2BP}.MAT.sam -g ${SampleID2BP}.MAT.fasta -o ${SampleID2BP}.uniq -t ${SampleID2BP}
