#!/bin/bash
# Sample batchscript to run a simple job on HPC
#SBATCH --partition=mghpcc-compute              # queue to be used
#SBATCH --account=bch-mghpcc                    # account name to be used
#SBATCH --time=5:01:00                         # Running time (in hours-minutes-seconds)
#SBATCH --job-name=DDR              # Job name
#SBATCH --output=output_%j.txt                  # Name of the output file
#SBATCH --nodes=1                               # Number of nodes needed
#SBATCH --ntasks=15                             # Number of CPUs needed
#SBATCH --mem=20G                                # Memory needed for the computing session


### Define parameters by hands
source /programs/biogrids.shrc


SampleID=WT-DSB-LessFive

SampleID1=JKM139-1_S2
SampleID2=JKM139-2_S3
SampleID3=JKM139-3_S4
SampleID4=JKM139-4_S5


in=/project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/20190616_220605_FirstDNA2/Update_V210924/


softpath=/scratch/ch220812/Software


Chr3M=${softpath}/iMMBIR/Database/Chr3.masked.fasta

genomeseq=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa

cd ${in}


mkdir ${SampleID}CombinedMATaInsertion
cd ${in}/${SampleID}CombinedMATaInsertion

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/20190616_220605_FirstDNA2/Update_V210924/TestShortMATa/MAT_detection_for_two_samples_LessFive.pl -i ${in}  -m ${SampleID1} -n ${SampleID2}  -p ${SampleID3}  -q ${SampleID4} -o ${SampleID}


#perl ${softpath}/iMMBIR/src/MAT_detection_for_four_samples_v2.pl -i ${in}/ -m ${SampleID1} -n ${SampleID2}  -p ${SampleID3}  -q ${SampleID4}  -o ${SampleID}

perl -ne '{chomp; my @array=split/\t/,$_; next if ($array[0] eq "MATsequence"); $n++;my $id=join ":",@array[1..6],$n; print ">$id\n$array[0]\n"}' ${SampleID}.MAT.txt >${SampleID}.MAT.fasta
#bedtools maskfasta -fi Chr3.fasta -bed Masked.bed -fo Chr3.masked.fasta
#bwa index Chr3.masked.fasta
bwa mem -A2 -E1 ${Chr3M} ${SampleID}.MAT.fasta >${SampleID}.MAT.sam

perl ${softpath}/iMMBIR/src/MAT_detection_for_deduplication_four.pl  -i ${SampleID}.MAT.sam -g ${SampleID}.MAT.fasta -o ${SampleID}.uniq -t ${SampleID}


