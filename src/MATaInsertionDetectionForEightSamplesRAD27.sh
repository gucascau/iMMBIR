#!/bin/bash
# Sample batchscript to run a simple job on HPC
#SBATCH --partition=mghpcc-compute              # queue to be used
#SBATCH --account=bch-mghpcc                    # account name to be used
#SBATCH --time=40:01:00                         # Running time (in hours-minutes-seconds)
#SBATCH --job-name=DDR              # Job name
#SBATCH --output=output_%j.txt                  # Name of the output file
#SBATCH --nodes=1                               # Number of nodes needed
#SBATCH --ntasks=15                             # Number of CPUs needed
#SBATCH --mem=20G                                # Memory needed for the computing session


# pipeline to detect the large insertion events of MATa
source /programs/biogrids.shrc
SampleID=RAD27

SampleID1=yYY398-A_S9
SampleID2=yYY398-B_S10
SampleID3=yYY398B-30C_S15
SampleID4=yYY398C-30C_S16
SampleID5=yYY398D-30C_S17
SampleID6=yYY398B-23C_S18
SampleID7=yYY398C-23C_S19
SampleID8=yYY398D-23C_S20


#in=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/Updated_V0903/UsedMutants/WT/Finished_blast

in=/project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2
softpath=/scratch/ch220812/Software


Chr3M=${softpath}/iMMBIR/Database/Chr3.masked.fasta

genomeseq=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa


cd ${in}


mkdir ${SampleID}CombinedMATaInsertionEight
cd ${in}/${SampleID}CombinedMATaInsertionEight

perl ${softpath}/iMMBIR/src/MAT_detection_for_eight_samples_v2.pl -i ${in} -m ${SampleID1} -n ${SampleID2}  -p ${SampleID3}  -q ${SampleID4} -f ${SampleID5} -g  ${SampleID6} -a  ${SampleID7} -b ${SampleID8} -o ${SampleID}

perl -ne '{chomp; my @array=split/\t/,$_; next if ($array[0] eq "MATsequence"); $n++;my $id=join ":",@array[1..10],$n; print ">$id\n$array[0]\n"}' ${SampleID}.MAT.txt >${SampleID}.MAT.fasta

#bedtools maskfasta -fi Chr3.fasta -bed Masked.bed -fo Chr3.masked.fasta
#bwa index Chr3.masked.fasta
bwa mem -A2 -E1 ${Chr3M} ${SampleID}.MAT.fasta >${SampleID}.MAT.sam

perl ${softpath}/iMMBIR/src/MAT_detection_for_deduplication_eight.pl  -i ${SampleID}.MAT.sam -g ${SampleID}.MAT.fasta -o ${SampleID}.uniq -t ${SampleID}



