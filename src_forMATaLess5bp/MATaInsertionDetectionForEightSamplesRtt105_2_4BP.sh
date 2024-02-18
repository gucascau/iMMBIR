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
SampleID=Rtt105_3_4BP

SampleID1=yYY591A1_S15
SampleID2=yYY591A2_S16
SampleID3=yYY591A_S7
SampleID4=yYY591B1_S17
SampleID5=yYY591B2_S18
SampleID6=yYY591B_S8
SampleID7=yYY591C1_S19
SampleID8=yYY591C2_S20


#in=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/Updated_V0903/UsedMutants/WT/Finished_blast

in=/project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/rtt105
softpath=/scratch/ch220812/Software


Chr3M=${softpath}/iMMBIR/Database/Chr3.masked.fasta

genomeseq=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa


cd ${in}


mkdir ${SampleID}CombinedMATaInsertionSix
cd ${in}/${SampleID}CombinedMATaInsertionSix

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/rtt105/MAT_detection_for_eight_samples_3_4bp.pl -i ${in} -m ${SampleID1} -n ${SampleID2}  -p ${SampleID3}  -q ${SampleID4} -f ${SampleID5} -g  ${SampleID6} -a  ${SampleID7} -b ${SampleID8} -o ${SampleID}

perl -ne '{chomp; my @array=split/\t/,$_; next if ($array[0] eq "MATsequence"); $n++;my $id=join ":",@array[1..10],$n; print ">$id\n$array[0]\n"}' ${SampleID}.MAT.txt >${SampleID}.MAT.fasta

#bedtools maskfasta -fi Chr3.fasta -bed Masked.bed -fo Chr3.masked.fasta
#bwa index Chr3.masked.fasta
bwa mem -A2 -E1 ${Chr3M} ${SampleID}.MAT.fasta >${SampleID}.MAT.sam

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/rtt105/MAT_detection_for_deduplication_eight_2bp_v3.pl  -i ${SampleID}.MAT.sam -g ${SampleID}.MAT.fasta -o ${SampleID}.uniq -t ${SampleID}



SampleID2BP=Rtt105_2BP

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/rtt105/MAT_detection_for_eight_samples_2bp.pl -i ${in} -m ${SampleID1} -n ${SampleID2}  -p ${SampleID3}  -q ${SampleID4} -f ${SampleID5} -g  ${SampleID6} -a  ${SampleID7} -b ${SampleID8} -o ${SampleID2BP}

perl -ne '{chomp; my @array=split/\t/,$_; next if ($array[0] eq "MATsequence"); $n++;my $id=join ":",@array[1..10],$n; print ">$id\n$array[0]\n"}' ${SampleID2BP}.MAT.txt >${SampleID2BP}.MAT.fasta

#bedtools maskfasta -fi Chr3.fasta -bed Masked.bed -fo Chr3.masked.fasta
#bwa index Chr3.masked.fasta
bwa mem -A2 -E1 ${Chr3M} ${SampleID2BP}.MAT.fasta >${SampleID2BP}.MAT.sam

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/CombinationWTRAD27DNA2/Sample_V0808/rtt105/MAT_detection_for_deduplication_eight_2bp_v3.pl  -i ${SampleID2BP}.MAT.sam -g ${SampleID2BP}.MAT.fasta -o ${SampleID2BP}.uniq -t ${SampleID2BP}


SampleIDF=Rtt105_2_4BP
cat ${SampleID}.uniq.final.txt ${SampleID2BP}.uniq.final.txt > ${SampleIDF}.uniq.final.txt
