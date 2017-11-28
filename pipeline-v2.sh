#!/bin/sh

# Thomas Pranzatelli
# 6/23/16
# NIDCR
# A master bash script for a footprinting pipeline.
#
# This part of the script accepts the four input variables:
# -d, -e, -i, -s

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -r|--reads)
    READS="$2"
    shift # past argument
    ;;
    -e|--encsr)
    ENCSR="$2"
    shift # past argument
    ;;
    -g|--genome)
    GENOM="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
    # unknown option
    ;;
esac
shift # past argument or value
done

EXPERIMENT=/data/ChioriniCompCor/metamachine/output/$ENCSR
DATA=/data/ChioriniCompCor/metamachine/ENCODE/$ENCSR
GENOME=/data/ChioriniCompCor/metamachine/genomes/$GENOM

# Here we load the different modules.
module load bowtie/2-2.2.9
module load samtools
module load homer/4.8.2
module load bedtools/2.25.0
module load rgt
module load bedops

if [[ $READS == 1 ]]
    then
    # .fastq files are concatenated.
    cat $DATA/*.fastq.gz > $EXPERIMENT/cat.fastq.gz
    echo "Single-end reads were concatenated."
    wc -l $EXPERIMENT/cat.fastq.gz
    # These files are passed to bowtie2 to be aligned to the given genome index.
    bowtie2 -x $GENOME/Bowtie2Index/genome -p 56 -t -q -U $EXPERIMENT/cat.fastq.gz -S $EXPERIMENT/bowtie2.sam
    echo "Single-end reads were aligned end-to-end."
    rm $EXPERIMENT/cat.fastq.gz
elif [[ $READS == 2 ]]
    then
    # .fastq files in each paired end folder are concatenated.
    cat $DATA/P1/*.fastq.gz > $EXPERIMENT/P1.fastq.gz
    cat $DATA/P2/*.fastq.gz > $EXPERIMENT/P2.fastq.gz
    echo "Paired-end reads were concatenated."
    wc -l $EXPERIMENT/P1.fastq.gz
    # These files are passed to bowtie2 to be aligned to the given genome index.
    bowtie2 -x $GENOME/Bowtie2Index/genome -p 56 -t -q -1 $EXPERIMENT/P1.fastq.gz -2 $EXPERIMENT/P2.fastq.gz -S $EXPERIMENT/bowtie2.sam
    echo "Paired-end reads were aligned end-to-end."
    rm $EXPERIMENT/P1.fastq.gz
    rm $EXPERIMENT/P2.fastq.gz
else
    echo "Use a read count argument of -r 1 for single-end reads and a read count argument of -r 2 for paired-end reads."
    exit
fi

# This alignment is converted to binary, and then sorted and indexed.
samtools view -u -b $EXPERIMENT/bowtie2.sam > $EXPERIMENT/view.bam
samtools sort -@ 55 $EXPERIMENT/view.bam -o $EXPERIMENT/picard.bam
samtools index $EXPERIMENT/picard.bam
rm $EXPERIMENT/bowtie2.sam
rm $EXPERIMENT/view.bam

# Peaks are estimated from the .bam files by HOMER.
makeTagDirectory $EXPERIMENT/tagDirectory/ $EXPERIMENT/picard.bam
findPeaks $EXPERIMENT/tagDirectory/ -region -size 50 -minDist 50 -o auto -tbp 0
echo "Open chromatin peaks found."
# These peaks are converted to .bed files and sorted and merged.
pos2bed.pl $EXPERIMENT/tagDirectory/peaks.txt > $EXPERIMENT/peaks.bed
echo "Peak .bed file produced."
bedtools sort -i $EXPERIMENT/peaks.bed > $EXPERIMENT/peaks.sorted.bed
echo "Peak .bed file sorted."
grep -v "#" $EXPERIMENT/tagDirectory/peaks.txt > $EXPERIMENT/tagDirectory/peaks2.txt
mv $EXPERIMENT/tagDirectory/peaks2.txt $EXPERIMENT/tagDirectory/peaks.txt
cut -f2-4,8 $EXPERIMENT/tagDirectory/peaks.txt > $EXPERIMENT/peaks.bedgraph
echo "Peak .bedgraph file produced."
rm $EXPERIMENT/peaks.bed
bedtools merge -i $EXPERIMENT/peaks.sorted.bed > $EXPERIMENT/peaks.merged.bed
rm $EXPERIMENT/peaks.sorted.bed
echo "Peaks merged and sorted for interval estimation."

# The .bed and .bam files are used by Wellington to estimate footprints.
cd $EXPERIMENT
cp /data/ChioriniCompCor/metamachine/patroloffice/experiment_matrix.bed .
rm -R $EXPERIMENT/footprints
mkdir /lscratch/$SLURM_JOBID/footprints
rgt-hint --output-location /lscratch/$SLURM_JOBID/footprints/ --estimate-bias-correction --organism $GENOM experiment_matrix.bed || fail "HINT-BC failed to run to completion."
cp -R /lscratch/$SLURM_JOBID/footprints $EXPERIMENT/footprints || fail "Transmitting the footprints from lscratch failed to run to completion."
echo "Footprint intervals produced."
rm -r $EXPERIMENT/tagDirectory
rm $EXPERIMENT/peaks.merged.bed
rm $EXPERIMENT/picard.bam
rm $EXPERIMENT/picard.bam.bai

# Produce the bedfiles of events out of PWMs for the ROC plot (if ChIP-seq
# validation data exists.)
if [[ -d $DATA/ChIP-validation ]]
    then cp -R $DATA/ChIP-validation $EXPERIMENT/ChIP-validation
    cp $EXPERIMENT/footprints/*.bed $EXPERIMENT/footprints.bed
    for chip in $EXPERIMENT/ChIP-validation/*
        do bedtools intersect -a $chip/segBarozzi.bed -b $EXPERIMENT/footprints.bed -wa -wb > $chip/positives.bed
        bedtools intersect -a $chip/segBarozzi.bed -b $EXPERIMENT/footprints.bed -v > $chip/negatives.bed
        bedtools intersect -a $chip/positives.bed -b $chip/sortedchip.bed -wa > $chip/true_positives.bed
        bedtools intersect -a $chip/positives.bed -b $chip/sortedchip.bed -v > $chip/false_positives.bed
        bedtools intersect -a $chip/negatives.bed -b $chip/sortedchip.bed -wa > $chip/false_negatives.bed
        bedtools intersect -a $chip/negatives.bed -b $chip/sortedchip.bed -v > $chip/true_negatives.bed
    done
    python /data/pranzatellitj/tools/ROC.py -e $EXPERIMENT -f -5
    echo "AUC of the ROC plot produced."
    rm -R $EXPERIMENT/ChIP-validation
else
    echo "No Chip-seq data found."
fi

# The footprints are mapped to PWMS, and then to the regulatory regions.
cp $EXPERIMENT/footprints/*.bed $EXPERIMENT/Wellington.bed
sort -k 1,1 -k2,2n $EXPERIMENT/Wellington.bed > $EXPERIMENT/sorted_Wellington.bed
bedmap --ec --multidelim '&' --skip-unmapped --sweep-all --bp-ovr 3 --echo --echo-map $EXPERIMENT/sorted_Wellington.bed $GENOME/*fimo_results.bed > $EXPERIMENT/bedmap.bed
python /data/pranzatellitj/tools/associate_footprints_hint.py -e $EXPERIMENT -g $GENOME
bedmap --ec --delim '@' --multidelim '%%' --skip-unmapped --sweep-all --echo --echo-map $GENOME/tss.tssified.bed $EXPERIMENT/associated_footprints.bed > $EXPERIMENT/association_matrix_p.bed
bedmap --ec --delim '@' --multidelim '%%' --skip-unmapped --sweep-all --echo --echo-map $GENOME/sorted_compcor_* $EXPERIMENT/associated_footprints.bed > $EXPERIMENT/association_matrix_e.bed
python /data/pranzatellitj/tools/package_everything_as_JSON.py -e $EXPERIMENT/association_matrix_e.bed -p $EXPERIMENT/association_matrix_p.bed -t $GENOME/tss.tssified.bed -o $EXPERIMENT
rm $EXPERIMENT/association_matrix_e.bed $EXPERIMENT/association_matrix_p.bed $EXPERIMENT/associated_footprints.bed $EXPERIMENT/sorted_Wellington.bed $EXPERIMENT/Wellington.bed $EXPERIMENT/bedmap.bed
