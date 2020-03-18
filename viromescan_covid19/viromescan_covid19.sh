#! /bin/bash

set -e

if [ $# -eq 0 ]
then
    echo -e "\n$0 -1 <INPUT_FASTQ_paired_end1> -2 [INPUT_FASTQ_paired_end2] -p [N_THREADS] -m [VIROMESCAN_PATH] -o <OUTPUT_DIR>
    
    -1/--input1: .fastq file containing the sequences (paired end 1)  (MANDATORY)
    -2/--input2: .fastq file containing the sequences (paired end 2) (if available)
    -p/--n_threads: number of threads to launch (default: 1)
    -m/--viromescan_path: pathway to viromescan folder (default: working directory)
    -o/--output: output directory (MANDATORY)
    "
    exit 1
fi


# -----------------------
INPUT_FASTQ_paired_end1=''; INPUT_FASTQ_paired_end2=''; N_THREADS=-1; VIROMESCAN_PATH=''; OUTPUT_DIR='';
 

OPTS=$(getopt -o 1:2:d:p:m:o: -l input1:,input2:,n_threads:,viromescan_path:,output:, -- "$@")
eval set -- "$OPTS"

while [ $# -gt 0 ] ; do
    case "$1"
    in
        -1|--input1) INPUT_FASTQ_paired_end1=$2; shift;;
        -2|--input2) INPUT_FASTQ_paired_end2=$2; shift;;
        -p|--n_threads) N_THREADS=$2; shift;;
	-m|--viromescan_path) VIROMESCAN_PATH=$2; shift;;
        -o|--output) OUTPUT_DIR=$2; shift;;
        (--) shift; break;;
    esac
    shift
done

# If the variable "N_THREADS" is not defined in input, we set the default as:
#if [ $N_THREADS -eq -1 ]
#then
  #N_THREADS=$(echo(1)|bc)
#fi


# Computing the analysis

mkdir $OUTPUT_DIR

if [ "$INPUT_FASTQ_paired_end2" == "" ]
then
bowtie2 -x $VIROMESCAN_PATH/viromescan/database/bowtie2/human_ALL+covid19 -q $INPUT_FASTQ_paired_end1 --sensitive-local --no-unal -S $OUTPUT_DIR/$OUTPUT_DIR.sam -p $N_THREADS
else
bowtie2 -x $VIROMESCAN_PATH/viromescan/database/bowtie2/human_ALL+covid19 -1 $INPUT_FASTQ_paired_end1 -2 $INPUT_FASTQ_paired_end2 --sensitive-local --no-unal -S $OUTPUT_DIR/$OUTPUT_DIR.sam -p $N_THREADS
fi

cat $OUTPUT_DIR/$OUTPUT_DIR.sam |grep -v ^@| awk '{print"@"$1"\n"$10"\n+\n"$11}' > $OUTPUT_DIR/$OUTPUT_DIR-virusnofiltr.fastq

mkdir $OUTPUT_DIR/tmp

bash $VIROMESCAN_PATH/viromescan/tools/bmtagger.sh -b $VIROMESCAN_PATH/viromescan/database/hg19/hg19reference.bitmask -x $VIROMESCAN_PATH/viromescan/database/hg19/hg19reference.srprism -T $OUTPUT_DIR/tmp/ -q1 -1 $OUTPUT_DIR/$OUTPUT_DIR-virusnofiltr.fastq -X -o $OUTPUT_DIR/$OUTPUT_DIR-nonhuman

java -jar $VIROMESCAN_PATH/viromescan/tools/picard-tools-1.71/FastqToSam.jar FASTQ=$OUTPUT_DIR/$OUTPUT_DIR-nonhuman.fastq QUALITY_FORMAT=Standard OUTPUT=$OUTPUT_DIR/$OUTPUT_DIR-nonhuman.bam SAMPLE_NAME=$OUTPUT_DIR-nonhuman

perl $VIROMESCAN_PATH/viromescan/tools/trimBWAstyle.usingBam.pl -f $OUTPUT_DIR/$OUTPUT_DIR-nonhuman.bam -q 3 -o 65

cat $OUTPUT_DIR/$OUTPUT_DIR-nonhuman.trimmed.1.fastq  $OUTPUT_DIR/$OUTPUT_DIR-nonhuman.trimmed.2.fastq $OUTPUT_DIR/$OUTPUT_DIR-nonhuman.trimmed.singleton.fastq > $OUTPUT_DIR/$OUTPUT_DIR-filter-human-quality.fastq

bash $VIROMESCAN_PATH/viromescan/tools/bmtagger.sh -b $VIROMESCAN_PATH/viromescan/database/Bacteria_custom/bacteria_custom.bitmask -x $VIROMESCAN_PATH/viromescan/database/Bacteria_custom/bacteria_custom.srprism -T $OUTPUT_DIR/tmp -q1 -1 $OUTPUT_DIR/$OUTPUT_DIR-filter-human-quality.fastq -X -o $OUTPUT_DIR/$OUTPUT_DIR-filter-human-quality-bacteria

bowtie2 -x $VIROMESCAN_PATH/viromescan/database/bowtie2/human_ALL+covid19 -q $OUTPUT_DIR/$OUTPUT_DIR-filter-human-quality-bacteria.fastq --very-sensitive-local --no-unal -S $OUTPUT_DIR/$OUTPUT_DIR-final.sam -p $N_THREADS

mkdir $OUTPUT_DIR/covid19

cat $OUTPUT_DIR/$OUTPUT_DIR-final.sam |grep -v ^@| grep "NC_002645.1" | awk '{print"@"$1"\n"$10"\n+\n"$11}' >> $OUTPUT_DIR/covid19/corona.fastq
cat $OUTPUT_DIR/$OUTPUT_DIR-final.sam |grep -v ^@| grep "NC_006577.2" | awk '{print"@"$1"\n"$10"\n+\n"$11}' >> $OUTPUT_DIR/covid19/corona.fastq
cat $OUTPUT_DIR/$OUTPUT_DIR-final.sam |grep -v ^@| grep "NC_005831.2" | awk '{print"@"$1"\n"$10"\n+\n"$11}' >> $OUTPUT_DIR/covid19/corona.fastq
cat $OUTPUT_DIR/$OUTPUT_DIR-final.sam |grep -v ^@| grep "NC_005147.1" | awk '{print"@"$1"\n"$10"\n+\n"$11}' >> $OUTPUT_DIR/covid19/corona.fastq
cat $OUTPUT_DIR/$OUTPUT_DIR-final.sam |grep -v ^@| grep "NC_012950.1" | awk '{print"@"$1"\n"$10"\n+\n"$11}' >> $OUTPUT_DIR/covid19/corona.fastq
cat $OUTPUT_DIR/$OUTPUT_DIR-final.sam |grep -v ^@| grep "NC_004718.3" | awk '{print"@"$1"\n"$10"\n+\n"$11}' >> $OUTPUT_DIR/covid19/corona.fastq
cat $OUTPUT_DIR/$OUTPUT_DIR-final.sam |grep -v ^@| grep "NC_045512.2" | awk '{print"@"$1"\n"$10"\n+\n"$11}' >> $OUTPUT_DIR/covid19/corona.fastq


samtools view -bS $OUTPUT_DIR/$OUTPUT_DIR-final.sam > $OUTPUT_DIR/$OUTPUT_DIR-final.bam

samtools sort $OUTPUT_DIR/$OUTPUT_DIR-final.bam -o $OUTPUT_DIR/$OUTPUT_DIR-final.sorted.bam

samtools index $OUTPUT_DIR/$OUTPUT_DIR-final.sorted.bam

samtools idxstats $OUTPUT_DIR/$OUTPUT_DIR-final.sorted.bam > $OUTPUT_DIR/final.genes.txt

mkdir $OUTPUT_DIR/results

cp $VIROMESCAN_PATH/viromescan/var/viromescan.R $OUTPUT_DIR

cp $VIROMESCAN_PATH/viromescan/var/HumanDNAcomplete.txt $OUTPUT_DIR 

cp $VIROMESCAN_PATH/viromescan/var/HumanALLcomplete.txt $OUTPUT_DIR 

cp $VIROMESCAN_PATH/viromescan/var/VirusALLcomplete.txt $OUTPUT_DIR 

cp $VIROMESCAN_PATH/viromescan/var/VirusDNAcomplete.txt $OUTPUT_DIR

cp $VIROMESCAN_PATH/viromescan/var/HumanALL+covid19complete.txt $OUTPUT_DIR

cp $VIROMESCAN_PATH/viromescan/var/covid19.R $OUTPUT_DIR/covid19

cd $OUTPUT_DIR/results

Rscript ../viromescan.R

cd ../covid19

bowtie2 -x $VIROMESCAN_PATH/viromescan/database/bowtie2/coronavirus -q corona.fastq --very-sensitive-local --no-unal -S coronavirus.sam -p $N_THREADS
samtools view -bS coronavirus.sam > coronavirus.bam
samtools sort coronavirus.bam -o coronavirus.sorted.bam
samtools index coronavirus.sorted.bam
samtools calmd -b coronavirus.sorted.bam --reference $VIROMESCAN_PATH/viromescan/database/bowtie2/coronavirus.fa > coronavirus.calmd.bam
samtools view -h coronavirus.calmd.bam > coronavirus.calmd.sam
samtools idxstats coronavirus.sorted.bam > coronavirus.genes.txt

Rscript covid19.R

cd ../../

mv $OUTPUT_DIR/results .
mv $OUTPUT_DIR/covid19 .

rm -rf $OUTPUT_DIR/*


mv results/* $OUTPUT_DIR/
mv covid19 $OUTPUT_DIR/

rm -rf results
rm -rf $OUTPUT_DIR/covid19/covid19.R
