#!/bin/bash

#### parse args
display_help() {
    echo "Usage: $0 --input_bam=<bam> --autosomes_bed=<bed> --n_regions_bed=<bed> --output_csv=<csv> --work_dir=<dir>" >&2
    echo
    exit 1
}

for i in "$@"
do
case $i in
    -i=* | --input_bam=*)
        INPUT_BAM="${i#*=}"
        shift
        ;;
    -a=* | --autosomes_bed=*)
        AUTOSOMES_BED="${i#*=}"
        shift
        ;;
    -n=* | --n_regions_bed=*)
        N_REGIONS_BED="${i#*=}"
        shift
        ;;
    -o=* | --output_csv=*)
        OUTPUT_CSV="${i#*=}"
        shift
        ;;
    -w=* | --work_dir=*)
        WORK_DIR="${i#*=}"
        shift
        ;;
    -h | --help)
        display_help
        exit 0
        ;;
    *)
        echo "Error: Unknown option: $1" >&2
        display_help
        exit 1
        ;;
esac
done

if [ -z "$INPUT_BAM" ] || [ -z "$AUTOSOMES_BED" ] || [ -z "$N_REGIONS_BED" ] || [ -z "$OUTPUT_CSV" ] || [ -z "$WORK_DIR" ]; then
  echo "Error: One or more variables are undefined"
  display_help
  exit 1
fi

SAMPLE_ID=$(echo $INPUT_BAM | awk -F '/' '{print $NF}' | sed 's/.bam//g')

echo "INPUT_BAM     = $INPUT_BAM"
echo "AUTOSOMES_BED = $AUTOSOMES_BED"
echo "N_REGIONS_BED = $N_REGIONS_BED"
echo "OUTPUT_CSV    = $OUTPUT_CSV"
echo "WORK_DIR      = $WORK_DIR"
echo "SAMPLE_ID     = $SAMPLE_ID"

#### run mosdepth
mosdepth --no-per-base --by 1000 --mapq 20 --threads 4 $WORK_DIR/$SAMPLE_ID $INPUT_BAM

#### filter outputs
# focus on autosomes
zcat $WORK_DIR/$SAMPLE_ID.regions.bed.gz | bedtools intersect -a stdin -b $AUTOSOMES_BED | gzip -9c > $WORK_DIR/$SAMPLE_ID.regions.autosomes.bed.gz
# exclude bins that overlap with N bases in ref
zcat $WORK_DIR/$SAMPLE_ID.regions.autosomes.bed.gz | bedtools intersect -v -a stdin -b $N_REGIONS_BED | gzip -9c > $WORK_DIR/$SAMPLE_ID.regions.autosomes_minus_n_bases.bed.gz

#### calculate metrics
BED="$WORK_DIR/$SAMPLE_ID.regions.autosomes_minus_n_bases.bed.gz"
mean_coverage=$(zcat $BED | datamash --round 6 mean 4)
sd_coverage=$(zcat $BED | datamash --round 6 sstdev 4)
median_coverage=$(zcat $BED | datamash --round 6 median 4)
mad_coverage=$(zcat $BED | datamash --round 6 madraw 4)
total_bases=$(zcat $BED | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
ge_1x_bases=$(zcat $BED | awk '$4>=1' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
ge_10x_bases=$(zcat $BED | awk '$4>=10' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
ge_15x_bases=$(zcat $BED | awk '$4>=15' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
ge_30x_bases=$(zcat $BED | awk '$4>=30' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
ge_40x_bases=$(zcat $BED | awk '$4>=40' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')

#### save output
header="mean_autosome_coverage,sd_autosome_coverage,median_autosome_coverage,mad_autosome_coverage,total_autosome_bases,ge_1x_autosome_bases,ge_10x_autosome_bases,ge_15x_autosome_bases,ge_30x_autosome_bases,ge_40x_autosome_bases"
row="$mean_coverage,$sd_coverage,$median_coverage,$mad_coverage,$total_bases,$ge_1x_bases,$ge_10x_bases,$ge_15x_bases,$ge_30x_bases,$ge_40x_bases"
echo "$header" > $OUTPUT_CSV
echo $row >> $OUTPUT_CSV

#### done
echo "DONE"
