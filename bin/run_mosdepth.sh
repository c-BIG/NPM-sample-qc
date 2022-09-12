#!/bin/bash

#### parse args
display_help() {
    echo "Usage: $0 --ref_fasta=<fasta> --input_bam_cram=<bam> --gap_regions=<gz> --output_csv=<csv> --work_dir=<dir>" >&2
    echo
    exit 1
}

for i in "$@"
do
case $i in
    -f=* | --ref_fasta=*)
	REF_FASTA="${i#*=}"
        shift
        ;;
    -i=* | --input_bam_cram=*)
        INPUT_BAM_CRAM="${i#*=}"
        shift
        ;;
    -n=* | --gap_regions=*)
        GAP_REGIONS="${i#*=}"
        shift
        ;;
    -o=* | --output_csv=*)
        OUTPUT_CSV="${i#*=}"
        shift
        ;;
    -s=* | --sample_id=*)
        SAMPLE_ID="${i#*=}"
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


#### run mosdepth
mosdepth --no-per-base --by 1000 --mapq 20 --threads 4 --fasta $REF_FASTA $SAMPLE_ID $INPUT_BAM_CRAM

#### filter outputs
# focus on autosomes
head -22 "$REF_FASTA.fai" |awk '{print $1"\t0""\t"$2}' > autosomes.bed
zcat $SAMPLE_ID.regions.bed.gz | bedtools intersect -a stdin -b autosomes.bed | gzip -9c > $SAMPLE_ID.regions.autosomes.bed.gz

# exclude bins that overlap with N bases in ref
zcat $GAP_REGIONS |cut -f2-4 -|egrep -v '_|-|X|Y'|sort -k1,1V -k2,2n > gap_regions.bed 
zcat $SAMPLE_ID.regions.autosomes.bed.gz | bedtools intersect -v -a stdin -b gap_regions.bed | gzip -9c > $SAMPLE_ID.regions.autosomes_minus_n_bases.bed.gz

#### calculate metrics
BED="$SAMPLE_ID.regions.autosomes_minus_n_bases.bed.gz"
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
