process mosdepth_datamash {

    tag { sample }

    input:
    tuple val(sample), path(regions)
    path autosomes_non_gap_regions

    output:
    tuple val(sample), path("${sample}.mosdepth.csv"), emit: coverage

    script:
    """
    # filter mosdepth outputs to write the bins that are overlap with autosomes non gap and N bases

    zcat "${sample}.regions.bed.gz" | bedtools intersect -a stdin -b ${autosomes_non_gap_regions} | gzip -9c > "${sample}.regions.autosomes_non_gap_n_bases.bed.gz"

    # calculate metrics
    BED="${sample}.regions.autosomes_non_gap_n_bases.bed.gz";
    mean_coverage=\$(zcat \$BED | datamash --round 6 mean 4);
    sd_coverage=\$(zcat \$BED | datamash --round 6 sstdev 4);
    median_coverage=\$(zcat \$BED | datamash --round 6 median 4);
    mad_coverage=\$(zcat \$BED | datamash --round 6 madraw 4);
    total_bases=\$(zcat \$BED | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}');
    ge_1x_bases=\$(zcat \$BED | awk '\$4>=1' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}');
    ge_10x_bases=\$(zcat \$BED | awk '\$4>=10' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}');
    ge_15x_bases=\$(zcat \$BED | awk '\$4>=15' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}');
    ge_30x_bases=\$(zcat \$BED | awk '\$4>=30' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}');
    ge_40x_bases=\$(zcat \$BED | awk '\$4>=40' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}');

    # save output
    header="mean_autosome_coverage,sd_autosome_coverage,median_autosome_coverage,mad_autosome_coverage,total_autosome_bases,ge_1x_autosome_bases,ge_10x_autosome_bases,ge_15x_autosome_bases,ge_30x_autosome_bases,ge_40x_autosome_bases";
    row="\$mean_coverage,\$sd_coverage,\$median_coverage,\$mad_coverage,\$total_bases,\$ge_1x_bases,\$ge_10x_bases,\$ge_15x_bases,\$ge_30x_bases,\$ge_40x_bases";
    echo "\$header" > "${sample}.mosdepth.csv";
    echo \$row >> "${sample}.mosdepth.csv"
    """
}
