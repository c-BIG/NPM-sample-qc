

=======================
NPM-sample-qc resources 
=======================

autosomes_non_gap_regions.bed
==============================

The white-listed region file 'autosomes_non_gap_regions.bed' derived from the following.

* Create autosomes bed file from reference genome index (fai)

``head -22 ${fai} |awk '{print \$1"\t0""\t"\$2}' > autosomes.bed``

* Download gap regions from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz

``zcat gap.txt.gz |cut -f2-4 -| egrep -w '^chr[1-9]|chr[1-2][0-9]' |sort -k1,1V -k2,2n > gap_regions.bed``
 
``bedtools subtract -a autosomes.bed -b gap_regions.bed >autosomes_non_gap_regions.bed``
