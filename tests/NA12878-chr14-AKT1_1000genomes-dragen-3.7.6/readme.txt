The VCF regions extracted from NA12878.hard-filtered.vcf.gz are AKT1 gene extended regions on both ends.. chr14:101769349-106795748

bcftools view -r chr14:101769349-106795748 NA12878.hard.filtered.vcf.gz >NA12878-chr14-AKT1.vcf
bcftools view NA12878-chr14-AKT1.vcf -Oz -o NA12878-chr14-AKT1.vcf.gz

tabix -p vcf NA12878-chr14-AKT1.vcf.gz
