import csv
import os, sys
import shutil
import subprocess
from pathlib import Path

# read each line from the csv file header coloumn sampleid,cram,vcf and copy the nextflow job script run.sh as sample specific
with open('input.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        file_name ='{0}.yml'.format(row['sampleid'])
        dir_name = (row['sampleid'])
        os.mkdir( dir_name )
        with open(os.path.join(dir_name, file_name), 'w') as f:
            f.write('sample_id: ' + row['sampleid'] + '\ncram: ' + row['cram'] + '\nvcf: ' + row['vcf']) # write sample id, bam/cram and vcf files path
            shutil.copy2('run.sh', dir_name) # copy nextflow job submission script
            file = Path(dir_name, 'run.sh')
            file.write_text(file.read_text().replace('NNNNN', row['sampleid'])) # replace the sample string with sample id
