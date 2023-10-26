process verifybamid2 {

    tag { sample }

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta
    path vbi2_ud
    path vbi2_bed
    path vbi2_mean

    output:
    tuple val(sample), path("${sample}.selfSM"), path("${sample}.Ancestry"), emit: verifybamid_bam_out
    // tuple val(sample), path("${sample}.selfSM"), emit: freemix
    // tuple val(sample), path("${sample}.Ancestry"), emit: ancestry

    script:
    def reference = ref_fasta ? /--Reference "${ref_fasta}"/ : ''

    """
   # run verifybamid2
   # get the percentage of cross-individual contamination rate 
    
    verifybamid2 --NumPC 4 --NumThread ${task.cpus*2} --SVDPrefix ${params.vbi2_svdprefix} ${reference} --BamFile ${bam} --Output ${sample}

    """
}

