#!/usr/bin/env nextflow
// BRCA SV and amplification timing: SVclone (per sample, SVs + SNVs, post-assigned to the SNV model)
// and AmplificationTimeR (all focal amplification events in one task).
nextflow.enable.dsl = 2

process SVCLONE {
    tag "$sample"
    publishDir "${params.outdir}/svclone/${sample}", mode: 'copy'

    input:
    tuple val(sample), val(sex), path('tumor.bam'), path('tumor.bam.bai'), path(svs), path(cnv), path(pp), path(snv_vcf)

    output:
    tuple val(sample), path("${sample}"), path('status.tsv'), path('svclone.ini'), emit: results

    script:
    """
    run_svclone.sh ${sample} ${sex} tumor.bam ${svs} ${cnv} ${pp} ${snv_vcf} ${task.cpus} ${params.snv_subsample}
    mkdir -p ${sample}
    """
}

process AMPTIMER {
    publishDir "${params.outdir}/amplificationtimer", mode: 'copy'

    input:
    path amp_inputs

    output:
    path 'out/*'

    script:
    """
    run_amptimer.R ${amp_inputs} out ${params.amptimer_commit}
    """
}

workflow {
    if (params.run_svclone) {
        Channel.fromPath(params.svclone_manifest, checkIfExists: true)
            .splitCsv(header: true)
            .filter { row -> !params.samples || params.samples.tokenize(',').contains(row.sample) }
            .map { row ->
                tuple(row.sample, row.sex,
                      file(row.bam, checkIfExists: true), file(row.bai, checkIfExists: true),
                      file(row.svs, checkIfExists: true), file(row.cnv, checkIfExists: true),
                      file(row.purity_ploidy, checkIfExists: true), file(row.snv_vcf, checkIfExists: true))
            }
            | SVCLONE
    }
    if (params.run_amptimer) {
        AMPTIMER(file(params.amptimer_inputs, checkIfExists: true, type: 'dir'))
    }
}
