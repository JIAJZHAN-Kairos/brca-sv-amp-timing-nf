#!/usr/bin/env nextflow
// BRCA SV and amplification timing: SVclone (per sample, SVs + SNVs, post-assigned to the SNV model)
// and AmplificationTimeR (selected focal amplification events in one task, or every gained CN segment
// of every sample in chunks of samples).
nextflow.enable.dsl = 2

process SVCLONE {
    tag "$sample"
    publishDir "${params.outdir}/svclone/${sample}", mode: 'copy'

    input:
    tuple val(sample), val(sex), val(mean_cov), path('tumor.bam'), path('tumor.bam.bai'), path(svs), path(cnv), path(pp), path(snv_vcf)

    output:
    tuple val(sample), path("${sample}"), path('status.tsv'), path('svclone.ini'), emit: results

    script:
    """
    run_svclone.sh ${sample} ${sex} tumor.bam ${svs} ${cnv} ${pp} ${snv_vcf} ${task.cpus} ${params.snv_subsample} ${mean_cov}
    mkdir -p ${sample}
    cp status.tsv ${sample}/svclone_status.tsv   # S3 work dirs drop empty directories
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

process AMPTIMER_GW {
    tag "chunk${chunk}"
    publishDir "${params.outdir}/amplificationtimer_genomewide/chunks", mode: 'copy'

    input:
    tuple val(chunk), val(samples), path(sample_files), path('events.tsv')

    output:
    path "chunk${chunk}*.tsv"

    script:
    def check = params.amptimer_check_full ? "AMPT_FULL_MUTS=1 run_amptimer_genomewide.R . chunk${chunk}.fullmuts.tsv ${params.amptimer_commit} ${samples.join(' ')}" : ''
    """
    run_amptimer_genomewide.R . chunk${chunk}.tsv ${params.amptimer_commit} ${samples.join(' ')}
    ${check}
    """

    stub:
    """
    ls *.cn.tsv | wc -l > chunk${chunk}.tsv; echo ${samples.join(',')} >> chunk${chunk}.tsv
    """
}

workflow {
    if (params.run_svclone) {
        Channel.fromPath(params.svclone_manifest, checkIfExists: true)
            .splitCsv(header: true)
            .filter { row -> !params.samples || params.samples.tokenize(',').contains(row.sample) }
            .map { row ->
                tuple(row.sample, row.sex, row.mean_cov,
                      file(row.bam, checkIfExists: true), file(row.bai, checkIfExists: true),
                      file(row.svs, checkIfExists: true), file(row.cnv, checkIfExists: true),
                      file(row.purity_ploidy, checkIfExists: true), file(row.snv_vcf, checkIfExists: true))
            }
            | SVCLONE
    }
    if (params.run_amptimer) {
        AMPTIMER(file(params.amptimer_inputs, checkIfExists: true, type: 'dir'))
    }
    if (params.run_amptimer_gw) {
        def gw = params.amptimer_gw_inputs
        Channel.fromPath("${gw}/samples.txt", checkIfExists: true)
            .splitText().map { it.trim() }.filter { it }
            .filter { s -> !params.samples || params.samples.tokenize(',').contains(s) }
            .collate(params.amptimer_gw_chunk)
            .toList().flatMap { chunks -> chunks.withIndex().collect { c, i -> [i + 1, c] } }
            .map { chunk, samples ->
                def files = samples.collectMany { s -> ['cn', 'mult', 'muts'].collect { k -> file("${gw}/${s}.${k}.tsv", checkIfExists: true) } }
                tuple(chunk, samples, files, file("${gw}/events.tsv", checkIfExists: true))
            }
            | AMPTIMER_GW
    }
}
