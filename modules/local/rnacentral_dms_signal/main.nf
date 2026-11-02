process RNACENTRAL_DMS_SIGNAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(treated_rc), path(untreated_rc, stageAs: 'untreated/*')

    output:
    tuple val(meta), path("${prefix}.dms_signal.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    rnacentral_dms_signal.py \
        --sample "${meta.id}" \
        --treated "${treated_rc}" \
        --untreated "${untreated_rc}" \
        --output "${prefix}.dms_signal.tsv" \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tskipped_sites\\tpct_a_muts\\tpct_c_muts\\tpct_g_muts\\tpct_u_muts\\n${meta.id}\\t0\\t30.00\\t30.00\\t20.00\\t20.00\\n' > ${prefix}.dms_signal.tsv
    """
}
