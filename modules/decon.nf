process ganonClassify {
    tag "Classify reads with Ganon"
    publishDir "./results/pre-assembly/decon", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(fastq)
    path(tsv)

    output:
    path "*.tre", emit: tre
    path "*.rep", emit: rep
    path "*.one", emit: one

    script:
    def database_id = tsv.baseName
    """
    ganon build-custom --input-target file --level file --taxonomy skip --input-file ${tsv} --db-prefix ${database_id} --threads ${task.cpus}
    ganon classify --db-prefix ${database_id} --single-reads ${fastq} --output-one --output-prefix ganonClassify --threads ${task.cpus}
    """
}

process decon {
    tag "Filter out contaminant reads"\
    publishDir "./results/pre-assembly/decon", mode: 'copy', overwrite: false, pattern: '**'

    input:
    val(ones)
    path(fastq)

    output:
    path "*_decontaminated.fastq", emit: deconFastq

    script:
    def one = ones.join{" "}
    def sample_id = fastq.baseName
    """
    cat ${one} | cut -f 1 | sort | uniq | seqkit grep -v -f - ${fastq} > ${sample_id}_decontaminated.fastq
    """
}