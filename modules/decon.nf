process ganonClassify {
    tag "Classify reads with Ganon"
    publishDir "./results/decon", mode: 'copy', overwrite: false, pattern: '**'

    input:
    val(fastqs)
    path(tsv)

    output:
    path "*.tre", emit: tre
    path "*.rep", emit: rep
    path "*.one", emit: one

    script:
    def database_id = tsv.baseName
    def fastq = fastqs.join(" ")
    """
    ganon build-custom --input-target file --level file --taxonomy skip --input-file ${tsv} --db-prefix ${database_id} --threads ${task.cpus}
    cat ${fastq} | ganon classify --db-prefix ${database_id} --single-reads - --output-one --output-prefix ganonClassify --threads ${task.cpus}
    """
}

process ganonReport {
    tag "Report Ganon's classification results"
    publishDir "./results/decon", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(tsv)
    path(classifyRep)

    output:
    path "*.tre", emit: tre

    script:
    def database_id = tsv.baseName
    """
    ganon report --db-prefix ${database_id} --input ${classifyRep} --output-prefix ganonReport --report-type reads
    """
}
