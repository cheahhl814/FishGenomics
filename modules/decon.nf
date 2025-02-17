process ganonClassify {
    tag "Classify reads with Ganon"
    publishDir "./results/decon", mode: 'copy', overwrite: false, pattern: '**'

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
    ganon classify --db-prefix ${database_id} --single-reads - --output-one --output-prefix ganonClassify --threads ${task.cpus}
    """
}

process decon {
    tag ""
}

cat ${one} | cut -f 1 | sort | uniq | seqkit grep -v -f - ${fastq} > ${decon_fastq}
