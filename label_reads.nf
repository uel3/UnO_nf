nextflow.enable.dsl=2

// Pipeline Input parameters

params.reads = "$HOME/coal_reads/*{1,2}_paired.fq.gz" //this allows me to assign the meta.id to all reads in the params.reads directory-this will be helpful 

ch_raw_short_reads = Channel
            .fromFilePairs(params.reads)
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
            .map { row ->
                        def sample = [:]
                        sample.id           = row[0]
                        return [ sample, row[1] ]
                }
ch_short_reads_grouped = ch_raw_short_reads //this would be the trimmed_reads_ch
        .map { sample, reads -> [sample.group, sample, reads] }
        .groupTuple(by: 0)
        .map { group, samples, reads ->
            def sample = [:]
            sample.id = "group-$group"
            sample.group = group
            def reads1 = reads.collect { it[0] }
            def reads2 = reads.collect { it[1] }
            [sample, reads1, reads2]
    }
process CHECK_READS{
    input:
    tuple val(sample), path( all_reads )
   
    output:
    stdout

    script:
    """
    echo ${sample.id}
    """
}
process COLLECT_READS{
    input:
    tuple val(sample), path(reads1), path(reads2)
    
    output:
    stdout
    
    script:
    def input = "-1 \"" + reads1.join(",") + "\" -2 \"" + reads2.join(",") + "\""
    """
    echo $input 
    """
    //this outputs only the *_1 reads which is what I want 
}
process CHECKR1ANDR2{
    input:
    tuple val(sample), path(reads1), path(reads2)
    output:
    stdout
    script:
    """
    echo $reads1 $reads2
    """
}
workflow {

        CHECK_READS( ch_raw_short_reads )
        COLLECT_READS(  ch_short_reads_grouped ) //this works to properly format for coassembly 
        CHECKR1ANDR2( ch_short_reads_grouped )
}