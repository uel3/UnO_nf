nextflow.enable.dsl=2

// Pipeline Input parameters

params.reads = "$HOME/coal_reads/*{1,2}_paired.fq.gz" //this allows me to assign the meta.id to all reads in the params.reads directory-this will be helpful 

ch_raw_short_reads = Channel
            .fromFilePairs(params.reads)
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
            .map { row ->
                        def sample = [:]
                        sample.id  = row[0]
                        sample.group = row[0]  // Use a string as the default group value
                        return [ sample, row[1] ]
                }
ch_short_reads_grouped = ch_raw_short_reads //this would be the trimmed_reads_ch
        .map { sample, reads -> [sample.group, sample, reads] }
        .groupTuple(by: 0)
        .map { group, samples, reads ->
            def groupedSample = [:]
            groupedSample.id = group //with this sample.id for reads in ch_short_reads_grouped is the SRRnum
            groupedSample.group = "paired" //with this the sample.group for reads in ch_short_reads_grouped is "parired" -I can call this?
            def reads1 = reads.collect { it[0] }
            def reads2 = reads.collect { it[1] }
            [groupedSample, reads1, reads2]
    }
process CHECK_READS{
    input:
    tuple val(sample), path( all_reads )
   
    output:
    stdout

    script:
    """
    echo 'this is from CHECK_READS' ${sample.id} ${sample.group}
    """
}
process COLLECT_READS{ //this is meant to represent read grouping and variables used beyond trimmmomatic  
    input:
    tuple val(sample), path(reads1), path(reads2)
    
    output:
    stdout
    
    script:
    def input = "-1 \"" + reads1.join(",") + "\" -2 \"" + reads2.join(",") + "\""
    """
    echo 'this is from COLLECT_READ' $sample.id $sample.group
    """
    //this outputs only the *_1 reads which is what I want 
    //reads1.baseName prints SRR_1_paired.fq
    //reads1.simpleName prints SRR_1_paired
}
process CHECKR1ANDR2{
    input:
    tuple val(sample), path(read)
    output:
    stdout
    //not ideal but this is how I got each read to be processes seperatley-not what I wanted...still need to ID a work around
    script:
    """
    echo 'this is from CHECKR1..' ${read[0]}
    echo 'this is from CHECHR2..' ${read[1]}
    """
}
workflow {

        CHECK_READS( ch_raw_short_reads )
        COLLECT_READS(  ch_short_reads_grouped ) //this works to properly format for coassembly 
        CHECKR1ANDR2( ch_raw_short_reads )
}