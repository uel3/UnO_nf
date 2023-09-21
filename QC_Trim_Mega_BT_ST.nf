/*
========================================================================================
   UnO Nextflow Workflow
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Pipeline Input parameters

params.outdir = 'results'
//params.genome = "${launchDir}/data/ref_genome/ecoli_rel606.fasta"
params.reads = "$HOME/coal_reads/*_{1,2}.fastq.gz"
//params.adapter = "$HOME/metagenome_practice/trimmomatic_output/TruSeq3-PE.fa"

println """\
         U n O - N F   P I P E L I N E
         ===================================
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

/*
========================================================================================
   Create Channels
========================================================================================
*/

//ref_ch = Channel.fromPath( params.genome, checkIfExists: true )  
reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true ) 
//trimmed_reads_ch = Channel.fromFilePairs(TRIMMOMATIC.out.trimmed_reads, checkIfExists: true )
//adapter_ch = Channel.fromPath( params.adapter, checkIfExists: true )

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

    raw_reads_ch = FASTQC_RAW( reads_ch )
    trimmed_reads_ch = TRIMMOMATIC( reads_ch )
    FASTQC_TRIMMED( trimmed_reads_ch.trimmed_reads )
    megahit_assembly_ch = MEGAHIT( trimmed_reads_ch.trimmed_reads )
    //MIDAS2_TRIMMED ( TRIMMOMATIC.out.trimmed_reads )
    bt2_index_ch = BOWTIE2_INDEX( megahit_assembly_ch.megahit_contigs ) // https://www.nextflow.io/docs/latest/process.html#understand-how-multiple-input-channels-work
    mapped_reads_ch = BOWTIE2_MAP_READS( bt2_index_ch.bowtie2_index, trimmed_reads_ch.trimmed_reads )
    //SAMTOOLS_SORT( mapped_reads_ch.aligned_sam )
    // Enter the rest of the processes for variant calling based on the bash script below

}

/*
========================================================================================
   Processes
========================================================================================
*/

/*
 * QC raw fastq reads.
 */
process FASTQC_RAW {
    tag{"FASTQC ${reads}"}
    label 'process_low'

    publishDir("${params.outdir}/fastqc_raw", mode: 'copy')
    
    input:
    tuple val( sample_id ), path( reads )

    output:
    path( "*_fastqc*" ), emit: raw_reads

    script:
    """
    fastqc ${reads}
    """

    stub:
    """
    touch ${sample_id}_1_fastqc.html
    touch ${sample_id}_1_fastqc.gz
    touch ${sample_id}_2_fastqc.html
    touch ${sample_id}_2_fastqc.gz
    """
}

/*
 * Trimming raw fastq reads.
 */
process TRIMMOMATIC {
    tag{"TRIMMOMATIC ${reads}"}
    label 'process_low'

    publishDir("${params.outdir}/trimmed_reads", mode: 'copy')
  
    input:
    tuple val( sample_id ), path( reads )

    output:
    tuple val( sample_id ), path( "*.trimmed.fq.gz" ), emit: trimmed_reads
    tuple val( sample_id ), path ("*.unpaired.fq.gz"), emit: unpaired


    script:
     """
     trimmomatic PE -threads 10 -phred33 ${reads} ${sample_id}_1.trimmed.fq.gz ${sample_id}_1.unpaired.fq.gz ${sample_id}_2.trimmed.fq.gz ${sample_id}_2.unpaired.fq.gz ILLUMINACLIP:TrueSeq3-PE.fa:20:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
     """

    stub:
     """
     touch ${sample_id}_1.trimmed.fq.gz
     touch ${sample_id}_2.trimmed.fq.gz
     touch ${sample_id}_1.unpaired.fq.gz
     touch ${sample_id}_2.unpaired.fq.gz
     """
}

/*
 * QC trimmed fastq reads.
 */
process FASTQC_TRIMMED {
    tag{"FASTQC ${reads_trimmed}"}
    label 'process_low'

    publishDir("${params.outdir}/fastqc_trimmed", mode: 'copy')
    
    input:
    tuple val( sample_id ), path( reads_trimmed )

    output:
    path( "*_fastqc*" )

    script:
    """
    fastqc ${reads_trimmed}
    """

    stub:
    """
    touch ${sample_id}_1_trimmed_fastqc.html
    touch ${sample_id}_2_trimmed_fastqc.html
    touch ${sample_id}_1_trimmed.fastqc.gz
    touch ${sample_id}_2_trimmed.fastqc.gz
    """
    
}

/*
 * Metagenomic co-assembly of reads with megahit.
 */
process MEGAHIT {
    tag{"MEGAHIT ${reads_trimmed}"}
    label 'process_low'
    
    publishDir("${params.outdir}/megahit_out", mode: 'copy')
    
    input:
    tuple val( sample_id ), path( reads_trimmed )

    output:
    tuple val( sample_id ), path("megahit_out/*.contigs.fa")                            , emit: megahit_contigs
    tuple val( sample_id ), path("megahit_out/intermediate_contigs/k*.contigs.fa")      , emit: megahit_k_contigs
    tuple val( sample_id ), path("megahit_out/intermediate_contigs/k*.addi.fa")         , emit: megahit_addi_contigs
    tuple val( sample_id ), path("megahit_out/intermediate_contigs/k*.local.fa")        , emit: megahit_local_contigs
    tuple val( sample_id ), path("megahit_out/intermediate_contigs/k*.final.contigs.fa"), emit: megahit_kfinal_contigs


    script:
    """
    megahit --presets meta-large -1 ${reads_trimmed[0]} -2 ${reads_trimmed[1]} -t 8 --out-prefix ${sample_id}
    """
    
    stub:
    """
    mkdir megahit_out
    touch megahit_out/stub.contigs.fa
    mkdir megahit_out/intermediate_contigs
    touch megahit_out/intermediate_contigs/kstub.contigs.fa
    touch megahit_out/intermediate_contigs/kstub.addi.fa
    touch megahit_out/intermediate_contigs/kstub.local.fa
    touch megahit_out/intermediate_contigs/kstub.final.contigs.fa
    """
}

/*
 * Creating an Index from the assembly to map reads to assembly.
 */
process BOWTIE2_INDEX {
    tag{"BOWTIE2_INDEX ${assembly}"}
    label 'process_low'
    
    publishDir("${params.outdir}/bowtie2_out", mode: 'copy')

    input:
    tuple val( sample_id ), path( assembly )

    output:
    path ( "${sample_id}_index*.bt2" ), emit: bowtie2_index 
    
    script:
    """
    bowtie2-build ${assembly} ${sample_id}_index
    touch ${sample_id}_index
    """

    stub:
    """
    mkdir bowtie2_out
    touch ${sample_id}_index.1.bt2
    touch ${sample_id}_index.2.bt2
    touch ${sample_id}_index.3.bt2
    touch ${sample_id}_index.4.bt2
    touch ${sample_id}_index.rev.1.bt2
    touch ${sample_id}_index.rev.2.bt2
    touch ${sample_id}_index
    """
}
/* 
* Mapping reads to bowtie2 index to evaluate coverage and for downstreaming binning tools.
*/
process BOWTIE2_MAP_READS{
    tag "BOWTIE2_MAP_READS ${index} ${reads_trimmed}" //"$assembler-$name"
    publishDir ("${params.outdir}/bowtie2_out/mapped", mode: 'copy') //Assembly/${assembler}/${name}_QC", mode: params.publish_dir_mode,
        //saveAs: {filename -> filename.indexOf(".bowtie2.log") > 0 ? filename : null}

    input:
    path( index )
    tuple val( sample_id ), path( reads_trimmed )
    //val   save_unaligned
    //val   sort_bam
    
    
    output:
    //path( "${sample_id}.sam" ), emit: aligned_sam
    path( "${sample_id}_sorted.bam*" ), emit: aligned_bam 
    //path( "${sample_id}.bowtie2.log" ), emit: align_log
    //path( "${sample_id}.bam.bai" ), emit: index_bam

    script:
    def idx = index[0].getBaseName(2)
    """
    bowtie2 -p 8 -x ${idx} -q -1 ${reads_trimmed[0]} -2 ${reads_trimmed[1]} --no-unal |samtools view -@ 2 -b -S -h | samtools sort -o ${sample_id}_sorted.bam 
    samtools index ${sample_id}_sorted.bam
    """
    //required the -h flag to make this work into view/sort commands
    //needed to dealre the correct output in declared output
    stub:
    """
    mkdir mapped
    touch ${sample_id}_sorted.bam
    touch ${sample_id}_sorted.bam.bai
    """
}
/*
 * MaxBin2 binning of assembled contigs
 */
process SAMTOOLS_SORT{
    tag "SAMTOOLS_SORT ${reads_mapped}"
    publishDir ("${params.outdir}/bowtie2_out/mapped", mode: 'copy')

    input:
    path( reads_mapped )

    output:
    path( "${reads_mapped.simpleName}_sorted.bam*" ), emit: aligned_bam  

    script:
    """
    samtools view -@ 2 -b -S -h ${reads_mapped}| samtools sort -o ${reads_mapped.simpleName}_sorted.bam 
    samtools index -b ${reads_mapped.simpleName}_sorted.bam
    """
    stub:
    """
    touch ${reads_mapped.simpleName}.bam
    touch ${reads_mapped.simpleName}.bam.bai
    """
}
/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
       Pipeline execution summary
       ---------------------------
       Completed at: ${workflow.complete}
       Duration    : ${workflow.duration}
       Success     : ${workflow.success}
       workDir     : ${workflow.workDir}
       exit status : ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
}
/*
========================================================================================
   THE END
========================================================================================
*/
