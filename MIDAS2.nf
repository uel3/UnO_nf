/*
========================================================================================
   MIDAS2 Nextflow Workflow
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Pipeline Input parameters

params.outdir = 'MIDAS2'
//params.genome = "${launchDir}/data/ref_genome/ecoli_rel606.fasta"
params.reads = "$HOME/coal_reads/*_{1,2}.fastq.gz"

println """\
         M I D A S 2- N F   P I P E L I N E
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

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

    midasdb_uhgg_ch = MIDAS2_DB_BUILD ()
    midas2_species_ch = MIDAS2_SPECIES( reads_ch )
    midas2_snps_ch = MIDAS2_SNPS ( reads_ch, midas2_species_ch.species_id )
    //MIDAS2_PARSE( midas_species_ch, midas_snps_ch )
    //Enter the rest of the processes for variant calling based on the bash script below

}

/*
========================================================================================
   Processes
========================================================================================
*/
/*
 * Download UHGG Database for MIDAS2.
 */
process MIDAS2_DB_BUILD {
    tag{"MIDAS2_DB_BUILD ughh_db"}
    label 'process_low'

    publishDir("${params.outdir}", mode: 'copy') //do I want to copy this--probably not--just link it
    conda '/scicomp/home-pure/uel3/.conda/envs/midas_changed'
    
    output:
    path( "my_midasdb_uhgg/*" ), emit: uhgg_db
    path( "my_midasdb_uhgg/metadata.tsv" ), emit: db_file

    script:
    """
    midas2 database --init --midasdb_name uhgg --midasdb_dir my_midasdb_uhgg
    """

    stub:
    """
    mkdir my_midasdb_uhgg
    touch my_midasdb_uhgg/stub
    touch my_midasdb_uhgg/metadata.tsv
    """
}
/*
 * MIDAS2 run species to get list of potential species in sample. 
 */
process MIDAS2_SPECIES {
    errorStrategy 'ignore'
    tag{"MIDAS2_SPECIES ${reads}"}
    label 'process_low'

    publishDir("${params.outdir}", mode: 'copy') 
    conda '/scicomp/home-pure/uel3/.conda/envs/midas_changed'

    input:
    tuple val( sample_id ), path( reads )
    //path( uhgg_db )

    output:
    path( "midas2_output/${sample_id}/species/log.txt" )
    path( "midas2_output/${sample_id}/species/species_profile.tsv" ), emit: species_id
    path( "midas2_output/${sample_id}/temp/*" ), optional: true //adding the optional: true keeps nf from throwing error

    script: //getting an error that midas2 cannot find hs-blastn but it is in the midas_changes env located :/scicomp/home-pure/uel3/.conda/envs/midas_changed/bin/hs-blastn
    //need to include -profile conda when running the script to activate the correct environment $nextflow run MIDAS2.nf -profile conda sge
    """
    midas2 run_species \
      --sample_name ${sample_id} \
      -1 ${reads[0]} \
      -2 ${reads[1]} \
      --midasdb_name uhgg \
      --midasdb_dir my_midasdb_uhgg \
      --num_cores 8 \
      midas2_output
    """

    stub:
    """
    mkdir midas2_output
    mkdir midas2_output/${sample_id}
    mkdir midas2_output/${sample_id}/species
    touch midas2_output/${sample_id}/species/log.txt
    touch midas2_output/${sample_id}/species/species_profile.tsv
    mkdir midas2_output/${sample_id}/temp
    touch midas2_output/${sample_id}/temp/stub
    """
    // a run through of this process resulted in a command error that stopped the process-this output was '[ScoreBlkKbpUngappedCalc] Warning: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence. Please verify the query sequence(s) and/or filtering options.' 
    //this type of error should not stop the process going to add an ignore error statement to see if it will work even with the warning 
    //adding the ignore statement allows the process to run but I am not getting the correct output-required me to restructure my outputs-since the outdir is called in the script, I needed to remove it from my publishDir call but also include it expected output
    //still not getting output to the outdir-it is in the owr==work directory 
}
/*
 * MIDAS2 run snps to get narrowed down list of potential species in sample. 
 */
process MIDAS2_SNPS {
    tag{"MIDAS2_SNPS ${reads}"}
    label 'process_low'

    publishDir("${params.outdir}", mode: 'copy')
    conda '/scicomp/home-pure/uel3/.conda/envs/midas_changed'
    
    input:
    tuple val( sample_id ), path( reads )
    path( species_profile )

    output:
    path( "midas2_output/${sample_id}/snps/log.txt" )
    path( "midas2_output/${sample_id}/snps/snps_summary.tsv"), emit: midas2_snps
    path( "midas2_output/${sample_id}/snps/*.snps.tsv.lz4" )
    path( "midas2_output/${sample_id}/bt2_indexes/snps/*" )
   
   //run_snps isn't working because it requires output from run_species-need to add a when statement to ensure run_snps occurs after run_species
   //when: midas2_output/${sample_id}/species/species_profile.tsv or species_id
    script:
   """
    midas2 run_snps \
      --sample_name ${sample_id} \
      -1 ${reads[0]} \
      -2 ${reads[1]} \
      --midasdb_name uhgg \
      --midasdb_dir my_midasdb_uhgg \
      --select_by median_marker_coverage,unique_fraction_covered \
      --select_threshold=2,0.5 \
      --num_cores 8 \
      midas2_output
    """
    when:
    file("MIDAS2/midas2_output/${sample_id}/species/species_profile.tsv").isFile() //this did not work

    stub:
    """
    mkdir midas2_output
    mkdir midas2_output/${sample_id}
    mkdir midas2_output/${sample_id}/bt2_indexes
    mkdir midas2_output/${sample_id}/bt2_indexes/snps
    touch midas2_output/${sample_id}/bt2_indexes/snps/stub
    mkdir midas2_output/${sample_id}/snps
    touch midas2_output/${sample_id}/snps/log.txt
    touch midas2_output/${sample_id}/snps/snps_summary.tsv
    touch midas2_output/${Sample_id}/snps/stub.snps.tsv.lz4
    """
}
/*
 * Parse MIDAS2 output to get readable list of potential species in sample. 
 */
process MIDAS2_PARSE {
    tag{"MIDAS2_PARSE "} //need to include variables for species and snps MIDAS2 outputs
    label 'process_low'

    publishDir("${params.outdir}/midas2_output", mode: 'copy') 

    input:
    tuple val( sample_id ), path( reads )

    output:
    path( "" ), emit: species_id_list

    script:
    """
    
    """

    stub:
    """
   
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
