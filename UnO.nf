/*
========================================================================================
   UnO Nextflow Workflow
========================================================================================
   Github   : github.com/uel3/UnO_nf
   Contact  :
   Commands to run: 
   $module load nextflow/22.10.6
   $nextflow run UnO.nf -profile conda, sge
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Pipeline Input parameters

params.outdir = 'results_multi'
//params.reads = "$HOME/coal_reads/subset_3/*_{1,2}.fastq.gz" //using a subet of 3 reads to test multiple reads on the pipeline-10 sets was too much 
params.reads = "$HOME/UnO_nf/test/*_{1,2}*.fastq.gz" //test data

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
 reads_ch = Channel
            .fromFilePairs(params.reads)
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
            .map { row ->
                        def sample = [:]
                        sample.id  = row[0]
                        sample.group = "reads"  // Ussing 0 for thr group to reflect nf-mag structure 
                        return [ sample, row[1] ]
                }

//adapter_ch = Channel.fromPath( params.adapter, checkIfExists: true )

/*
========================================================================================
   Import Modules
========================================================================================
*/
//include { MIDAS2 as MIDAS2 } from '.../UnO_nf' //need to incorporate MIDAS2.nf into pipeline 
/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {
    raw_reads_ch = FASTQC_RAW( reads_ch )
    //checkpoint
    trimmed_reads_ch = TRIMMOMATIC( reads_ch )
        grouped_reads_ch = TRIMMOMATIC.out.trimmed_reads //this would be the trimmed_reads_ch
        .map { sample, reads -> [sample.group, sample, reads] }
        .groupTuple(by: 0)
        .map { group, samples, reads ->
            def groupedSample = [:]
            groupedSample.id = "grouped_$group" //with this sample.id for reads in ch_short_reads_grouped is the "paired_reads"
            groupedSample.group = group //with this the sample.group for reads in ch_short_reads_grouped is "reads" -I can call this?
            def reads1 = reads.collect { it[0] }
            def reads2 = reads.collect { it[1] }
            [groupedSample, reads1, reads2]
    }
    FASTQC_TRIMMED( trimmed_reads_ch.trimmed_reads )
    megahit_assembly_ch = MEGAHIT( grouped_reads_ch )
    //checkpoint
    //MIDAS2_TRIMMED ( TRIMMOMATIC.out.trimmed_reads )
    bt2_index_ch = BOWTIE2_INDEX( megahit_assembly_ch.megahit_contigs )
    trimmed_reads_out = TRIMMOMATIC.out.trimmed_reads
    bt2_reads_input = trimmed_reads_out.map{ sample, reads -> [sample.group, sample, reads]}
    bt2_map_input_ch = BOWTIE2_INDEX.out.bowtie2_index //Formatting added to combine BAMs from all reads into JGISUMMARIZE process based on groups--maintains std formatting
        .map {sample, assembly, index -> [ sample.group, sample, assembly, index] }
        .combine(bt2_reads_input, by:0)
        .map {group, assembly_sample, assembly, index, sample, reads -> [assembly_sample, assembly, index, sample, reads] }
    mapped_reads_ch = BOWTIE2_MAP_READS( bt2_map_input_ch )
        .groupTuple(by:0)
        .map { assembly_sample, assembly, bams, bais -> [ assembly_sample, assembly.sort()[0], bams, bais ] }
    text_file_ch = grouped_reads_ch
        .map { sample, reads1, reads2 -> 
        reads1.join("\n") + "\n" + reads2.join("\n") 
        }
        .collectFile(name: 'grouped_reads.txt') //may need to use the mag functionality of converting maxbin2 depth to abundance file instead of this...
    max_bin_ch = MAXBIN2_BIN( megahit_assembly_ch.megahit_contigs, text_file_ch )
    jgi_mapped_reads_ch = mapped_reads_ch.map{ assembly_sample, assembly, bams, bais -> [assembly_sample, bams, bais]}
    bam_contig_depth_ch = METABAT2_JGISUMMARIZECONTIGDEPTHS( jgi_mapped_reads_ch )
    ch_metabat_depths = METABAT2_JGISUMMARIZECONTIGDEPTHS.out.bam_contig_depth
        .map { assembly_sample, depths ->
            def meta_new = assembly_sample + [binner: 'MetaBAT2']
            [ meta_new, depths ]
        }
    consolidated_metabat_input = mapped_reads_ch
        .map { assembly_sample, assembly, bams, bais ->
            def meta_new = assembly_sample + [binner: 'MetaBAT2']
            [meta_new, assembly, bams, bais]
        }
        .join( ch_metabat_depths, by: 0 )
        .map { assembly_sample, assembly, bam, bai, depths ->
            [ assembly_sample, assembly, depths ]}
    metabat_bins_ch = METABAT2_BIN( consolidated_metabat_input ) //formatting for consolidating BAMS to generate depth and then Bin based on those combinded depths 
    metaquast_ch = METAQUAST_EVAL( megahit_assembly_ch.megahit_contigs )
    //checkpoint
    contig2bin_tsv_ch = DASTOOL_CONTIG2BIN( max_bin_ch.binned_fastas, metabat_bins_ch.binned_fastas)
    refined_dastool_bins_ch = DASTOOL_BINNING(contig2bin_tsv_ch.maxbin2_fastatocontig2bin, contig2bin_tsv_ch.metabat2_fastatocontig2bin, megahit_assembly_ch.megahit_contigs)
    refined_bin_channel_check = DASTOOL_BINNING.out.bins
    bin_evaluation_ch = CHECKM_REFINED(refined_dastool_bins_ch.bins)
    //checkpoint
    prodigal_bin_input_ch = DASTOOL_BINNING.out.bins
        .map {sample, bins -> [bins]}
    flatten_refined_dastool_bins_ch = prodigal_bin_input_ch.flatten() //this flattens the channel so each bin file can be seperately processed by prodigal 
    prodigal_gene_prediction_ch = PRODIGAL_ANON( flatten_refined_dastool_bins_ch )
    //gene_annotation = GENEANNOTATIONTOOL( )
    //taxonomic_classification = SOMETAXTOOLS( )
    //mag_abundace_estimation =MAGABUNDANCETOOL( )
    concatenated_bins_ch = CONCAT_REFINED_BINS( refined_dastool_bins_ch.bins )
    concatenated_refined_bins_index_ch = BT2_INDEX_CONCAT_BIN_FILE( concatenated_bins_ch.concatenated_bins_file )
    concat_bin_index_ch = BT2_INDEX_CONCAT_BIN_FILE.out.bowtie2_concat_bin_index
    //concat_bin_index_ch.view() //shows sample and index
    read_to_bin_mapping_input_ch = concat_bin_index_ch
        .map { sample, index -> [sample.group, index]}
        .combine(bt2_reads_input, by:0)
        .map {group, index, sample, reads -> [index, sample, reads]}
    //read_to_bin_mapping_input_ch.view() //shows index, sample, reads
    reads_mapped_to_bins_ch = MAP_READS_CONCAT_BINS( read_to_bin_mapping_input_ch )
    //reads_mapped_to_bins_ch.view()
    //read_coverage_bins_ch = JGISUMMARIZE_BINS_DEPTH( reads_mapped_to_bins_ch.mapped_reads_bam )
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
    label 'UnO'
    publishDir("${params.outdir}/fastqc_raw", mode: 'copy')
    
    input:
    tuple val( sample ), path( reads )

    output:
    path( "*_fastqc*" ), emit: raw_reads

    script:
    """
    fastqc ${reads}
    """

    stub:
    """
    touch ${sample.id}_1_fastqc.html
    touch ${sample.id}_1_fastqc.gz
    touch ${sample.id}_2_fastqc.html
    touch ${sample.id}_2_fastqc.gz
    """
}

/*
 * Trimming raw fastq reads.
 */
process TRIMMOMATIC {
    tag{"TRIMMOMATIC ${reads}"}
    label 'UnO'
    publishDir("${params.outdir}/trimmed_reads", mode: 'copy')
  
    input:
    tuple val( sample ), path( reads )

    output:
    tuple val( sample ), path( "*.trimmed.fq.gz" ), emit: trimmed_reads
    tuple val( sample ), path ("*.unpaired.fq.gz"), emit: unpaired


    script:
     def prefix = "${sample.id}"
     """
     trimmomatic PE -threads 10 -phred33 ${reads} ${prefix}_1.trimmed.fq.gz ${prefix}_1.unpaired.fq.gz ${prefix}_2.trimmed.fq.gz ${prefix}_2.unpaired.fq.gz ILLUMINACLIP:TrueSeq3-PE.fa:20:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
     """

    stub:
     """
     touch ${sample.id}_1.trimmed.fq.gz
     touch ${sample.id}_2.trimmed.fq.gz
     touch ${sample.id}_1.unpaired.fq.gz
     touch ${sample.id}_2.unpaired.fq.gz
     """

}

/*
 * QC trimmed fastq reads.
 */
process FASTQC_TRIMMED {
    tag{"FASTQC ${reads_trimmed}"}
    label 'UnO'
    publishDir("${params.outdir}/fastqc_trimmed", mode: 'copy')
    
    input:
    tuple val( sample ), path( reads_trimmed )

    output:
    path( "*_fastqc*" )

    script:
    """
    fastqc ${reads_trimmed}
    """

    stub:
    """
    touch ${sample.id}_1_trimmed_fastqc.html
    touch ${sample.id}_2_trimmed_fastqc.html
    touch ${sample.id}_1_trimmed.fastqc.gz
    touch ${sample.id}_2_trimmed.fastqc.gz
    """
    
}

/*
 * Metagenomic co-assembly of reads with megahit.
 */
process MEGAHIT {
    tag{"MEGAHIT ${reads1} ${reads2}"}
    label 'UnO'
    publishDir("${params.outdir}/megahit_out", mode: 'copy')
    
    input:
    tuple val( sample ), path( reads1 ), path( reads2 )

    output:
    tuple val( sample ), path("megahit_out/*.contigs.fa")                            , emit: megahit_contigs
    tuple val( sample ), path("megahit_out/intermediate_contigs/k*.contigs.fa")      , emit: megahit_k_contigs
    tuple val( sample ), path("megahit_out/intermediate_contigs/k*.addi.fa")         , emit: megahit_addi_contigs
    tuple val( sample ), path("megahit_out/intermediate_contigs/k*.local.fa")        , emit: megahit_local_contigs
    tuple val( sample ), path("megahit_out/intermediate_contigs/k*.final.contigs.fa"), emit: megahit_kfinal_contigs


    script:
    def prefix = "${sample.group}"
    def input = "-1 \"" + reads1.join(",") + "\" -2 \"" + reads2.join(",") + "\""
    """
    megahit --presets meta-large $input -t 12 --mem-flag 2 --out-prefix "MEGAHIT-${prefix}"
    """
    //adding MEGAHIT to denote assembler used per mag
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
    label 'UnO'
    publishDir("${params.outdir}/bowtie2_out", mode: 'copy')

    input:
    tuple val( sample ), path( assembly )

    output:
    tuple val( sample ), path( assembly ), path ( "${sample.group}_index*.bt2" ), emit: bowtie2_index //adding path(assembly) per mag 

    script:
    def prefix = "${sample.group}"
    """
    bowtie2-build ${assembly} ${prefix}_index
    touch ${prefix}_index
    """

    stub:
    """
    mkdir bowtie2_out
    touch ${sample.group}_index.1.bt2
    touch ${sample.group}_index.2.bt2
    touch ${sample.group}_index.3.bt2
    touch ${sample.group}_index.4.bt2
    touch ${sample.group}_index.rev.1.bt2
    touch ${sample.group}_index.rev.2.bt2
    touch ${sample.group}_index
    """
}
/* 
* Mapping reads to bowtie2 index to evaluate coverage and for downstreaming binning tools.
*/
process BOWTIE2_MAP_READS {
    tag "BOWTIE2_MAP_READS ${index} ${reads_trimmed}" //"$assembler-$name"
    label 'UnO'
    publishDir ("${params.outdir}/bowtie2_out/mapped", mode: 'copy') //Assembly/${assembler}/${name}_QC", mode: params.publish_dir_mode,
        //saveAs: {filename -> filename.indexOf(".bowtie2.log") > 0 ? filename : null}

    input:
    tuple val( assembly_sample ), path( assembly ), path( index ), val( reads_sample ), path( reads_trimmed ) //things to consider--this is using unformatted output from trimmomatic. It contains reads information-unlike the grouping steps -look into prefix set up--
    //val   save_unaligned
    //val   sort_bam
    
    
    output:
    tuple val( assembly_sample ), path( assembly ), path( "${reads_sample.id}_sorted.bam" ), path( "${reads_sample.id}_sorted.bam.bai"), emit: aligned_bam
    //path( "${sample_id}.bowtie2.log" ), emit: align_log

    script:
    //def idx = index[0].getBaseName(2)
    def prefix = "${reads_sample.id}"
    def input = "-1 \"${reads_trimmed[0]}\" -2 \"${reads_trimmed[1]}\""
    //need to look into getting these formatted accordingly...how does this work with multiple sets of reads...the sample.ids are for each reads sample.group is for the set of group reads 
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2l" -o -name "*.rev.1.bt2" | sed 's/.rev.1.bt2l//' | sed 's/.rev.1.bt2//'`
    bowtie2 -p 8 -x \$INDEX -q $input --no-unal |samtools view -@ 2 -b -S -h | samtools sort -o ${prefix}_sorted.bam 
    samtools index ${prefix}_sorted.bam
    """
    //required the -h flag to make this work into view/sort commands
    //needed to dealre the correct output in declared output
    // will need this to accomodate multiple reads 
    stub:
    """
    mkdir mapped
    touch ${sample.id}_sorted.bam
    touch ${sample.id}_sorted.bam.bai
    """
}
/*
 * MaxBin2 binning of assembled contigs
 */
process MAXBIN2_BIN {
    tag "MAXBIN2_BIN ${assembly} ${readstext}"
    label 'UnO'
    publishDir ("${params.outdir}/MaxBin2", mode: 'copy')

    input:
    tuple val( sample ), path( assembly )
    file(readstext) 

    output:
    tuple val( sample ), path( "*.fasta" )   , emit: binned_fastas
    tuple val( sample ), path( "*.summary" ) , emit: summary
    tuple val( sample ), path( "*.log" )     , emit: log
    tuple val( sample ), path( "*.marker" )  , emit: marker_counts
    tuple val( sample ), path( "*.noclass" ) , emit: unbinned_fasta
    tuple val( sample ), path( "*.tooshort" ), emit: tooshort_fasta
    tuple val( sample ), path( "*.abund*" )  , emit: abundance, optional: true
    tuple val( sample ), path( "*_bin.tar.gz" ) , emit: marker_bins , optional: true
    
    script:
    """
    run_MaxBin.pl -thread 8 -contig ${assembly} -out MaxBin2_${assembly.simpleName} -reads_list ${readstext}
    """
    //getting an inital error of run_MaxBin.pl: command not found-looking into it-mamba install maxbin2 
    //mag has additonal code to zip the fasta bins--might consider to cut down on space
    stub:
    """
    mkdir MaxBin2
    touch MaxBin2.fasta
    touch MaxBin2.summary
    touch MaxBin2.log
    touch MaxBin2.marker
    touch MaxBin2.noclass
    touch MaxBin2.tooshort
    touch MaxBin2.abundstub
    touch MaxBin2_bin.tar.gz
    """
}
/*
 * JGIsummarizeBAMcontigdepths script within metabat2 
 */
 process METABAT2_JGISUMMARIZECONTIGDEPTHS {
    tag "METABAT2_JGISUMMARIZECONTIGDEPTHS ${aligned_bam}"
    label 'UnO'
    publishDir ("${params.outdir}/MetaBat2", mode: 'copy')
    ignoreEmpty = true
    input:
    tuple val( assembly_sample ), path( aligned_bam ), path( aligned_bai )
    
    output:
    tuple val( assembly_sample ), path( "*_depth.txt" ) , emit: bam_contig_depth

    script:
    def prefix = "${assembly_sample.id}"
    """
    jgi_summarize_bam_contig_depths --outputDepth ${prefix}_depth.txt ${aligned_bam}
    """

    stub:
    """
    mkdir MetaBat2
    touch stub.depth.txt
    """
 }
 /*
 * MetaBat2 binning of assembled contig
 */
process METABAT2_BIN {
    tag "METABAT2_BIN ${assembly} ${depth}"
    label 'UnO'
    publishDir ("${params.outdir}/MetaBat2", mode: 'copy')
    ignoreEmpty = true //this makes MetaBat2 wait to receive all input files from JGISUMMARIZE
    input:
    tuple val( assembly_sample ), path( assembly ), path( depth )
    
    output:
    //path( "${assembly}.depth.txt" )                                             , emit: metabat2_depth
    //path( "${assembly}.unbinned.fa" )                                           , emit: unbinned
    //path( "${assembly}.paired.txt" )                                            , emit: summary
    tuple val( assembly_sample ), path( "*.fa" )                                  , emit: binned_fastas //changing the struture of this to allow me to call the directory later 
        
    script:
    """
    metabat2 -i ${assembly} -a ${depth} -o MetaBat2_${assembly.baseName}
    """
    //mag has additonal code to zip the fasta bins--might consider to cut down on space
    //formatted bin outdir correctly for legibility 
    stub:
    """
    touch ${assembly_sample}.paired.txt
    touch ${assembly_sample}.unbinned.fa
    touch ${assembly_sample}.fa
    """
}
//adding quast to environment via mamba 
/*
*Metaquast assembly evaluation
*/
process METAQUAST_EVAL {
    tag "METAQUAST_EVAL ${assembly}"
    label 'UnO'
    publishDir ("${params.outdir}/Metaquast_out", mode: 'copy')

    input:
    tuple val( sample ), path( assembly )

    output:
    path( "${assembly.simpleName}/*" )                   , emit: quast_qc
    //path( "${assembly.simpleName}/icarus.html" )                    , emit: icarus
    //path( "${assembly.simpleName}/metaquast.log" )                  , emit: metaquast_log
    //path( "${assembly.simpleName}/report.html" )                    , emit: metaquast_report
    //path( "${assembly.simpleName}/combined_reference/*" )           , emit: combined_references
    //path( "${assembly.simpleName}/icarus_viewers/*" )               , emit: icarus_viewers
    //path( "${assembly.simpleName}/krona_charts" )                   , emit: krona_charts
    //path( "${assembly.simpleName}/not_aligned/*" )                  , emit: metaquast_unaligned
    //path( "${assembly.simpleName}/quast_downloaded_references*" )   , emit: quast_references
    //path( "${assembly.simpleName}/runs_per_reference/*" )           , emit: runs_reference
    //path( "${assembly.simpleName}/summary*" )                       , emit: summary
   
    script:
    """
    metaquast.py -o ${assembly.simpleName} -t 8 ${assembly}
    """

    stub:
    """
    mkdir Metaquast_out
    mkdir ${assembly.simpleName}
    touch ${assembly.simpleName}/*
    """
    //touch ${assembly.simpleName}/metaquast.log
    //touch ${assembly.simpleName}/report.html
    //mkdir combined_reference
    //touch ${assembly.simpleName}/combined_reference/*
    //mkdir icarus_viewers
    //touch ${assembly.simpleName}/icarus_viewers/*
    //mkdir krona_charts
    //touch ${assembly.simpleName}/krona_charts
    //mkdir not_aligned
    //touch ${assembly.simpleName}/not_aligned/*
    //mkdir quast_downloaded_references
    //touch ${assembly.simpleName}/quast_downloaded_references*
    //mkdir runs_per_reference
    //touch ${assembly.simpleName}/runs_per_reference/*
    //mkdir summary 
    //touch ${assembly.simpleName}/summary*

}
/*
*Converting contigs bins into TSV files for DAStool
*/
process DASTOOL_CONTIG2BIN { //needed to add a conda profile for das_tool environment-set for this process specifically 
    tag "DASTOOL_CONTIG2BIN ${maxbin_fasta} ${metabat_fasta}"
    label 'dastool'

    publishDir ("${params.outdir}/DASTool_out", mode: 'copy')
    
    //needed to add a conda profile for das_tool environment-set for this process specifically  
    input:
    tuple val( sample ), path( maxbin_fasta )
    tuple val( assembly_sample ), path( metabat_fasta )


    output:
    path("*maxbin.contigs2bin.tsv"), emit: maxbin2_fastatocontig2bin
    path("*metabat.contigs2bin.tsv"), emit: metabat2_fastatocontig2bin
    
    script: //shortened output names to see if output is not longer empty -didn't work adding . as the input worked will need to include better names than just maxbin/metabat.contigs2bin.tsv but works for now 
    """
    Fasta_to_Contig2Bin.sh -i . -e fasta > maxbin.contigs2bin.tsv 
    Fasta_to_Contig2Bin.sh -i . -e fa > metabat.contigs2bin.tsv

    """

    stub:
    """
    mkdir DASTool_out
    touch DASTool_out/stub_maxbin.contigs2bin.tsv
    touch DASTool_out/stub_metabat.contigs2bin.tsv
    """

}
/*
*Refining MaxBin2 and MetaBat2 bins with DASTool 
*/
process DASTOOL_BINNING { //can provide link to conda environemnt with yml file. This may be helpful for dastool and MIDAS2 
    tag "DASTOOL_BINNING ${maxbin_tsv} ${metabat_tsv} ${assembly}"
    label 'dastool'
    publishDir ("${params.outdir}/DASTool_out", mode: 'copy')

    input:
    path( maxbin_tsv )
    path( metabat_tsv )
    tuple val( sample ), path( assembly )

    output: //updating the output based on DASTool NF module output-this did nothing to change the error-per DASTool github this is an error with command line syntax
    tuple val( sample ), path("*.log")                                      , emit: log
    tuple val( sample ), path("*_summary.tsv")              , optional: true, emit: summary
    tuple val( sample ), path("*_DASTool_contig2bin.tsv")   , optional: true, emit: contig2bin
    tuple val( sample ), path("*.eval")                     , optional: true, emit: eval
    tuple val( sample ), path("*_DASTool_bins/*.fa")        , optional: true, emit: bins
    tuple val( sample ), path("*.pdf")                      , optional: true, emit: pdfs
    tuple val( sample ), path("*.candidates.faa")           , optional: true, emit: fasta_proteins
    tuple val( sample ), path("*.faa")                      , optional: true, emit: candidates_faa
    tuple val( sample ), path("*.archaea.scg")              , optional: true, emit: fasta_archaea_scg
    tuple val( sample ), path("*.bacteria.scg")             , optional: true, emit: fasta_bacteria_scg
    tuple val( sample ), path("*.b6")                       , optional: true, emit: b6
    tuple val( sample ), path("*.seqlength")                , optional: true, emit: seqlength
    
    script:
    """
    DAS_Tool -i ${maxbin_tsv},${metabat_tsv}\
    -c ${assembly}\
    -o ${assembly.simpleName}_refined_bins\
    -l MaxBin2,MetaBat2\
    --write_bins\
    --debug
    """

    stub:
    """
    mkdir DASTool_out
    touch DASTool_out/${assembly.simpleName}_refined_bins_DASTool_contig2bin.tsv
    touch DASTool_out/${assembly.simpleName}_refined_bins_DASTool_summary.tsv
    touch DASTool_out/${assembly.simpleName}_refined_bins_DASTool.log
    touch DASTool_out/${assembly.simpleName}_refined_bins_proteins.faa
    touch DASTool_out/${assembly.simpleName}_refined_bins_proteins.faa.all.b6
    touch DASTool_out/${assembly.simpleName}_refined_bins_proteins.faa.archaea.scg
    touch DASTool_out/${assembly.simpleName}_refined_bins_proteins.faa.bacteria.scg
    touch DASTool_out/${assembly.simpleName}_refined_bins_proteins.faa.findSCG.b6
    touch DASTool_out/${assembly.simpleName}_refined_bins_proteins.faa.scg.candidates.faa
    touch DASTool_out/${assembly.simpleName}_refined_bins.seqlength
    mkdir DASTool_out/${assembly.simpleName}_refined_bins_DASTool_bins
    touch DASTool_out/${assembly.simpleName}_refined_bins_DASTool_bins/Binner.fa
    """
}
/*
*Evaluating Bin quality with CheckM 
*/
process CHECKM_REFINED { /*adding checkM to the environment downgrades some packages The following NEW packages will be INSTALLED:
  htslib                                    1.18-h81da01d_0 --> 1.17-h81da01d_2
  libdeflate                                1.19-hd590300_0 --> 1.18-h0b41bf4_0
  libtiff                                  4.6.0-h29866fb_1 --> 4.6.0-h8b53f26_0*/
    tag "DASTOOL_BINNING ${refined_bins}"
    label 'UnO'
    publishDir ("${params.outdir}", mode: 'copy') //need to add a way to include naming information for labeling purposes
    
    input:
    tuple val( sample ), path( refined_bins )
    
    output:
    path("CheckM/bins/*")                                  ,  emit: checkm_bins
    path("CheckM/checkm.log")                              ,  optional: true, emit: checkm_log
    path("CheckM/lineage.ms")                              ,  emit: marker_file
    path("CheckM/*.tsv")                                   ,  emit: checkm_output_qa 
    path("CheckM/storage/*.tsv")                           ,  emit: checkm_stats
    path("CheckM/storage/aai_qa/*")                        , optional: true, emit: aai_qa
    path("CheckM/storage/*_info.pkl.gz")                   , optional: true, emit: checkm_info 
    path("CheckM/storage/tree/*")                          , optional: true, emit: tree

    script:
    """
    checkm lineage_wf -t 10 -x fa -f CheckM/checkm_qa.tsv --tab_table . CheckM
    """
    
    stub:
    """
    mkdir CheckM
    touch CheckM/lineage.ms
    mkdir CheckM/bins
    touch CheckM/bins/stub.bin
    mkdir CheckM/storage
    touch CheckM/storage/bin_stats.analyze.tsv
    touch CheckM/storage/bin_stats_ext.tsv
    touch CheckM/storage/bin_stats_tree.tsv
    touch CheckM/storage/marker_gene_stats.tsv
    touch CheckM/storage/stub_info.pkl.gz
    mkdir CheckM/storage/aai_qa
    touch CheckM/storage/aai_qa/stub
    """
}
//there should be a step that only proceeds with Prodigal Gene prediction for HQ or MQ bins
/*
Predicting genes from refined DASTool Bins with Prodigal
*/
process PRODIGAL_ANON{
    tag "PRODIGAL_ANON ${refined_bins}" //prodigal is already in UnO.yml from previous package
    label 'UnO'
    publishDir ("${params.outdir}/Prodigal", mode: 'copy')

    input:
    path( refined_bins )
    
    output:
    path("${refined_bins.baseName}.gff"),                 emit: gene_annotations //might need to change these to ${sample} instead of the refined.baseName
    path("${refined_bins.baseName}.fna"),                 emit: nucleotide_fasta
    path("${refined_bins.baseName}.faa"),                 emit: amino_acid_fasta
    path("${refined_bins.baseName}_all.txt"),             emit: all_gene_annotations
    
    script:
    """
    prodigal \\
        -i $refined_bins \\
        -f gff \\
        -d ${refined_bins.baseName}.fna \\
        -o ${refined_bins.baseName}.gff \\
        -a ${refined_bins.baseName}.faa \\
        -s ${refined_bins.baseName}_all.txt
    """
    stub:
    """
    """
}
/*
Creating and Bowtie2 index of refined DASTool bins to map trimmmed to
*/
process BOWTIE2_INDEX_BINS{
    tag "BOWTIE2_INDEX_BINS ${refined_bins}"
    label 'UnO'
    publishDir ("${params.outdir}/bowtie2_out/bin_index", mode: 'copy')

    input:
    path( refined_bins )

    output:
    path( "${refined_bins.baseName}_index*.bt2" ), emit: bowtie2_bin_index 

    script:
    def prefix = "${refined_bins.baseName}"
    """
    bowtie2-build ${refined_bins} ${prefix}_index
    touch ${prefix}_index
    """

    stub:
    """
    mkdir bowtie2_out
    mkdir bowtie2_out/bin_index
    touch ${sample}_index.1.bt2
    touch ${sample}_index.2.bt2
    touch ${sample}_index.3.bt2
    touch ${sample}_index.4.bt2
    touch ${sample}_index.rev.1.bt2
    touch ${sample}_index.rev.2.bt2
    touch ${sample}_index
    """
}
/* 
Mapping reads to Bins with Bowtie2
*/
process BOWTIE2_MAP_BINS{
    tag "BOWTIE2_MAP_BINS ${index} ${reads_trimmed}" //"$assembler-$name"
    label 'UnO'
    publishDir ("${params.outdir}/bowtie2_out/mapped_bins", mode: 'copy') 

    input:
    path( index )
    tuple val( sample ), path( reads_trimmed )  
    //val   save_unaligned
    //val   sort_bam
    
    output:
    path( "${sample.id}_sorted.bam" ), emit: aligned_bam_to_bins
    path( "${sample.id}_sorted.bam.bai"), emit: aligned_bam_to_bins_index
    //path( "${sample_id}.bowtie2.log" ), emit: align_log

    script:
    def prefix = "${sample.id}"
    def input = "-1 \"${reads_trimmed[0]}\" -2 \"${reads_trimmed[1]}\""
    //need to look into getting these formatted accordingly...how does this work with multiple sets of reads...the sample.ids are for each reads sample.group is for the set of group reads 
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2l" -o -name "*.rev.1.bt2" | sed 's/.rev.1.bt2l//' | sed 's/.rev.1.bt2//'`
    bowtie2 -p 8 -x \$INDEX -q $input --no-unal |samtools view -@ 2 -b -S -h | samtools sort -o ${prefix}_sorted.bam 
    samtools index ${prefix}_sorted.bam
    """
    //required the -h flag to make this work into view/sort commands
    //needed to dealre the correct output in declared output
    // will need this to accomodate multiple reads 
    stub:
    """
    mkdir mapped_bins
    touch ${sample.id}_sorted.bam
    touch ${sample.id}_sorted.bam.bai
    """
}

/* 
Including alternate read mapping approach for co-assembled + refined bins
*/
process CONCAT_REFINED_BINS{
    tag "CONCAT_REFINED_BINS ${refined_bins}"
    label 'UnO'
    publishDir ("${params.outdir}/DASTool_out", mode: 'copy')

    input:
    tuple val( sample ), path( refined_bins )

    output:
    tuple val( sample ), path("DAStool_bins_concatenated.fna"), emit: concatenated_bins_file

    script:
    """
    cat ${refined_bins} > DAStool_bins_concatenated.fna
    """
    stub:
    """
    mkdir DASTool_out
    touch DAStool_bins_concatenated.fna
    """

}
/*
Generate index file for concatenated bin file
*/
process BT2_INDEX_CONCAT_BIN_FILE{
    tag "BT2_INDEX_CONCAT_BIN_FILE ${concat_bins_file}" // need to check inputs and outputs
    label 'UnO'
    publishDir ("${params.outdir}/bowtie2_out/concat_bin_index", mode: 'copy')
    
    input:
    tuple val(sample), path( concat_bins_file )

    output:
    tuple val( sample), path( "${concat_bins_file.baseName}_index*.bt2" ), emit: bowtie2_concat_bin_index 

    script:
    def prefix = "${concat_bins_file.baseName}"
    """
    bowtie2-build ${concat_bins_file} ${prefix}_index
    touch ${prefix}_index
    """
    stub:
    """
    mkdir bowtie2_out
    mkdir concat_bin_index
    touch ${concat_bins_file.baseName}_index.1.bt2
    touch ${concat_bins_file.baseName}_index.2.bt2
    touch ${concat_bins_file.baseName}_index.3.bt2
    touch ${concat_bins_file.baseName}_index.4.bt2
    touch ${concat_bins_file.baseName}_index.rev.1.bt2
    touch ${concat_bins_file.baseName}_index.rev.2.bt2
    touch ${concat_bins_file.baseName}_index
    """
}
/*
Mapping reads to bin with bowtie index 
*/
process MAP_READS_CONCAT_BINS{
    tag "MAP_READS_CONCAT_BINS ${index} ${reads_trimmed}" //need to double check on inputs and outputs 
    label 'UnO'
    publishDir ("${params.outdir}/bowtie2_out/mapped_bins", mode: 'copy') 

    input:
    tuple path(index), val(sample), path(reads_trimmed)

    output:
    tuple val( sample ), path( "${sample.id}_sorted.bam" ), emit: mapped_reads_bam
    tuple val( sample ), path( "${sample.id}_sorted.bam.bai"), emit: mapped_reads_bam_index
    //path( "${sample_id}.bowtie2.log" ), emit: align_log

    script:
    def prefix = "${sample.id}"
    def input = "-1 \"${reads_trimmed[0]}\" -2 \"${reads_trimmed[1]}\""
    //need to look into getting these formatted accordingly...how does this work with multiple sets of reads...the sample.ids are for each reads sample.group is for the set of group reads 
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2l" -o -name "*.rev.1.bt2" | sed 's/.rev.1.bt2l//' | sed 's/.rev.1.bt2//'`
    bowtie2 -p 8 -x \$INDEX -q $input --no-unal |samtools view -@ 2 -b -S -h | samtools sort -o ${prefix}_sorted.bam 
    samtools index ${prefix}_sorted.bam
    """

    stub:
    """
    mkdir mapped_bins
    touch ${sample.id}_sorted.bam
    touch ${sample.id}_sorted.bam.bai
    """
}
/*
Using MetaBat2 to generate covereage information from reads mapped to index
*/
process JGISUMMARIZE_BINS_DEPTH{
    tag "JGISUMMARIZE_BINS_DEPTH ${mapped_reads_bam}" //need to check that these work together 
    label 'UnO'
    publishDir ("${params.outdir}/Read_Coverage", mode: 'copy')
    ignoreEmpty = true
    
    input:
    tuple val( sample ), path( mapped_reads_bam ) //all bam files need to be put into 1 depth file
    
    output:
    tuple val( sample ), path( "*_bins_cov_table.txt" ) , emit: bam_contig_depth_mapped_reads //this sample value will need to be the name of the refined bins.fna file

    script:
    def prefix = "DasTool_refined_bins"
    """
    jgi_summarize_bam_contig_depths --outputDepth ${prefix}_bins_cov_table.txt ${mapped_reads_bam}
    """

    stub:
    """
    mkdir Read_Coverage
    touch ${prefix}_bins_cov_table.txt
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
