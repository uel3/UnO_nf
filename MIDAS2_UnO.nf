midas2 run_species \
      --sample_name ${sample_name} \
      -1 /scicomp/home-pure/uel3/coal_reads/${sample_name}_1.fastq \ #this is path to reads
      -2 /scicomp/home-pure/uel3/coal_reads/${sample_name}_2.fastq \
      --midasdb_name uhgg \
      --midasdb_dir my_midasdb_uhgg \ #there also needs to be a MIDAS db name assigned
      --num_cores 8 \
      midas2_output #this needs to be named by user
  done
  
  midas2 run_snps \ #need to have run run_species first
      --sample_name ${sample_name} \
      -1 reads/${sample_name}_1.fastq \ #need path to reads
      -2 reads/${sample_name}_2.fastq \
      --midasdb_name uhgg \
      --midasdb_dir my_midasdb_uhgg \ #need to provide name for db
      --select_by median_marker_coverage,unique_fraction_covered \ # need to optimize which parameters give best results
      --select_threshold=2,0.5 \ #same, what are the thresholds we want
      --num_cores 4 \
      midas2_output #this needs to provided in run_species

MIDAS2_output_parsing
#!/bin/bash
get_species_from_snps(){
pathtosnps="$1"
pathtomidasdbdir="$2"
awk 'NR==FNR{a[$1]=$0;next}$1 in a{print a[$1], $18, $19}' "$pathtosnps"/snps/snps_summary.tsv "$pathtomidasdbdir"/metadata.tsv > Species_SNP_ID.txt
echo "Output Saved to Species_SNP_ID.txt"
}
get_species_from_snps $1 $2 