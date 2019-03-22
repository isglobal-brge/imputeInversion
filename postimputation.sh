#!/bin/bash

#######################################
# Post process imputation data
# Arguments:
#	  $1: Chromosome to be phased
#   $2: folder path
#   $3: base name
#   $4: prefix to name inversions
#   $5: number of CPUs 
#   $6: bcftools reference
# Returns:
#   None
#######################################

postimpute(){
  chr=$1
  dir=$2
  base=$3
  prefix=$4
  cpus=$5
  ref=$6

  if [ ! -d $dir/postimputation_files/$prefix ]; then
    mkdir $dir/postimputation_files/$prefix # Create file with the chromosome number to store the files separated by chromosomes
  fi
  
  if [ ! -d $dir/${base}_imputed_files/$prefix ]; then
    mkdir $dir/${base}_imputed_files/$prefix
  fi
  
  if [[ $chr == X ]] # If it is the chromosome X
  then
    # Compress the male imputed file with bgzip to avoid problems using bcftools afterwards
    zcat $dir/pimputed_files/$prefix/${prefix}_${base}_males_imputed.dose.vcf.gz | bgzip -c > $dir/postimputation_files/$prefix/${prefix}_${base}_males_imputed_bgzip.vcf.gz && tabix $dir/postimputation_files/$prefix/${prefix}_${base}_males_imputed_bgzip.vcf.gz
    # Store header in a text file to modify it
    bcftools view -h $dir/postimputation_files/$prefix/${prefix}_${base}_males_imputed_bgzip.vcf.gz > $dir/postimputation_files/hdrmales.txt
    # Add lines to the header conatining information about ID=GENOTYPED (avoid problems in the next steps)
    sed -i '/^##FORMAT=<ID=GT/i ##FILTER=<ID=GENOTYPED,Description="Marker was genotyped AND imputed">\n##FILTER=<ID=GENOTYPED_ONLY,Description="Marker was genotyped but NOT imputed">' $dir/postimputation_files/hdrmales.txt
    # Create a temporal file with the new header and the data from the imputed file
    bcftools reheader -h $dir/postimputation_files/hdrmales.txt --output $dir/postimputation_files/$prefix/${prefix}_${base}_males_hdr_imputed_bgzip.vcf.gz $dir/postimputation_files/$prefix/${prefix}_${base}_males_imputed_bgzip.vcf.gz
    tabix -p vcf $dir/postimputation_files/$prefix/${prefix}_${base}_males_hdr_imputed_bgzip.vcf.gz
    # Compress the female imputed file with bgzip to avoid problems using bcftools afterwards
    zcat $dir/pimputed_files/$prefix/${prefix}_${base}_females_imputed.dose.vcf.gz | bgzip -c > $dir/postimputation_files/$prefix/${prefix}_${base}_females_imputed_bgzip.vcf.gz && tabix $dir/postimputation_files/$prefix/${prefix}_${base}_females_imputed_bgzip.vcf.gz
    # Store header in a text file to modify it
    bcftools view -h $dir/postimputation_files/$prefix/${prefix}_${base}_females_imputed_bgzip.vcf.gz > $dir/postimputation_files/hdrfemales.txt
    # Add lines to the header conatining information about ID=GENOTYPED (avoid problems in the next steps)
    sed -i '/^##FORMAT=<ID=GT/i ##FILTER=<ID=GENOTYPED,Description="Marker was genotyped AND imputed">\n##FILTER=<ID=GENOTYPED_ONLY,Description="Marker was genotyped but NOT imputed">' $dir/postimputation_files/hdrfemales.txt
    # Create a temporal file with the new header and the data from the imputed file
    bcftools reheader -h $dir/postimputation_files/hdrfemales.txt --output $dir/postimputation_files/$prefix/${prefix}_${base}_females_hdr_imputed_bgzip.vcf.gz postimputation_files/$prefix/${prefix}_${base}_females_imputed_bgzip.vcf.gz
    tabix -p vcf $dir/postimputation_files/$prefix/${prefix}_${base}_females_hdr_imputed_bgzip.vcf.gz
    # Merge files containing the female and the male individuals 
    bcftools merge --force-samples --output-type z --output $dir/postimputation_files/TEMP1.vcf.gz $dir/postimputation_files/$prefix/${prefix}_${base}_males_hdr_imputed_bgzip.vcf.gz $dir/postimputation_files/$prefix/${prefix}_${base}_females_hdr_imputed_bgzip.vcf.gz
    tabix -p vcf $dir/postimputation_files/TEMP1.vcf.gz
    rm $dir/postimputation_files/TEMP* # Remove temporal files created during the process
    rm $dir/postimputation_files/hdr* # Remove files containing the headers
    
  else # For the rest of the chromosomes
    #only to be executed in case the file exists
    if [ -f $dir/pimputed_files/$prefix/${prefix}_${base}_imputed.dose.vcf.gz ]
    then
        # Compress the imputed files with bgzip to avoid problems using bcftools afterwards
        zcat $dir/pimputed_files/$prefix/${prefix}_${base}_imputed.dose.vcf.gz | bgzip -c > $dir/postimputation_files/$prefix/${prefix}_${base}_imputed_bgzip.vcf.gz && tabix $dir/postimputation_files/$prefix/${prefix}_${base}_imputed_bgzip.vcf.gz    
        # Store header in a text file to modify it
        bcftools view -h $dir/postimputation_files/$prefix/${prefix}_${base}_imputed_bgzip.vcf.gz > $dir/postimputation_files/hdr.txt
        # Add lines to the header conatining information about ID=GENOTYPED (avoid problems in the next steps) 
        sed -i '/^##FORMAT=<ID=GT/i ##FILTER=<ID=GENOTYPED,Description="Marker was genotyped AND imputed">\n##FILTER=<ID=GENOTYPED_ONLY,Description="Marker was genotyped but NOT imputed">' $dir/postimputation_files/hdr.txt
        # Create a temporal file with the new header and the data from the imputed file
        bcftools reheader -h $dir/postimputation_files/hdr.txt --output $dir/postimputation_files/TEMP1.vcf.gz $dir/postimputation_files/$prefix/${prefix}_${base}_imputed_bgzip.vcf.gz
        tabix -p vcf $dir/postimputation_files/TEMP1.vcf.gz
    fi
  fi
  # Change ids by chr:position:ref:alt prior to set the rsnumber ids
  bcftools annotate --remove ID --set-id +'%CHROM:%POS:%REF:%ALT' --threads $cpus --output-type z --output $dir/postimputation_files/TEMP2.vcf.gz $dir/postimputation_files/TEMP1.vcf.gz
  tabix -p vcf $dir/postimputation_files/TEMP2.vcf.gz
  # Set rsnumber ids in the final file, using as a reference a file containing information about all the SNPs (in hg19, chr:pos:ref:alt and the rsnumber id)
  bcftools annotate --annotations $ref --columns ID --threads $cpus --output-type z --output $dir/${base}_imputed_files/$prefix/${prefix}_${base}_imputed_final.vcf.gz $dir/postimputation_files/TEMP2.vcf.gz
  tabix -p vcf $dir/${base}_imputed_files/$prefix/${prefix}_${base}_imputed_final.vcf.gz
  rm $dir/postimputation_files/TEMP* # Remove temporal files created during the process
  rm $dir/postimputation_files/hdr* # Remove files containing the headers
}


