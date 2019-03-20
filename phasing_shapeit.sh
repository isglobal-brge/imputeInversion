#!/bin/bash

## Stop on an error
set -e

#######################################
# Phase data with shapeit
# Arguments:
#	  $1: Chromosome to be phased
#   $2: folder path
#   $3: base name
#   $4: path for shapeit reference
#   $5: number of CPUs used in shapeit
#   $6: family argument
# Returns:
#   None
#######################################
phase(){
  
  chr=$1
  prefix=$2
  base=$3
  ShapeRef=$4
  cpus=$5
  family=$6
  
  if [ ! -d $prefix/phased_files/$chr ]; then
    mkdir $prefix/phased_files/$chr # Create file with the chromosome number to store the files separated by chromosomes
  fi
  if [[ $chr == X ]] # If it is the chromosome X
    then
      # Run SHAPEIT to phase chromosome X (specifying the '--chrX' flag). SHAPEIT will phase female individuals and will set NA to haploid heterozygous males
      shapeit -B $prefix/prephasing_files/$chr/$chr_${base}_filtered \
            -M $ShapeRef/genetic_map_chrX_nonPAR_combined_b37.txt \
            -O $prefix/phased_files/$chr/$chr_${base}_filtered_phased \
            --chrX \
            --thread $cpus \
  	       -L $prefix/phased_files/${chr}/shapeit_${chr}.log
      # Convert output files from SHAPEIT (.haps) to vcf format      
      shapeit -convert \
            --input-haps $prefix/phased_files/${chr}/${chr}_${base}_filtered_phased \
            --output-vcf $prefix/phased_files/${chr}/${chr}_${base}_filtered_phased.vcf \
  	       -L $prefix/phased_files/${chr}/shapeit_conv_${chr}.log
    else
      # Run SHAPEIT to phase the chromosome (no need to indicate anything else, just the input file containing a single chromosome)
      if  [[ $family == No ]] ||  [[ $family == NO ]] ||  [[ $family == no ]] #If user indicated no relatedness
      then
           echo -e "Performing shapeit with NO relatedness"
           shapeit -B $prefix/prephasing_files/${chr}/${chr}_${base}_filtered \
  	      -M $ShapeRef/genetic_map_chr${chr}_combined_b37.txt \
  	      -O $prefix/phased_files/${chr}/${chr}_${base}_filtered_phased \
  	      --thread $cpus \
  	      -L $prefix/phased_files/${chr}/shapeit_${chr}.log
      # Convert output files from SHAPEIT (.haps) to vcf format
           shapeit -convert \
            --input-haps $prefix/phased_files/${chr}/${chr}_${base}_filtered_phased \
            --output-vcf $prefix/phased_files/${chr}/${chr}_${base}_filtered_phased.vcf \
  	       -L $prefix/phased_files/${chr}/shapeit_conv_${chr}.log
      elif [[ $family == Yes ]] || [[ $family == YES ]] || [[ $family == yes ]]
      then
           echo -e "Performing shapeit taking into account relatedness"
  	     shapeit -B $prefix/prephasing_files/${chr}/${chr}_${base}_filtered \
             -M $ShapeRef/genetic_map_chr${chr}_combined_b37.txt \
             --duohmm \
             -W 5 \
             --force \
             -O $prefix/phased_files/${chr}/${chr}_${base}_filtered_phased \
             --thread $cpus \
  	       -L $prefix/phased_files/${chr}/shapeit_${chr}.log
           # Convert output files from SHAPEIT (.haps) to vcf format
           shapeit -convert \
             --input-haps $prefix/phased_files/${chr}/${chr}_${base}_filtered_phased \
             --output-vcf $prefix/phased_files/${chr}/${chr}_${base}_filtered_phased.vcf \
  	       -L $prefix/phased_files/${chr}/shapeit_conv_${chr}.log
      else
          echo -e "Shapeit stopped due to unknown family parameter"
          exit		
      fi
    fi
}
