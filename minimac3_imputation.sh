#!/bin/bash

#######################################
# Impute data with minimac
# Arguments:
#   $1: Chromosome to be phased
#   $2: folder path
#   $3: base name
#   $4: prefix to name inversions
#   $5: inversion start
#   $6: inversion end
#   $7: Path to reference folder
#   $8: number of CPUs used in minimac
# Returns:
#   None
#######################################
impute(){
  
  chr=$1
  dir=$2
  base=$3
  prefix=$4
  start=$5
  end=$6
  pathRef=$7
  cpus=$8

  if [ ! -d $dir/pimputed_files/$prefix ]; then
    mkdir $dir/pimputed_files/$prefix # Create file with the chromosome number to store the files separated by chromosomes
  fi
  
  if [[ $chr == X ]] # If it is the chromosome X
  then
    sed 's/^23/X/' $dir/phased_files/${chr}/${chr}_${base}_filtered_phased.vcf  > $dir/phased_files/${chr}/${chr}_${chr}_filtered_phased_X.vcf # Change number of the chromosome (it was indicated with 23, change to X)
    awk '$5 == 2 { print $2 }' < $dir/${base}.fam > $dir/females_${base}.txt # Store ids of the females individuals in a text file
    awk '$5 == 1 { print $2}'  < $dir/${base}.fam > $dir/males_${base}.txt   # Store ids of the males individuals in a text file
    # Separate the data by sex using the above files. Store males individuals in a separated file
    vcftools --vcf $dir/phased_files/${chr}/${chr}_${base}_filtered_phased_X.vcf \
            --keep $dir/males_${base}.txt  \
            --recode \
            --out $dir/phased_files/${chr}/${chr}_${base}_males_filtered_phased_X
    # Store females in another vcf file
    vcftools --vcf $dir/phased_files/${chr}/${chr}_${base}_filtered_phased_X.vcf \
            --keep $dir/females_${base}.txt \
            --recode \
            --out $dir/phased_files/${chr}/${i}_${base}_females_filtered_phased_X
    # Impute the inversion region of the chromosome X (indicating by --start and --end) for the males
    Minimac3-omp --refHaps $pathRef/invX_006.vcf.gz \
            --haps $dir/phased_files/${chr}/${chr}_${base}_males_filtered_phased_X.recode.vcf \
            --rsid \
            --format GT,GP \
            --chr X \
            --start $start \
            --end $end \
            --prefix $prefix/${prefix}_${base}_males_imputed \
            --cpus $cpus
    # Impute the inversion region of the chromosome X (indicated by --start and --end) for the females
    Minimac3-omp --refHaps $pathRef/invX_006.vcf.gz \
            --haps $dir/phased_files/${chr}/${chr}_${base}_females_filtered_phased_X.recode.vcf \
            --rsid \
            --format GT,GP \
            --chr X \
            --start $start \
            --end $end \
            --prefix ${prefix}_${base}_females_imputed \
            --cpus $cpus
          
  else # For the rest of the chromosomes
    # Impute the inversion region (indicated by --start and --end)
    Minimac3-omp --refHaps $pathRef/${prefix}.vcf.gz \
            --haps $dir/phased_files/${chr}/${chr}_${base}_filtered_phased.vcf \
            --rsid \
            --format GT,GP \
            --chr $chr \
            --start $start \
            --end $end \
            --prefix ${prefix}_${base}_imputed \
            --cpus $cpus
  fi    
 
  mv $base* $dir/pimputed_files/$prefix/
}
