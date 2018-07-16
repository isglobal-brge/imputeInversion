#!/bin/bash


mkdir phased_files # Folder to store all the files generated during phasing process

for i in ${unique_chr[@]} # For each single chromosome 
do
  mkdir phased_files/${i} # Create file with the chromosome number to store the files separated by chromosomes
  if [[ $i == X ]] # If it is the chromosome X
  then
    # Run SHAPEIT to phase chromosome X (specifying the '--chrX' flag). SHAPEIT will phase female individuals and will set NA to haploid heterozygous males
    shapeit -B prephasing_files/${i}/${i}_${data}_filtered \
          -M ~/data/PublicData/REFERENCES/reference_panels/genetic_map_chrX_nonPAR_combined_b37.txt \
          -O phased_files/${i}/${i}_${data}_filtered_phased \
          --chrX \
          --thread $cpus
    # Convert output files from SHAPEIT (.haps) to vcf format      
    shapeit -convert \
          --input-haps phased_files/${i}/${i}_${data}_filtered_phased \
          --output-vcf phased_files/${i}/${i}_${data}_filtered_phased.vcf
  else
    # Run SHAPEIT to phase the chromosome (no need to indicate anything else, just the input file containing a single chromosome)
    shapeit -B prephasing_files/${i}/${i}_${data}_filtered \
          -M ~/data/PublicData/REFERENCES/reference_panels/genetic_map_chr${i}_combined_b37.txt \
          -O phased_files/${i}/${i}_${data}_filtered_phased \
          --thread $cpus
    # Convert output files from SHAPEIT (.haps) to vcf format    
    shapeit -convert \
          --input-haps phased_files/${i}/${i}_${data}_filtered_phased \
          --output-vcf phased_files/${i}/${i}_${data}_filtered_phased.vcf
  fi
done 

mv shapeit_* phased_files/ # Move shapeit .log files to the phased_files folder (to keep all the files in the same folder in case the user wants to keep the intermediate files)

# By default, the intermediate files will be deleted
if [[ -z "$keep_files" ]] || [[ $keep_files != Yes ]] && [[ $keep_files != YES ]] && [[ $keep_files != yes ]] # If user did not indicate 'YES' to keep the intermediate files 
then
  rm -rf prephasing_files # Remove folder containing pre-phasing files (remove files at this point to avoid having to many files at the end of the process if we want to impute all the chromosomes and cause potential memory problems) 
fi

. minimac3_imputation.sh # Call imputation script (keeping the variables)



