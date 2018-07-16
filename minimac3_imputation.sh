#!/bin/bash


mkdir pimputed_files # Folder to store all the files generated during phasing process

counter=0 # Start counter (will be used for the files names, starting and ending positions for each imputation)
for i in ${chr[@]} # For each element in 'chr' (there might be repited chromomes: more than one invresion per chromosome to be imputed)
do
  mkdir pimputed_files/${prefix[$counter]} # Create file with the chromosome number to store the files separated by chromosomes
  if [[ $i == X ]] # If it is the chromosome X
  then
    sed 's/^23/X/' phased_files/${i}/${i}_${data}_filtered_phased.vcf  > phased_files/${i}/${i}_${data}_filtered_phased_X.vcf # Change number of the chromosome (it was indicated with 23, change to X)
    awk '$5 == 2 { print $2 > "females.txt"}' ${data}.fam # Store ids of the females individuals in a text file
    awk '$5 == 1 { print $2 > "males.txt"}' ${data}.fam # Store ids of the males individuals in a text file
    # Separate the data by sex using the above files. Store males individuals in a separated file
    vcftools --vcf phased_files/${i}/${i}_${data}_filtered_phased_X.vcf \
            --keep males.txt\
            --recode \
            --out phased_files/${i}/${i}_${data}_males_filtered_phased_X
    # Store females in another vcf file
    vcftools --vcf phased_files/${i}/${i}_${data}_filtered_phased_X.vcf \
            --keep females.txt\
            --recode \
            --out phased_files/${i}/${i}_${data}_females_filtered_phased_X
    # Impute the inversion region of the chromosome X (indicating by --start and --end) for the males
    Minimac3-omp --refHaps ~/data/PublicData/REFERENCES/reference_panels/ALL.chrX.Non.Pseudo.Auto.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz \
            --haps phased_files/${i}/${i}_${data}_males_filtered_phased_X.recode.vcf \
            --rsid \
            --format GT,DS,GP \
            --chr X \
            --start ${start[$counter]} \
            --end ${end[$counter]} \
            --prefix pimputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_imputed \
            --cpus $cpus
    # Impute the inversion region of the chromosome X (indicated by --start and --end) for the females
    Minimac3-omp --refHaps ~/data/PublicData/REFERENCES/reference_panels/ALL.chrX.Non.Pseudo.Auto.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz \
            --haps phased_files/${i}/${i}_${data}_females_filtered_phased_X.recode.vcf \
            --rsid \
            --format GT,DS,GP \
            --chr X \
            --start ${start[$counter]} \
            --end ${end[$counter]} \
            --prefix pimputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_imputed \
            --cpus $cpus
          
  else # For the rest of the chromosomes
    # Impute the inversion region (indicated by --start and --end)
    Minimac3-omp --refHaps ~/data/PublicData/REFERENCES/reference_panels/ALL.chr${i}.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz \
            --haps phased_files/${i}/${i}_${data}_filtered_phased.vcf \
            --rsid \
            --format GT,DS,GP \
            --chr $i \
            --start ${start[$counter]} \
            --end ${end[$counter]} \
            --prefix pimputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed \
            --cpus $cpus
  fi    
  counter=$counter+1
done  

# By default, the intermediate files will be deleted
if [[ -z "$keep_files" ]] || [[ $keep_files != Yes ]] && [[ $keep_files != YES ]] && [[ $keep_files != yes ]] # If user did not indicate 'YES' to keep the intermediate files
then
  rm -rf phased_files # Remove folder containing phased files 
fi

. postimputation.sh # Call imputation script (keeping the variables)