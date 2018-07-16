#!/bin/bash


mkdir postimputation_files # Folder to store all the files generated during the postimputation process
mkdir ${data}_imputed_files # Folder to store the final files (this one will NOT be removed)

counter=0 # Start counter (will be used for the files names, starting and ending positions for each imputation)
for i in ${chr[@]} # For each element in 'chr' (there might be repited chromomes: more than one invresion per chromosome to be imputed)
do
  mkdir postimputation_files/${prefix[$counter]} # Create file with the chromosome number to store the files separated by chromosomes
  mkdir ${data}_imputed_files/${prefix[$counter]}
  if [[ $i == X ]] # If it is the chromosome X
  then
    # Compress the male imputed file with bgzip to avoid problems using bcftools afterwards
    zcat pimputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_imputed.dose.vcf.gz | bgzip -c > postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_imputed_bgzip.vcf.gz && tabix postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_imputed_bgzip.vcf.gz
    # Store header in a text file to modify it
    bcftools view -h postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_imputed_bgzip.vcf.gz > postimputation_files/hdrmales.txt
    # Add lines to the header conatining information about ID=GENOTYPED (avoid problems in the next steps)
    sed -i '/^##FORMAT=<ID=GT/i ##FILTER=<ID=GENOTYPED,Description="Marker was genotyped AND imputed">\n##FILTER=<ID=GENOTYPED_ONLY,Description="Marker was genotyped but NOT imputed">' postimputation_files/hdrmales.txt
    # Create a temporal file with the new header and the data from the imputed file
    bcftools reheader -h postimputation_files/hdrmales.txt --output postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_hdr_imputed_bgzip.vcf.gz postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_imputed_bgzip.vcf.gz
    tabix -p vcf postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_hdr_imputed_bgzip.vcf.gz
    # Compress the female imputed file with bgzip to avoid problems using bcftools afterwards
    zcat pimputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_imputed.dose.vcf.gz | bgzip -c > postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_imputed_bgzip.vcf.gz && tabix postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_imputed_bgzip.vcf.gz
    # Store header in a text file to modify it
    bcftools view -h postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_imputed_bgzip.vcf.gz > postimputation_files/hdrfemales.txt
    # Add lines to the header conatining information about ID=GENOTYPED (avoid problems in the next steps)
    sed -i '/^##FORMAT=<ID=GT/i ##FILTER=<ID=GENOTYPED,Description="Marker was genotyped AND imputed">\n##FILTER=<ID=GENOTYPED_ONLY,Description="Marker was genotyped but NOT imputed">' postimputation_files/hdrfemales.txt
    # Create a temporal file with the new header and the data from the imputed file
    bcftools reheader -h postimputation_files/hdrfemales.txt --output postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_hdr_imputed_bgzip.vcf.gz postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_imputed_bgzip.vcf.gz
    tabix -p vcf postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_hdr_imputed_bgzip.vcf.gz
    # Merge files containing the female and the male individuals 
    bcftools merge --force-samples --output-type z --output postimputation_files/TEMP1.vcf.gz postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_hdr_imputed_bgzip.vcf.gz postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_hdr_imputed_bgzip.vcf.gz
    tabix -p vcf postimputation_files/TEMP1.vcf.gz
    
  else # For the rest of the chromosomes
    # Compress the imputed files with bgzip to avoid problems using bcftools afterwards
    zcat pimputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed.dose.vcf.gz | bgzip -c > postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed_bgzip.vcf.gz && tabix postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed_bgzip.vcf.gz    
    # Store header in a text file to modify it
    bcftools view -h postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed_bgzip.vcf.gz > postimputation_files/hdr.txt
    # Add lines to the header conatining information about ID=GENOTYPED (avoid problems in the next steps) 
    sed -i '/^##FORMAT=<ID=GT/i ##FILTER=<ID=GENOTYPED,Description="Marker was genotyped AND imputed">\n##FILTER=<ID=GENOTYPED_ONLY,Description="Marker was genotyped but NOT imputed">' postimputation_files/hdr.txt
    # Create a temporal file with the new header and the data from the imputed file
    bcftools reheader -h postimputation_files/hdr.txt --output postimputation_files/TEMP1.vcf.gz postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed_bgzip.vcf.gz
    tabix -p vcf postimputation_files/TEMP1.vcf.gz
  fi
  # Change ids by chr:position:ref:alt prior to set the rsnumber ids
  bcftools annotate --remove ID --set-id +'%CHROM:%POS:%REF:%ALT' --threads $cpus --output-type z --output postimputation_files/TEMP2.vcf.gz postimputation_files/TEMP1.vcf.gz
  tabix -p vcf postimputation_files/TEMP2.vcf.gz
  # Set rsnumber ids in the final file, using as a reference a file containing information about all the SNPs (in hg19, chr:pos:ref:alt and the rsnumber id)
  bcftools annotate --annotations ~/data/PublicData/REFERENCES/reference_panels/All_20180418.vcf.gz --columns ID --threads $cpus --output-type z --output ${data}_imputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed_final.vcf.gz postimputation_files/TEMP2.vcf.gz
  tabix -p vcf ${data}_imputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed_final.vcf.gz
  
  counter=$counter+1
done

rm postimputation_files/TEMP* # Remove temporal files created during the process
rm postimputation_files/hdr* # Remove files containing the headers

# By default, the intermediate files will be deleted
if [[ -z "$keep_files" ]] || [[ $keep_files != Yes ]] && [[ $keep_files != YES ]] && [[ $keep_files != yes ]] # If user did not indicate 'YES' to keep the intermediate files
then
  rm -rf pimputed_files # Remove folder containing imputed files
  rm -rf postimputation_files # Remove folder containing postimputed files (intermediate files between the imputed ones and the final ones with the correct ids)
fi

