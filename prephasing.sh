#!/bin/bash

#######################################
# Prephase data
# Arguments:
#   $1: path to PLINK data
#	$2: names of chromosomes to process
#	$3: References file for SNP testing
# Returns:
#   None
#######################################
function prephase(){
	
	local data=$1
	local unique_chr=$2
	local ref=`echo $3`
	
	local dir=`dirname $data`
	local base=`basename $data`
	local path=`dirname "${BASH_SOURCE[0]}"`

	
	## Create files for imputation (Imputation preparation)
	if [ ! -d $dir/prephasing_files ]; then
		mkdir $dir/prephasing_files # Folder to store all the files generated during pre-phasing process
	fi


	plink --bfile $data --freq --out $data # Generate freq file
	perl $path/HRC-1000G-check-bim.pl -b $data.bim -f $data.frq -r $ref -g 
	# Run script following the steps described in: http://www.well.ox.ac.uk/~wrayner/tools/
	# Following lines, part of the Run-plink.sh script that the above script generates (what this does is well described in the section 'HRC or 1000G Imputation preparation and checking' in http://www.well.ox.ac.uk/~wrayner/tools/)
	plink --bfile $data --exclude $dir/Exclude-$base-1000G.txt --make-bed --set-hh-missing --out $dir/TEMP1
	plink --bfile $dir/TEMP1 --update-map $dir/Chromosome-$base-1000G.txt --update-chr --make-bed --out $dir/TEMP2
	plink --bfile $dir/TEMP2 --update-map $dir/Position-$base-1000G.txt --make-bed --out $dir/TEMP3
	plink --bfile $dir/TEMP3 --flip $dir/Strand-Flip-$base-1000G.txt --make-bed --out $dir/TEMP4
	plink --bfile $dir/TEMP4 --reference-allele $dir/Force-Allele1-$base-1000G.txt --make-bed --out $dir/prephasing_files/${base}_updated
	rm $dir/TEMP* # Remove temporal files
	rm $dir/Run-plink.sh # Remove the script generated by the process above

	for i in ${unique_chr[@]} # For each single chromosome 
	do
		if [ ! -d $dir/prephasing_files/${i} ]; then
		  	  mkdir $dir/prephasing_files/${i}  # Create file with the chromosome number to store the files separated by chromosomes
	  fi
	  if [[ $i == X ]] # If it is the chromosome X
	  then
		# Select chromosome '23' and fitler it: --geno filters out all variants with missing call rates exceeding the provided value (0.05 to avoid problems with SHAPEIT) to be removed, while --mind does the same for samples
		plink --bfile $dir/prephasing_files/${base}_updated \
			--reference-allele $dir/Force-Allele1-$base-1000G.txt \
			--make-bed \
			--mind 0.05 \
			--geno 0.05 \
			--chr 23 \
			--out $dir/prephasing_files/${i}/${i}_${base}_filtered
	  else
		# Separate chromosome $i and fitler it: --geno filters out all variants with missing call rates exceeding the provided value (0.05 to avoid problems with SHAPEIT) to be removed, while --mind does the same for samples
		plink --bfile $dir/prephasing_files/${base}_updated \
			--reference-allele $dir/Force-Allele1-$base-1000G.txt \
			--make-bed \
			--mind 0.05 \
			--geno 0.05 \
			--chr $i \
			--out $dir/prephasing_files/${i}/${i}_${base}_filtered
	  fi
	done  


  # Move all files created during the pre-phasing process to the folder (to keep all the files in the same folder in case the user wants to keep the intermediate files)
  mv $dir/*-1000G.txt $dir/prephasing_files/
  mv $data.frq $dir/prephasing_files/
  mv $data.hh $dir/prephasing_files/
  mv $data.log $dir/prephasing_files/
}
