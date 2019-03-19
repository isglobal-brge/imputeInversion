#!/bin/bash

## Exit if any line fails
set -e


## Load scripts
path=`dirname "${BASH_SOURCE[0]}"`
source $path/prephasing.sh

## Load config file
source $path/config

# Define usage(). It will show the options of the function and the valid inputs when the user introduces a wrong option, does not introduce a necessary input or uses the -h|--help for help
usage() { echo -e "\nUsage: $0 \n\n\t-d,--data <FileName>\n\t\tFile name without extension (format must be PLINK - bed/bim/fam)\n\n\t-i,--inversions <File> \n\t\tFile containing inversion ranges | Inversion range | Inversion name | 'All'\n\n\t--family <Yes/No>\n\t\tIndicator for family relatedness\n\n\t--threads <Integer>\n\t\tNumber of threads to be used while phasing with SHAPEIT and while imputing with Minimac3\n\n\t-k,--keep-files <Yes/No>\n\t\tIntermediate files will be removed by default. To keep them, write 'Yes', 'YES' or 'yes'\n\n\t-h,--help \n\t\tShows this message\n\nValid inversion names: ${SORTEDKEYS[@]}\n" 1>&2; exit 1; }

# Define an associative array between the inversion names and their ranges
declare -A RANGES=( ["inv8_001"]="8:8055789-11980649" ["inv17_007"]="17:43661775-44372665" ["inv7_005"]="7:54290974-54386821" ["invX_006"]="X:72215927-72306774" 
                    ["inv12_004"]="12:47290470-47309756" ["inv7_011"]="7:70426185-70438879" ["inv7_003"]="7:31586765-31592019" ["inv11_001"]="11:41162296-41167044" 
                    ["inv2_013"]="2:139004949-139009203" ["inv6_006"]="6:130848198-130852318" ["inv3_003"]="3:162545362-162547641" ["inv7_014"]="7:151010030-151012107" 
                    ["inv11_004"]="11:66018563-66019946" ["inv1_008"]="1:197756784-197757982" ["inv21_005"]="21:28020653-28021711" ["inv12_006"]="12:71532784-71533816" 
					["inv6_002"]="6:31009222-31010095" ["inv14_005"]="14:65842304-65843165" ["inv1_004"]="1:92131841-92132615" 
                    ["inv2_002"]="2:33764554-33765272" )

# Define the reverse associative array above, to use the names of the inversions in the ouput files                    
declare -A REVRANGES=( ["8:8055789-11980649"]="inv8_001" ["17:43661775-44372665"]="inv17_007" ["7:54290974-54386821"]="inv7_005" ["X:72215927-72306774"]="invX_006" 
                       ["12:47290470-47309756"]="inv12_004" ["7:70426185-70438879"]="inv7_011" ["7:31586765-31592019"]="inv7_003" ["11:41162296-41167044"]="inv11_001" 
                       ["2:139004949-139009203"]="inv2_013" ["6:130848198-130852318"]="inv6_006" ["3:162545362-162547641"]="inv3_003" ["7:151010030-151012107"]="inv7_014" 
                       ["11:66018563-66019946"]="inv11_004" ["1:197756784-197757982"]="inv1_008" ["21:28020653-28021711"]="inv21_005" ["12:71532784-71533816"]="inv12_006" 
					   ["6:31009222-31010095"]="inv6_002" ["14:65842304-65843165"]="inv14_005" ["1:92131841-92132615"]="inv1_004" ["2:33764554-33765272"]="inv2_002" )

# Define array with the keys of the RANGES associative array                    
KEYS=(${!RANGES[@]})
# Define an array containing the keys (inversion names) sorted, to show them in order when the usage message is printed (showing the valid inversion names)
IFS=$'\n' SORTEDKEYS=($(sort -t v -k 2 -g<<<"${KEYS[*]}"))
unset IFS

# Assign the arguments to variables
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -d|--data)
    data="$2"
    shift # past argument
    shift # past value
    ;;
    -i|--inversions)
    inversion="$2"
    shift # past argument
    shift # past value
    ;;
    -t|--threads)
    cpus="$2"
    shift # past argument
    shift # past value
    ;;
    --family)
    family="$2"
    shift # past argument
    shift # past value
    ;;
    -k|--keep-files)
    keep_files="$2"
    shift # past argument
    shift # past value
    ;;
    -h|--help|*)
    usage
    shift # past argument
    ;;
esac
done

## Exit if no data provided
if [[ -z "$data" ]] # If $data is empty (the name of the files was not provided)
then
  echo -e "\nNo input data to be imputed" # Show a message indicating the error
  usage # Exit and show the usage message
fi

## Create variables to keep the coordinates and the names of the inversion(s) to impute
if [[ -z "$inversion" ]] # If $inversion is empty (no input indicating inversion to impute)
then
  echo -e "\nNo input indicating inversion(s) region(s)" # Show a message indicating the error
  usage # Exit and show the usage message
elif [[ $inversion == *.txt ]] # If the input indicating the inversion(s) is a text file (should contain the coordinates of the inversions as indicated in the manual)
then 
  chr=( $(awk -F '[ \t:-]+' '{print $1}' $inversion) ) # Keep the first column (columns can be separated by spaces, tabs, : or -) of the file (chromosome) in the array 'chr' 
  start=( $(awk -F '[ \t:-]+' '{print $2}' $inversion) ) # Keep the second column of the file (inversion starting position) in the array 'start'
  end=( $(awk -F '[ \t:-]+' '{print $3}' $inversion) ) # Keep the third column of the file (inversion ending position) in the array 'end'
  counter=0 # Start counter to give names to the inversions
  getnames=() # Create array that stores the inversions to be imputed in the format chr:start-end (will be used later to get the names from the REVERSERANGES)
  for i in ${chr[@]}
  do
    getnames+=( "${chr[$counter]}:${start[$counter]}-${end[$counter]}" )
    counter=$counter+1
  done
elif [[ " ${RANGES[*]} " == *" ${RANGES[$inversion]} "* ]] # If the input indicating the inversion is the name of the inversion (name of the inversion is a key in the associative array RANGES)
then
  chr=${RANGES[$inversion]%:*} # Create variable 'chr' with the chromosome number (placed before the ':' in the associative array RANGES)
  start=${RANGES[$inversion]##*:} # Create variable 'start' with the inversion starting position (first, select the part after the ':' in the associative array RANGES)
  start=${start%-*} # Select the part before the '-' of the start variable (that had the start-end positions)
  end=${RANGES[$inversion]##*-} # Create variable 'end' with the inversion ending position (placed after the '-' in the associative array RANGES)
  prefix=$inversion # Keep the name of the inversion in the variable 'prefix' to use it in the final output files
elif [[ " ${RANGES[*]} " == *" $inversion "* ]] # If the input indicating the inversion is the coordinates of the inversion
then
  chr=${inversion%:*} # Keep the chromosome number in 'chr' (which is before the ':' in the coordinates)
  start=${inversion##*:} # Keep what is after the ':' in the coordinates
  start=${start%-*} # Keep the starting position in 'start' (which is before the '-' in the start variable above)
  end=${inversion##*-} # Keep the ending position in 'end' (which is after the '-' in the coordinates)
  for i in "${!RANGES[@]}" # For all the keys in the associative array RANGES...
  do
    if [[ " ${RANGES[$i]} " == " $inversion " ]] # ...check if that key is associated to the coordinates indicated in the input...
    then
      prefix=$i # ...and store that key (inversion name) in the variable 'prefix' to use it in the final output files
    fi
  done
elif [[ $inversion == all ]] || [[ $inversion == All ]] || [[ $inversion == ALL ]] # If the input indicating the inversions to impute is ALL|All|all
then
  chr=() # Create array 'chr'
  start=() # Create array 'start'
  end=() # Create array 'end'
  for i in ${RANGES[@]} # For every element in the associative array RANGES
  do
    chr+=( ${i%:*} ) # Keep the chromosome number in 'chr' (which is before the ':' in the coordinates)
    middle=${i##*:} # Keep what is after the ':' in the coordinates
    start+=( ${middle%-*} ) # Keep the starting position in 'start' (which is before the '-' in the 'middle' variable above)
    end+=( ${i##*-} ) # Keep the ending position in 'end' (which is after the '-' in the coordinates)
  done
  counter=0 # Start counter to give names to the inversions
  getnames=() # Create array that stores the inversions to be imputed in the format chr:start-end (will be used later to get the names from the REVERSERANGES)
  for i in ${chr[@]}
  do
    getnames+=( "${chr[$counter]}:${start[$counter]}-${end[$counter]}" )
    counter=$counter+1
  done
else # If the input does not correspond with any of the previous options, it is not valid 
  echo -e "\nInput indicating inversion(s) region(s) is NOT valid" # Show a message indicating the error
  usage # Exit and show the usage message
fi

# Store names if there are multiple inversions to impute
if [[ -z "$prefix" ]] # If $prefix is empty (only if there is more than one inversion to impute, following the code above)
then
  prefix=()
  for i in "${getnames[@]}" # For each inversion to impute...
  do
    if [[ " ${REVRANGES[*]} " == *" ${REVRANGES[$i]} "* ]]  # ...if it is a key in the REVRANGES associative array...
    then
      prefix+=( ${REVRANGES[$i]} ) # ...store its value of the associative array REVRANGES
    fi
  done
fi

unique_chr=( $(echo "${chr[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') ) # Create an array with uniques chromosome numbers, this will be used for prephasing and phasing (there are more than one inversion in some chromosomes) 

## Relatedness between samples
if [[ -z "$family" ]] # If no relatedness was set ($family empty)
then
  echo -e "\n Samples are not family related\n\n"
else
  echo -e "\n Some of the samples are family related\n\n"
fi

## Number of threads (By default: $cpus = 1)
if [[ -z "$cpus" ]] # If the number of threads to use is not set by the user ($cpus empty)
then
  echo -e "\nNumber of threads was not indicated. Value set to 1" # Show a message indicating that this value is set to one
  cpus=1 # Set to one the number of threads to be used by the program
elif ! [[ "$cpus" =~ ^[0-9]+$ ]] # If the input to indicate number of threads is not a number
then
  echo -e "\nInput for number of threads not valid. Value set to 1" # Show a message indicating that this value is set to one
  cpus=1 # Set to one the number of threads to be used by the program
fi
echo -e "Prephasing starting\n\n"
prephase $data $unique_chr $HRCref
. ~/data/software/imputeInversion-master/prephasing_v3.sh # Call the prephasing script to proceed with the process keeping the variables
. ~/data/software/imputeInversion-master/phasing_shapeit_v2.sh # Call the phasing script (keeping the variables)
. ~/data/software/imputeInversion-master/minimac3_imputation_v3.sh # Call imputation script (keeping the variables)
. ~/data/software/imputeInversion-master/postimputation_v2.sh # Call imputation script (keeping the variables)

# By default, the intermediate files will be deleted
if [[ -z "$keep_files" ]] || [[ $keep_files != Yes ]] && [[ $keep_files != YES ]] && [[ $keep_files != yes ]] # If user did not indicate 'YES' to keep the intermediate files 
then
  rm -rf prephasing_files # Remove folder containing pre-phasing files (remove files at this point to avoid having to many files at the end of the process if we want to impute all the chromosomes and cause potential memory problems) 
  rm -rf phased_files # Remove folder containing phased files 
  rm -rf pimputed_files # Remove folder containing imputed files
  rm -rf postimputation_files # Remove folder containing postimputed files (intermediate files between the imputed ones and the final ones with the correct ids)
fi