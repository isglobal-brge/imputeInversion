#!/bin/bash

#######################################
## Run minimac Imputation alone
#######################################

## Load scripts
path=`dirname "${BASH_SOURCE[0]}"`
source $path/minimac3_imputation.sh

## Load config file
source $path/config

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
elif [[ " ${RANGES[*]} " == *" ${RANGES[$inversion]} "* ]] # If the input indicating the inversion is the name of the inversion (name of the inversion is a key in the associative array RANGES)
then
  chr=${RANGES[$inversion]%:*} # Create variable 'chr' with the chromosome number (placed before the ':' in the associative array RANGES)
  start=${RANGES[$inversion]##*:} # Create variable 'start' with the inversion starting position (first, select the part after the ':' in the associative array RANGES)
  start=${start%-*} # Select the part before the '-' of the start variable (that had the start-end positions)
  end=${RANGES[$inversion]##*-} # Create variable 'end' with the inversion ending position (placed after the '-' in the associative array RANGES)
  prefix=$inversion # Keep the name of the inversion in the variable 'prefix' to use it in the final output files
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
  family="NO"
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

dir=`dirname $data`
base=`basename $data`

impute $chr $dir $base $prefix $start $end $minimacRefs $cpus
