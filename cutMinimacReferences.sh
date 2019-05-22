#'#################################################################################
#'#################################################################################
#' Make minimac references specific per inversion and references for annotation
#'#################################################################################
#'#################################################################################

## Run in ~/PublicData/REFERENCES/reference_panels
mkdir imputeInversionRefsMinimac

### Index Refs
for i in `ls ALL*.vcf.gz`
do 
  tabix -p vcf $i
done

## Create file with ranges (R) --> imputeInversion folder (~/data/software/imputeInversion)
library(scoreInvHap)
library(GenomicRanges)

seqlevelsStyle(inversionGR) <- "NCBI"
df <- data.frame(as.character(inversionGR  + 5e5), stringsAsFactors = FALSE)
df$chr <- as.character(seqnames(inversionGR))

write.table(df, row.names = TRUE, col.names = FALSE, 
            quote = FALSE, "scoreInvHapRanges.txt")

invRanges=~/data/software/imputeInversion/scoreInvHapRanges.txt

while read -r LINE 
do
  inv=`echo $LINE | cut -d' ' -f1`
  range=`echo $LINE | cut -d' ' -f2`
  chr=`echo $LINE | cut -d' ' -f3`

  bcftools view -r $range ALL.chr${chr}.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz -o imputeInversionRefsMinimac/$inv.vcf.gz -O z 
done < $invRanges

## Inversion X
bcftools view -r X:71715927-72806774 ALL.chrX.Non.Pseudo.Auto.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz -o imputeInversionRefsMinimac/invX_006.vcf.gz -O z 

## Inversion 10
bcftools view -r 10:26720925-28156433 ALL.chr10.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz -o imputeInversionRefsMinimac/inv10_005.vcf.gz -O z 


## Annotation References
bcftools view -r `cut -d' ' -f2 $invRanges | paste -d, -s` All_20180418.vcf.gz -o imputeInversionAnnotRef.vcf.gz -O z 
