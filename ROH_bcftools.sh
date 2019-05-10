#!/bin/bash
#PBS -N ROH
#PBS -l walltime=16:00:00
#PBS -l select=1:ncpus=2:mem=50gb:scratch_local=80gb
#PBS -j oe

#supply vcf in -v vcf=jklfj argument the script will calculte various vcf statistics and outputs it into a folder with the vcfstats_* name.
# also supply dp as a depth mask
# qsub -v 'vcf=snp_complete_final.vcf.passed.gz,dp=vartables/complete_negative_depth.mask.bed' ../../../../../scripts/ROH_bcftools.sh 
nt='1'
#vcf="missing0.85_depth8to150.recode.vcf" # vcf file after filtering
#ref="Toxostoma_redivivum.genome.fa"
wd=$PBS_O_WORKDIR
sn="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/sample_lists/renaming.old.new.melan.txt" # reheadering does not work if used with all names
batch=melan
dp_path=$wd/$dp
dp_base=`basename $dp`
vcf_base=`basename $vcf`
#module add gatk-3.7
module add bcftools-1.6
module add vcftools-0.1.15
module add htslib-1.6
module add debian7-compat
module add debian8-compat

trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
cp $wd/$vcf $SCRATCHDIR # copy vcf file from V3 to scratch
cp $dp_path $SCRATCHDIR
#cp /storage/brno3-cerit/home/janav/scaffolds_demography_27M.txt $SCRATCHDIR
cd $SCRATCHDIR
echo "files copied to scratch at `date`"
ls
base=`echo $vcf_base | cut -d"." -f 1`


mkdir $batch # in the end you will move everything into this folder

###rename vcf
bcftools reheader -s $sn -o ${base}.renamed.vcf.gz $vcf_base #| bgzip -c > ${base}.renamed.vcf.gz
vcf=${base}.renamed.vcf.gz


vcftools --gzvcf $vcf_base --not-chr Mmel_ChrZ  --exclude-bed $dp_base --minDP 7 --recode --max-missing 0.75 --recode-INFO-all --stdout | bgzip -c > ${base}.masked.mindp7.0.75missing.noZ.vcf.gz # export in plink
vcf=${base}.masked.mindp7.0.75missing.noZ.vcf.gz
tabix -p vcf $vcf

#vcftools --vcf $vcf $chrlist --LROH --out allchr1Msubset_depth10to30

#bcftools needs frequency of nonref allele
bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' $vcf | bgzip -c > freqs.tab.gz
tabix -s1 -b2 -e2 freqs.tab.gz

bcftools roh  --AF-file freqs.tab.gz -o ${vcf}.roh  $vcf
#process it in R, generate graphs
rm $vcf_base
cp freqs.tab.gz ${batch}/
mv * ${batch}/
cp -r ${batch} $wd  || export CLEAN_SCRATCH=false
echo "filtering of $vcf done"
