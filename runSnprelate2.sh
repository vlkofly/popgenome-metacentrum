#! /bin/bash
#PBS -N snprelate_final
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=1:mem=100gb:scratch_local=100gb
#PBS -j oe

# script adjusted in july 2018, to take -v vcf=vcf_to_analyse as a parameter and then output results into a folder snprelate
# qsub -v 'vcf=mocker_to_toxo_snps_fin.vcf.passed.gz' ../../../scripts/runSnprelate.sh
# supply dp mask in variable dp give relative path from wd

# adjusted in august 15 from runSnprelate.sh > runSnprelate2.sh
# the idea is that I will do final filtering for neutral sites and mask Z chromosome
# then use the GDS file to built a phylogenetic tree by snphylo
# will take filtered vcf file from ROH



wd=$PBS_O_WORKDIR
scriptdir="/storage/brno3-cerit/home/vlkofly/scripts/snp_pca.R"
script=`basename $scriptdir`
batchdir="snprelate_chr9"
sn="/storage/brno3-cerit/home/janav/renaming.old.new.txt"
#scaf="/storage/brno3-cerit/home/janav/scaffolds_demography_27M.txt"
#Zchrom="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/reference/toxostoma/blast_sex/z_scaffolds_in_toxo.list.txt" # scaffolds to remove due to their affinity with sex Z chromosome
#dp_mask="/storage/brno3-cerit/home/janav/vcf_toxo/fin/strict_biallelic/vartables/depth.blacklist.mask.bed"
#dp_mask=${wd}/${dp}
coding="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/mapped_to_melanotis/snp_results/neutral.1.bed" #positive mask with putatively netural sites in bed format, do some bedjobs
vcfout="neutral_phasedchr9.vcf.gz"


module add vcftools-0.1.15
module add bcftools-1.6
module add debian7-compat
module add debian8-compat
module add R-3.4.3-gcc

trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
cp $wd/$vcf $SCRATCHDIR # copy vcf file from V3 to scratch
cp $coding $SCRATCHDIR # copy vcf file from V3 to scratch
cp $scriptdir $SCRATCHDIR
cp /storage/brno3-cerit/home/vlkofly/mockingbird_genomics/sample_lists/samplelist_parvulus_new.txt $SCRATCHDIR


cd $SCRATCHDIR
echo "files copied at `date` and filltering follows:"
echo "$vcf taken as input and filtered vcf produced: $vcfout"
ls
#Zchrom=`basename $Zchrom`
coding=`basename $coding`
vcf=`basename $vcf`

# run it for a particular population set
for f  in `cat samplelist_parvulus_new.txt`
do
	ind+=" --indv $f"
done

#vcftools --vcf $vcf --not-chr Mmel_ChrZ --max-missing 1 --minDP 7 --exclude-bed $coding --recode --stdout > $vcfout
#vcftools --gzvcf $vcf --not-chr Mmel_ChrZ  $ind --max-missing 1 --minDP 7 --exclude-bed $coding --recode --stdout > $vcfout

#bcftools view -Oz -T ^$Zchrom -e 'GT="./."' $vcf | bcftools view -Ov -T ^$coding -o $vcfout 

inputNo=`zcat $vcf grep -v "#" | wc -l` #some chekpoint of vcf filtering
finterNo=`grep -v "#" $vcfout | wc -l`


echo "input size: $inputNo filter size: $finterNo filterd at `date` rscript follows"

module rm debian7-compat
module rm debian8-compat

Rscript $script $vcf 
ls
rm $vcf
rm $script
mkdir $batchdir
mv * $batchdir/
cp -r $batchdir $wd/ || export CLEAN_SCRATCH=false
echo "snprelate analysis of $vcf done at `date`"
