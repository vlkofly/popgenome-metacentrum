#!/bin/bash
#PBS -N msmc0_split_scaffolds
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=16:mem=100gb:scratch_local=80gb
#PBS -j oe

####written in July 2018
# the script takes vcf with nonvariant and snp sites and scaff as arguments and outputs multihetsep files that are taken by msmc2_msmc.sh script as input
# all following variables need to be adjusted according individual user
# qsub -v 'vcf=merged_filtered/snp_complete_final.dphetmask.maxdp250.mindp8.merged.filtered.vcf.gz,scaff=../../../msmc_workdir/toxo_large_scaffoldsaa.intervals' ~/scripts/msmc0_split_scaff_depth_filt.sh

### inscript variables:
wd=$PBS_O_WORKDIR
vcf_base=`basename $vcf` # vcf file after filtering
scaff_base=`basename $scaff`
outdir=$wd/msmc # where to send the output full path
msmc_all_parser="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/msmc-tools-master/vcfAllSiteParser.py" # download from https://github.com/stschiff/msmc-tools
msmc_multihetsep="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/msmc-tools-master/generate_multihetsep.py"
#sample_names="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/sample_lists/sample.names.new.pol.txt" # have to modify this as well, get them automatically from vcf
batch_dir=msmc_complet
#module add gatk-3.7
module add bcftools-1.6
module add vcftools-0.1.15
module add htslib-1.6
module add parallel-20160622
module add debian7-compat
module add debian8-compat

#check if input arguments were supplied correctly
if [ -z "$vcf" ]; then
        echo "Error! Specify vcf file!"
        echo "qsub -v 'vcf=blablabla.vcf.gz'"
        exit 1
        fi

if [ -z "$scaff" ]; then
        echo "Error! Specify scaff file!"
        echo "qsub -v 'scaff=list.of.scaffolds.txt'"
        exit 1
        fi


trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 


cp $wd/$vcf* $SCRATCHDIR # copy vcf file  to scratch
cp $wd/$scaff $SCRATCHDIR # here you will supply the interval list of intervals larger than 10kb
#cp $sample_names $SCRATCHDIR # here you will supply the interval list of intervals larger than 10kb
cd $SCRATCHDIR

# apart from separating by chrom I will need to separate by sample as well probably by bcftools
sample_names=samples.txt
zcat $vcf_base| head -n 200 | grep CHROM | cut -f 10- | tr \\t \\n > $sample_names

mkdir $batch_dir # create a folder where separated chroms will go

echo "working directory created and filtering started at `date`"

ls

cd $batch_dir # the sample folders will be in this folder
echo "we are in `pwd`"
#this crazy loop/parallel/pipe separates a vcf file by individual and scaffold, filters out min and max depth and exports map and biallelic snps 
while read s
do
echo "processing sample $s"
mkdir ${s}_dir
ls
parallel -a ../$scaff_base -j 8 "bcftools view -Oz -s $s --trim-alt-alleles -r {} ../$vcf_base | bcftools view -Ov -e 'GT=\"./.\"' - | $msmc_all_parser {} ${s}_dir/${s}_{}.mask.bed.gz | gzip -c > ${s}_dir/${s}_{}.vcf.gz" # min and max depth values removed

done < ../$sample_names

echo "filtering done at `date`\n\n"

###so now do generate_multigetsep.py part, I will put it into this script with an option to run different parts of the scripts


while read s # for every sample generate hetsep files
do
	echo "processing sample $s \n\n"
	ls
	mkdir ${s}_dir/multihetsep
	parallel -a ../$scaff_base -j 8 " $msmc_multihetsep ${s}_dir/${s}_{}.vcf.gz --mask ${s}_dir/${s}_{}.mask.bed.gz > ${s}_dir/multihetsep/${s}_{}.hetsep.txt"
done < ../$sample_names

echo "multihetsep done at `date`"
cd ../ # go out of batchdir
rm $scaff_base
rm $sample_names
rm $vcf_base*
cp -r $batch_dir $wd/  || export CLEAN_SCRATCH=false
echo "first part of msmc analysis done for $vcf done"
