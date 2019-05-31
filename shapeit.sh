#!/bin/bash
#PBS -N shapeit_dpmiss
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=1:mem=220gb:scratch_local=150gb
#PBS -j oe

# read backed phasing by shapeit
# example of running the script: qsub -v 'vcf=scaffold_1.vcf.gz,scf=1,outdir=blabla' ~/scripts/shapeit.sh

#scf #  main variable argument in the script the script will take plink data based on that and map; integer only

echo "started at `date`"

#inscript variables definition segment:
wd=$PBS_O_WORKDIR
nt=1 # when low number of samples do not multithread
sbin=/storage/brno3-cerit/home/filip_kolar/programs/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit
rbin=/storage/brno3-cerit/home/filip_kolar/programs/extractPIRs.v1.r68.x86_64/extractPIRs
recrate="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/mapped_to_melanotis/phasing/genetic_map/genetic.map.melanotis.tsv" #full path to a file with recombination rate map 
recmap="-M recmap.${scf}.tsv" # if recombination map available
#recmap="--rho 0.0001" # uniform recombination rate
sambam="sample_bam.txt" # relative path from wd to file defining sample names and corresponding bam files
#sample_name\t fullpath_to_bam
#a2th2_01	/storage/plzen1/home/vlkofly/rgadd2/a2th2_01.final.bam



module add parallel-20160622
module add htslib-1.9
module add vcftools-0.1.15
module add bcftools-1.9
trap 'clean_scratch' TERM EXIT

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi

cp ${wd}/$vcf $SCRATCHDIR
cp ${wd}/${vcf}.tbi $SCRATCHDIR
cp $recrate $SCRATCHDIR  #copy genetic map
cp ${wd}/${sambam} $SCRATCHDIR  #copy the bam files that will be used for read backed phasing

cd $SCRATCHDIR
recrate=`basename $recrate`
vcf=`basename $vcf`

cut -f 2 sample_bam.txt | parallel -j 1 "cp {} $SCRATCHDIR/; cp {.}.bai $SCRATCHDIR/" # copy all the bam files into SCRATCH 



echo "data imported succesfully at `date`:"
du -ha *.ba*

sample_num=`cat sample_bam.txt | wc -l`
echo "number of samples $sample_num"

#parse bam information file for a given chromosome, basename of path and add chrom identifier
cut -f 2 sample_bam.txt| sed 's!.*/!!' > path.tmp
yes ${scf} | head -n $sample_num > chr.tmp
cut -f 1 sample_bam.txt > id.tmp

paste id.tmp path.tmp chr.tmp > bam_define.txt

head bam_define.txt # check how the 
rm *tmp

# the problem with duplication of sites, try to solve it by splitting the vcf
tabix -h $vcf $scf | bgzip -c > ${scf}.vcf.gz
vcf=${scf}.vcf.gz

#### the missingness is also a problem, in case of GM it is sample I370 that has high missingness that is not allowed in shapeit
echo "missingness before filtering:"
vcftools --gzvcf $vcf --missing-indv --stdout | tee raw.mis.tmp # print missingness per individual
vcftools --gzvcf $vcf --max-missing 0.9 --recode --stdout | bgzip -c > ${scf}.0.9miss.vcf.gz # to solve the problem with missingness per site
bcftools view -Oz -e 'FORMAT/DP[4] < 6' -o ${scf}.lowmis.vcf.gz ${scf}.0.9miss.vcf.gz # remove sites that have depth in sample 4 (I370) lower than 6 to solve missingness per individual TODO generalize it.
vcf=${scf}.lowmis.vcf.gz
echo "missingness after filtering:"
vcftools --gzvcf $vcf --missing-indv --stdout | tee filt.mis.tmp # print missingness per individual befor

N_raw=`tail -n 1 raw.mis.tmp | cut -f 2`
N_filt=`tail -n 1 filt.mis.tmp | cut -f 2`

printf "\nfiltstats\t%s\t%d\t%d\n" $scf $N_raw $N_filt
printf "\nPIRextraction follows\n"
# extract PIRs
$rbin --bam bam_define.txt --vcf $vcf --out pirlist --base-quality 20 --read-quality 20 # increased from 20 to 25 for scaffold_3 where I was getting that sequencing underflow error.

 printf "\npirlist generated at `date`\nhead of pirlist follows:"

head -n 50 pirlist
np=`wc -l pirlist `

printf "\nNO_pirs:\t%s\n" $np


#parse genetic map for a given chromosome

grep -P "${scf}\t" $recrate | cut -f 2- > recmap.${scf}.tsv
head -n 1 $recrate | cut -f 2- > header.map.tsv
cat header.map.tsv recmap.${scf}.tsv | sponge recmap.${scf}.tsv

printf "\nPhasing itself follows:\n"

#rm header.map.tsv

# shapeit with readback # when running for arabidopsis:

#$sbin -assemble --input-vcf $vcf --input-pir pirlist  \
#$recmap \
#-O ${scf}.PIR.phased \
#-T $nt -L PIR_based.log \
#--window 0.5 \
#--effective-size 240000 \
#--states 350 \
#--burn 20 \
#--prune 20 \
#--main 80

# shapeit with readback when run for mockers
# use just one thread for mockers as there are low number of individuals
# the effective size of all the populations piled atop?
$sbin -assemble --input-vcf $vcf --input-pir pirlist  \
$recmap \
-O ${scf}.PIR.phased \
-T $nt -L ${scf}.PIR_based.log \
--window 1 \
--effective-size 10000 \
--states 350 \
--burn 30 \
--prune 30 \
--main 100


$sbin -convert --input-haps ${scf}.PIR.phased --output-vcf ${scf}.phased.vcf.gz # convert the phasing into vcf

rm *bam
rm *bai

mkdir $outdir 
mv * $outdir/

cp -r $outdir ${wd}/ || export CLEAN_SCRATCH=false

#cp $SCRATCHDIR/* $wd/filtervcf_toxo/ || export CLEAN_SCRATCH=false
echo "phasing of $scf done and results exported at `date`"

