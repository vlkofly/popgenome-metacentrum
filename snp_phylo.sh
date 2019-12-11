#! /bin/bash
#PBS -N snp_phylo
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=1:mem=100gb:scratch_local=100gb
#PBS -j oe

### Program description:
# input GDS file as argument gds
# shpPhylo program infers phylogeny given the data supplied
# http://chibba.pgml.uga.edu/snphylo/


wd=$PBS_O_WORKDIR
bin="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/SNPhylo/snphylo.sh"
batchdir="snp_phylo"

#module add vcftools-0.1.15
#module add bcftools-1.6
#module add debian7-compat
#module add debian8-compat
module add R-3.4.3-gcc
module add muscle-3.8.31

trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
cp $wd/$gds $SCRATCHDIR # copy vcf file from V3 to scratch

cd $SCRATCHDIR
echo "files copied at `date` and filltering follows:"
echo "$gds taken as input and results send to $batchdir"

gds=`basename $gds`

$bin -d $gds -r -l 0.2 -b 1000 -o F_280 -P final_1000bootstrp  # later change for higher bootstrap

mkdir $batchdir
mv * $batchdir/
cp -r $batchdir $wd/ || export CLEAN_SCRATCH=false
echo "phylotree analysis of $gds done at `date`"
