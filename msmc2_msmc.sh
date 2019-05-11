#!/bin/bash
#PBS -N msmc2 
#PBS -l walltime=16:00:00
#PBS -l select=1:ncpus=3:mem=20gb:scratch_local=50gb
#PBS -j oe

###this script will take an s argument representing sample name

#s argument
echo "processing sample $s at date `date`"
##you have to submit the job from a parent directory created by script msmc0 in current case it is msmc_inv
mtpath=${PBS_O_WORKDIR}/${s}_dir/multihetsep # fullpath to multihetsep
#msmc="/storage/brno3-cerit/home/janav/programs/msmc-master/msmc"
msmc="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/msmc2/build/release/msmc2" # changed for msmc2
bs="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/msmc-tools-master/multihetsep_bootstrap.py" # bootstrap directory
nt=3
outdir=msmc2_15_segments
echo $mtpath
module add gcc-5.3.0
module add debian7-compat
module add debian8-compat
module add R-3.4.3-gcc
module add gsl-1.16-gcc

###copy multihetsep files to scratch

cp ${mtpath}/*hetsep.txt $SCRATCHDIR
cd $SCRATCHDIR
echo "hetsep copied at `date`"
ls

echo "delete following empty hetsep:"
find . -name '*hetsep.txt' -empty -type f
find . -name '*hetsep.txt' -empty -type f -delete
printf "\n\nbegin msmc\n"
mkdir $outdir
#####msmc_line#####
$msmc -t $nt -p 1*2+15*1+1*2 -o ${outdir}/${s} *hetsep.txt # --fixedRecombination 1*3+20*1+1*3  1*2+15*1+1*2 

####bootstrap####
# need to get number of chromosomes
nchr=`ls *hetsep.txt | wc -l`
$bs --nr_chromosomes $nchr ${outdir}/${s} *hetsep.txt # this command only generates bootstrapped hetseps
# now you have to run msmc for each hetsep question is if I should parallelize it somehow, there are 20 bootstraps
parallel -j 2 "$msmc -t $nt -p 1*2+15*1+1*2 -o ${outdir}/${s}.bs.{} ${outdir}/${s}_{}/*txt" ::: `seq 1 20` > ${outdir}/bootstrap_msmc.log

# concatenate the bootstraps
cd $outdir
cat *bs.*final.txt | grep -v time > ${s}.bs.final.concat.txt
###plot the Ne line#####
module rm gsl-1.16-gcc
module rm debian7-compat
module rm debian8-compat
Rscript /storage/brno3-cerit/home/vlkofly/scripts/msmc3_plot.R ${s}.final.txt ${s}.bs.final.concat.txt #modify to include the bootstrap
cd ../

cp -r $outdir $PBS_O_WORKDIR/${s}_dir/  || export CLEAN_SCRATCH=false
echo "msmc of $s done"
