#!/bin/bash
#PBS -N salescan
#PBS -l walltime=15:00:00
#PBS -l select=1:ncpus=10:mem=120gb:scratch_local=50gb
#PBS -j oe

# follow up script after shapeit supply single argument vcf  

echo "started at `date`"
wd=$PBS_O_WORKDIR
outdir=selscan_result
nt=9
#vcf="mocker_to_toxo_1.merged.vcf.gz" # supply vcf as an argument
bin=/storage/brno3-cerit/home/filip_kolar/programs/selscan/bin/linux/selscan
norm=/storage/brno3-cerit/home/filip_kolar/programs/selscan/bin/linux/norm
maxnsl="300 900 1200"


module add parallel-20160622
trap 'clean_scratch' TERM EXIT

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi

#cp ${wd}/${plinkdir}/${ident}.bed $SCRATCHDIR
#cp ${wd}/${plinkdir}/${ident}.bim $SCRATCHDIR
#cp ${wd}/${plinkdir}/${ident}.fam $SCRATCHDIR
cp ${wd}/$vcf $SCRATCHDIR
#cp ${wd}/$map $SCRATCHDIR
cp /storage/brno3-cerit/home/filip_kolar/scripts/plot_selscan.R $SCRATCHDIR

cd $SCRATCHDIR


echo "data imported succesfully at `date`:"
ls

vcf=`basename $vcf`

parallel "$bin --nsl --max-extend-nsl {} --keep-low-freq --trunc-ok --vcf $vcf --out ${vcf}.selscan.window{}" ::: $maxnsl  # crucial is the max extend nsl - how far should I look for homozygosity

echo "selscan done at `date`"
ls

#scaffold_1.PIR.phased.BAB.final.vcf.selscan.window100.nsl.out
#parallel "$norm --nsl --files ${vcf}.selscan.window{}.nsl.out --bins 50 --log ${vcf}.50norm.wind.{}.log " ::: 100 # this changes nothing
parallel "$norm --nsl --files ${vcf}.selscan.window{}.nsl.out --bins 100 --log ${vcf}.50norm.wind.{}.log " ::: $maxnsl

#parse vcf format for shapeit that there are snp ids and chrom is integer, maybe not neccessary, but have to keep the consistency with map.


#parallel "Rscript plot_shapeit.R ${vcf}.selscan.window{}.nsl.out.50bins.norm" ::: 100
parallel "Rscript plot_selscan.R ${vcf}.selscan.window{}.nsl.out.100bins.norm" ::: $maxnsl


rm $vcf
#rm $map

mkdir $outdir 
mv * $outdir/

cp -r $outdir ${wd}/ || export CLEAN_SCRATCH=false

#cp $SCRATCHDIR/* $wd/filtervcf_toxo/ || export CLEAN_SCRATCH=false
echo "prep of $scf done and results exported at `date`"

