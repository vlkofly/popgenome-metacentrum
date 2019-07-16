#!/bin/bash
#PBS -N lastz_align 
#PBS -l select=1:ncpus=9:mem=250gb:scratch_local=100gb
#PBS -l walltime=40:00:00
#PBS -j oe

# align two genomes/chromosomes/contigs in fasta: f1 - target f2 - query folder with splitted fasta files  (supplied as arguments in pbs manners -v 'f1=xx.fa,f2=xy.fa')


wd=$PBS_O_WORKDIR
export bin="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/bin/lastz"
outd="ficedula_sturnus_perchr"


trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi

cp -r $wd/$f1 $SCRATCHDIR/
cp -r $wd/$f2 $SCRATCHDIR/

cd $SCRATCHDIR
ext=${f1##*.}

#name=melanotis_ficedula.maf

#### split fasta files into separate chromosome files

# done beforehand


# now align each to each by parallel

laf () {
f1=$1
f2=$2
f1base=`basename $f1`
f2base=`basename $f2`
name=${f1base}_${f2base}.maf

$bin ${f1}[nameparse=full] ${f2}[nameparse=full]\
       	--notransition --gfextend --hspthresh=2200  --chain --gapped --ydrop=3400 --gappedthresh=6000 --inner=2000 \
	--rdotplot=${name/maf/rdotplot} --progress=1000 --format=maf --output=${name}

echo "alignment of $name finished at `date`"
}

export -f laf

f1list=`ls ${f1}/*.fa`
f2list=`ls ${f2}/*.fa`

echo $f1list

parallel -j 8 laf ::: $f1list ::: $f2list

find . -empty -delete

# here you have to add a find command that will remove small files


du -ha *maf | sort -hr > sorted.size.mafs.list


ls *rdotplot | parallel -j 8 Rscript ~/scripts/plot_rdotplot.R -i {} -o {}.png 
rm -r $f1
rm -r $f2
mkdir $outd
mv * $outd
cp -r $outd $wd/ || export CLEAN_SCRATCH=false

