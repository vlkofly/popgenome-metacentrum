#!/bin/bash
#PBS -N multiz_align 
#PBS -l select=1:ncpus=2:mem=10gb:scratch_local=100gb
#PBS -l walltime=10:00:00
#PBS -j oe

# align two pairwise alignments f1 f2  (supplied as arguments in pbs manners -v 'f1=xx.maf,f2=xy.maf')
# before this step you have to choose good alignments and concatenate all alignments of particular chrom of reference

# another input argument vcf 
# think about which vcf to take for polarisation the one without excdepth and repeats and without Z with some filter for missingness?
# no I will leave the vcf for some next script

#module add parallel

wd=$PBS_O_WORKDIR
outd="multiz_results_chr1"
export TMPDIR=$SCRATCHDIR # this did not solve the problem

trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi
cp -r $wd/$f1 $SCRATCHDIR/
cp -r $wd/$f2 $SCRATCHDIR/

cd $SCRATCHDIR
ext=${f1##*.}

f1=`basename $f1`
f2=`basename $f2`

# alter the names so multiz understands see maffilter
sed -e 's/^s >\([ficedulamnotsr]\+\):\(.\+\):1:+:[0-9]*\s/s \1.\2 /g' $f1 | sponge $f1
sed -e 's/^s >\([ficedulamnotsr]\+\):\(.\+\):1:+:[0-9]*\s/s \1.\2 /g' $f2 | sponge $f2
#sed -e 's/^s >\([ficedulamnotsr]\+\):\(.\+\):1:+:[0-9]*\s/s \1:\2 /g' $f2 | sponge $f2
#sed -e 's/^s >\([ficedulamnotsr]\+\):\(.\+\):1:+:[0-9]*\s/s \1:\2 /g' $f2 | sponge $f2

###possibly when there are multiple sources of target I should unify the name?

#sed -e 's/^s sturnus.[A-Za-z0-9_.]\+\s/s sturnus.Chr28 /g' $f2 | sponge $f2 # no difference

# ordering?
maf_order $f1 melanotis ficedula | sponge $f1
maf_order $f2  melanotis sturnus | sponge $f2

#f1=ordered.$f1
#f2=ordered.$f2

export path="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/polarisation/tba_pipeline"

# the orginal fasta are needed so this command will fech unique scaffold names in maf
grep "s " $f1 |tr -s " " | cut -f 2 -d " " | sort | uniq | cut -f 2 -d "." > f1fa
grep "s " $f2|tr -s " " | cut -f 2 -d " " | sort | uniq | cut -f 2 -d "." > f2fa

# and download to scratch
parallel "find ${path}/sturnus_split/ ${path}/ficedula_split/ ${path}/melanotis_split/ -name '*{}.*' -exec cp \{\} $SCRATCHDIR/ \\;" :::: f1fa
parallel "find ${path}/sturnus_split/ ${path}/ficedula_split/ ${path}/melanotis_split/ -name '*{}.*' -exec cp  \{\} $SCRATCHDIR/ \\;" :::: f2fa
echo "fasta files copied:"
ls *fa
# the tba pipeline has some naming constrains that fasta has to be named by the species without suffix
cat melanotis*fa > melanotis
cat sturnus*fa > sturnus
cat ficedula*fa > ficedula

# actually I forgot to rename the files before alignments in melanotis and sturnus, but you can keep the renaming of fasta sequences in here
module add python27-modules-gcc 

python ~/scripts/tbarename_fasta.py --fasta melanotis
python ~/scripts/tbarename_fasta.py --fasta sturnus 
mv melanotis.lenad melanotis
mv sturnus.lenad sturnus
echo "fasta renamed:"
ls

f1s=${f1/.maf/.sing.maf} # the sing extension is another caprice of TBA
f2s=${f2/.maf/.sing.maf}
f1s=`echo $f1s | sed 's/melanotis_[a-zA-Z0-9]*_ficedula/melanotis.ficedula/g'`
f2s=`echo $f2s | sed 's/melanotis_[a-zA-Z0-9]*_sturnus/melanotis.sturnus/g'`
# single cov
single_cov2 $f1 F=${f1/.maf/.single_cov_deleted.maf} > $f1s
single_cov2 $f2 F=${f2/.maf/.single_cov_deleted.maf} > $f2s
#sed 's/##eof maf//g' $f1s | sponge $f1s # maybe that is another problem? 
#sed 's/##eof maf//g' $f2s | sponge $f2s 

#output all blocks.
out=${f1/.maf/}_${f2}

tba + "((melanotis sturnus) ficedula)" $f1s $f2s $out 

wc -l *maf

# get the lengths of alignments
printf "\n"
grep "s melanotis" $f1 | tr -s " " | cut -f 4 -d " " | awk -v n="$f1" '{s+=$1}END{print s, n}'
grep "s melanotis" $f2 | tr -s " " | cut -f 4 -d " " | awk -v n="$f2" '{s+=$1}END{print s, n}'
grep "s melanotis" $f1s | tr -s " " | cut -f 4 -d " " | awk -v n="$f1s" '{s+=$1}END{print s, n}'
grep "s melanotis" $f2s | tr -s " " | cut -f 4 -d " " | awk -v n="$f2s" '{s+=$1}END{print s, n}'

grep "s melanotis" $out | tr -s " " | cut -f 4 -d " " | awk -v n="$out" '{s+=$1}END{print s, n}'

####now parse maf to bed with WGAbed https://henryjuho.github.io/WGAbed/
# wga needs the chromom of the reference

refchr=`grep "s melanotis" $f2|tr -s " " | cut -f 2 -d " " | sort | uniq | cut -f 2 -d "."` 

maf_to_bed.py -i $out -r melanotis -c $refchr > ${out/maf/bed}

# I will have to remove indels as Ben is doing moreover check the outfile of this 


# do not forget to remove auxiliary files
rm $f1
rm $f2
rm *fa
#rm melanotis
#rm sturnus
#rm ficedula
mkdir $outd
mv * $outd/
cp -r $outd $wd/ || export CLEAN_SCRATCH=false

