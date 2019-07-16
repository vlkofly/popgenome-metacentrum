#!/bin/bash
#PBS -N polarise 
#PBS -l select=1:ncpus=1:mem=100gb:scratch_local=100gb
#PBS -l walltime=4:00:00
#PBS -j oe

# supply argument chr=Chr1 & config & outd 
# polarise per chromosome, because to put whole genome on one pile would be crazy
# or do it in one piece?
# example of running the script qsub -v 'chr=Chr1' ~/scripts/polarisation3_est_sfs.sh

# you have to modify on line 80 -100 how to parse the outgroups
module add gsl-2.1-gcc
module add bedtools-2.26.0
export CPATH=/software/gsl/2.1/include/:$CPATH
export LIBRARY_PATH=/software/gsl/2.1/lib/:$LIBRARY_PATH

wd=$PBS_O_WORKDIR
bin="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/est-sfs-release-2.03/est-sfs"
beddir="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/polarisation/tba_pipeline/multiz/multiz_results/final_beds" # the output from tba pipeline
#focfreq="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/mapped_to_melanotis/snp_results/focal_freq_sfs_renamed.txt" # filtered vcf
#focfreq="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/polarisation/est-sfs/polyglottos_outgroup/focal_freq_sfs_polyglottos_est-sfs_input_renamed.txt"
###check if the optional parameters were set, if not set default
if [ -z "$outd" ]; then
        echo "outidr not specified usig default est_sfs_results"
	outd="est_sfs_results"
fi

if [ -z "$config" ]; then
        echo "configuration file not specified usig default config_est_sfs_kimura.txt"
	config="config_est_sfs_kimura.txt"
fi

if [ -z "$focfreq" ]; then
       echo "You need to supply focal freq file if you do not want default:"
       focfreq=$wd/focal_freq_sfs.txt
fi



# the vcf was processed this way to match est-sfs format which is:
# 20,0,0,0        0,0,0,1 0,0,0,1 0,0,0,1 # ACGT counts
# vcftools --vcf ROH/snp_gala.masked.mindp7.0.75missing.vcf --max-missing 1 --counts
# sed 's/:/\t/g' out.frq.count | tail -n +2 | awk '{a=0;t=0;g=0;c=0;if ($5 == "A") a=$6; if ($5 == "T") t=$6; \
# if ($5 == "G") g=$6; if ($5 == "C") c=$6;if ($7 == "A") a=$8; if ($7 == "T") t=$8; if ($7 == "G") g=$8; \
# if ($7 == "C") c=$8; printf "%s\t%d\t%d\t%s,%s,%s,%s\n" ,$1,$2-1,$2,a,c,g,t}' | sed 's/Mmel_//g' > focal_freq_sfs_renamed.txt

# for polyglottos the input was generated from vcf file: /storage/brno3-cerit/home/vlkofly/mockingbird_genomics/mapped_to_melanotis/vcf_melanotis/all_polyglot/all_polyglottos.vcf.gz
# with use of two parsing scripts: in folder /polarisation/est-sfs/polyglottos_outgroup/
# parse_polyglottos.sh intersect_polyglottos.sh  


trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi

cp  $beddir/*_${chr}_*.bed $SCRATCHDIR/
cp  $focfreq $SCRATCHDIR/
cp $wd/$config $SCRATCHDIR/ # est_sfs config file should be in the working directory
cp $wd/seedfile.txt $SCRATCHDIR/ # the same as seedfile.txt
cd $SCRATCHDIR

bed=`ls *bed`
focfreq=`basename $focfreq`

# prepare focfreq (get only sites for a given chromosome?) test it on one chromosome and then do it for all

grep -P "${chr}\t" $focfreq > ${focfreq/.txt/${chr}.txt} # pearl grep expr 
focfreq=${focfreq/.txt/${chr}.txt}

# than you can probably run for loop for all chromosomes and do the polarisation alltogether, with chrom28 it runs just 3 mins

# remove indels in the bed alignment when sorting do not forget to redirect tmp
awk '{if ($8 ~ /.,.,./) print $0}' $bed | sort -T $SCRATCHDIR -k1,1 -k2,2n > ${bed/.bed/noindel.bed}
bed=${bed/.bed/noindel.bed}
# let's assume the both beds are sorted and do the intersection when you want to keep all sites that are in vcf
# because est-sfs allows missingness in outgroups

intersectBed -sorted -loj -a $focfreq -b $bed > ${chr}.intersected

# parse the intersected file so it conforms to the format of est-sfs
# get the column with outgroup values corresponds to melanotis,ficedula,sturnus

### this option for outgroup melanotis,sturnus,ficedula
#paste <(cut -f 4 ${chr}.intersected) <(cut -f 12 ${chr}.intersected | awk -F, '{am=0;tm=0;gm=0;cm=0;as=0;ts=0;gs=0;cs=0;af=0;tf=0;gf=0;cf=0; \
#if($1 ~ /[aA]/)  am=1;if($1 ~ /[tT]/) tm=1;if($1 ~ /[gG]/) gm=1;if($1 ~ /[cC]/) cm=1; if($3 ~ /[aA]/) as=1; if($3 ~ /[tT]/) ts=1; \
#if($3 ~ /[gG]/) gs=1;if($3 ~ /[cC]/) cs=1;if($2 ~ /[aA]/) af=1;if($2 ~ /[tT]/) tf=1;if($2 ~ /[gG]/) gf=1;if($2 ~ /[cC]/) cf=1; \
#printf "%d,%d,%d,%d %d,%d,%d,%d %d,%d,%d,%d\n", am,cm,gm,tm,as,cs,gs,ts,af,cf,gf,tf }') > ${chr}.est_sfs.input
#missing=`grep "0,0,0,0 0,0,0,0 0,0,0,0" ${chr}.est_sfs.input | wc -l`

### this option for sturnus and ficedula:
paste <(cut -f 4 ${chr}.intersected) <(cut -f 12 ${chr}.intersected | awk -F, '{am=0;tm=0;gm=0;cm=0;as=0;ts=0;gs=0;cs=0;af=0;tf=0;gf=0;cf=0; \
if($1 ~ /[aA]/)  am=1;if($1 ~ /[tT]/) tm=1;if($1 ~ /[gG]/) gm=1;if($1 ~ /[cC]/) cm=1; if($3 ~ /[aA]/) as=1; if($3 ~ /[tT]/) ts=1; \
if($3 ~ /[gG]/) gs=1;if($3 ~ /[cC]/) cs=1;if($2 ~ /[aA]/) af=1;if($2 ~ /[tT]/) tf=1;if($2 ~ /[gG]/) gf=1;if($2 ~ /[cC]/) cf=1; \
printf "%d,%d,%d,%d %d,%d,%d,%d\n", as,cs,gs,ts,af,cf,gf,tf }') > ${chr}.est_sfs.input

### this option for polyglottos,sturnus,ficedula
# one more column was added
#paste <(cut -f 4,5 ${chr}.intersected) <(cut -f 13 ${chr}.intersected | awk -F, '{am=0;tm=0;gm=0;cm=0;as=0;ts=0;gs=0;cs=0;af=0;tf=0;gf=0;cf=0; \
#if($1 ~ /[aA]/)  am=1;if($1 ~ /[tT]/) tm=1;if($1 ~ /[gG]/) gm=1;if($1 ~ /[cC]/) cm=1; if($3 ~ /[aA]/) as=1; if($3 ~ /[tT]/) ts=1; \
#if($3 ~ /[gG]/) gs=1;if($3 ~ /[cC]/) cs=1;if($2 ~ /[aA]/) af=1;if($2 ~ /[tT]/) tf=1;if($2 ~ /[gG]/) gf=1;if($2 ~ /[cC]/) cf=1; \
#printf "%d,%d,%d,%d %d,%d,%d,%d\n", as,cs,gs,ts,af,cf,gf,tf }') > ${chr}.est_sfs.input


# how many sites are not covered by any outgroup?
missing=`grep "0,0,0,0 0,0,0,0 0,0,0,0" ${chr}.est_sfs.input | wc -l`

# now run the est-sfs
$bin $config ${chr}.est_sfs.input seedfile.txt ${chr}.sfs.output.txt ${chr}.site.output.txt


# parse the output get sites to be switched the same as polarisation key that I have for lyrata
# check the distribution of probabilities and based on that derive thresholds
# there will be some sites that will be clear adepts for repolarisation, approx when the prob (column 3 in sfs output) of major allele being ancestral will be <0.4
# than some gray zone in betwee 0.4-0.6 sites that cannot be clearly polarised and then >0.6 sites where major is ancestral (that will be majority of sites)

#paste <(cut -f 1-3,13 ${chr}.intersected) ${chr}.est_sfs.input <(tail -n +9 ${chr}.site.output.txt) |  tr -s " " \\t  > ${chr}.site.complete.output.bed #to assign the prob estimates with chrom positions
paste <(cut -f 1-3,12 ${chr}.intersected) ${chr}.est_sfs.input <(tail -n +9 ${chr}.site.output.txt) |  tr -s " " \\t  > ${chr}.site.complete.output.bed # in the case without polyglottos as outgroup

# add header in the Rscript


# check also if you change the outgroups what will happen with sfs and polarisation key

cat ${chr}.sfs.output.txt | tr "," \\n | sponge ${chr}.sfs.output.txt

# now write an R script that would summarise the polarisation and pick the polarisation key
Rscript ~/scripts/plot_polarisation.R ${chr}.site.complete.output.bed

# I will need to annotate alleles in vcf which is derived and which ancestral.


wc -l *
echo "while there is $missing sites without reference"
#rm *bed
#rm focal_freq*
mkdir $outd
mv * $outd
cp -r $outd $wd/ || export CLEAN_SCRATCH=false

