#!/bin/bash
#PBS -N prepare_selscan
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=1:mem=120gb:scratch_local=300gb
#PBS -j oe

# follow up script after shapeit that prepares inputs for selscan 

echo "started at `date`"

#inscript variables definition segment:
wd=$PBS_O_WORKDIR
outdir=selscan_prepared
#vcf="mocker_to_toxo_1.merged.vcf.gz" # supply vcf as an argument
#scf #  main variable argument integer only
ident=scaffold_${scf}
gmap=/storage/brno3-cerit/home/filip_kolar/programs/predictGMAP/src/predictGMAP # program to generate genetic map https://github.com/szpiech/predictGMAP
popmap=/storage/brno3-cerit/home/filip_kolar/phasing_halleri/popmap.txt # one pop every row separated by tab from comma separated list of samples POP	s1,s2,s3

module add bcftools-1.6
module add htslib-1.6

trap 'clean_scratch' TERM EXIT

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi

cp ${wd}/$vcf $SCRATCHDIR
cp /storage/brno3-cerit/home/filip_kolar/phasing_arenosa/rec_map/gmap_reference.txt $SCRATCHDIR  #copy genetic mapp_reference.txt
cp $popmap $SCRATCHDIR  #copy genetic map
cp /storage/brno3-cerit/home/filip_kolar/scripts/polarize_phased.py $SCRATCHDIR
cp /storage/brno3-cerit/home/filip_kolar/phasing_arenosa/polarisation/repolarized.lookupKey.perSpeciesThreshold.txt $SCRATCHDIR

cd $SCRATCHDIR


echo "data imported succesfully at `date`:"
ls

# check if vcf is gunzipped and in that case unzip


vcf=`basename $vcf`

if [[ $vcf == *.gz ]]; then
	gunzip $vcf
	vcf=${vcf/.gz/}
fi

## polarize
python3 polarize_phased.py -vcf $vcf -poltab repolarized.lookupKey.perSpeciesThreshold.txt -scf $ident > ${ident}.repolarisation.log

echo "repolarisation done at `date`"
## remove missigness in complete dataset or after population split?


## genetic map
# get query sites, extract sites from vcf
grep -v "#" $vcf | cut -f 2 > query.sites.txt 
grep chr${scf} gmap_reference.txt | sponge gmap_reference.txt # select only a given chrom map

$gmap --ref gmap_reference.txt --query query.sites.txt --max-gap 4000000 --out selscan_map_${ident}.txt

echo "map done at `date`"
#####separate by pop and parse vcf to shapeit format, it assumes 8 inds per population, but should not be hard to rewrite it
bgzip ${vcf/.vcf/}.polarised.vcf
tabix -p vcf ${vcf/.vcf/}.polarised.vcf.gz
popmap=`basename $popmap`
while read pop samp
do
pvcf=${vcf/.vcf/}.${pop}.vcf
pvcf_fin=${vcf/.vcf/}.${pop}.final.vcf

echo $pop
echo $samp

bcftools view -Ov -s $samp --trim-alt-alleles --min-ac 1 -e 'AF=1' -o $pvcf ${vcf/.vcf/}.polarised.vcf.gz


#grep "#" $pvcf > vcf.header
#grep -v "#" $pvcf | awk '{print gsub("scaffold_","",$1),$2, "chr"$1":"$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' > vcf.body.txt
#cat vcf.header vcf.body.txt > $pvcf_fin 
#parse vcf format for shapeit that there are snp ids and chrom is integer, maybe not neccessary, but have to keep the consistency with map.
# actually I without this and it works as well
done < $popmap




mkdir $outdir 
mv * $outdir/

cp -r $outdir ${wd}/ || export CLEAN_SCRATCH=false

#cp $SCRATCHDIR/* $wd/filtervcf_toxo/ || export CLEAN_SCRATCH=false
echo "prep of $scf done and results exported at `date`"

