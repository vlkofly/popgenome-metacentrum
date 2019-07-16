#!/bin/bash
#PBS -N twisst 
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=6:mem=30gb:scratch_local=80gb
#PBS -j oe

# this is a wrapper script written for Metacentrum that wraps TWISST analysis of Simon Martin
# Written on May 3 2019 for the purpose of analysing Galapagos mockingbird genomics
# The command line variable will be vcf (I will run parallel for each chromosomes)
# The vcf must be phased and in this case it is result of script ~/../filip_kolar/scripts/shapeit.sh 

# Inscript variables section:
nt='6'
wd=$PBS_O_WORKDIR
parsevcf="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/bin/parseVCF.py" # binary dir
phyml="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/genomics_general/phylo/phyml_sliding_windows.py" # binary dir
twisst="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/twisst/run_twisst_parallel.py" # binary dir
phyml_bin="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/phyml/src/phyml"
groups=${wd}/gala.groups.txt # taxon groups must be present at working directory
#outdir=twisst_results/shorter_snp
#outdir=windsize100chr9.method.complete.diffinputgroups # outdir as online variable

echo "output goes into $outdir"

### load modules section:
module add python27-modules-gcc # the scripts needs numpy and others
module add phyml-3.0
module add ete-toolkit-3.3.0

### copy files to scratch section:
trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
cp $wd/$vcf $SCRATCHDIR # copy vcf file to scratch 
# to have the script super kosher you should upload the programs to scratch as well
cp $parsevcf $SCRATCHDIR
cp $phyml $SCRATCHDIR
cp $twisst $SCRATCHDIR
cp $groups $SCRATCHDIR
cp /storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/genomics_general/genomics.py $SCRATCHDIR ####in phyml you have to load the genomics.py
cp /storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/twisst/twisst.py $SCRATCHDIR
cp /storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/twisst/plot_twisst.R $SCRATCHDIR
#cp /storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/twisst/example_plot.R $SCRATCHDIR




parsevcf=`basename $parsevcf`
phyml=`basename $phyml`
twisst=`basename $twisst`
groups=`basename $groups`
cd $SCRATCHDIR
echo "files copied at `date` :"
ls
vcf=`basename $vcf`
geno=${vcf/vcf.gz/geno.gz}
base=${vcf/.vcf.gz/}

###commands begin here
# parse VCF to geno 
python $parsevcf -i $vcf | gzip > $geno  
echo "vcf turned into geno at `date`"
ls

# do sliding windows trees
# there are many parameter to play with... lets do it for 50 snps per window
#default command all individuals
#python $phyml --phyml $phyml_bin -p $base -g $geno -T $nt --windType sites --verbose --windSize 100  -Mi 50  --tmp . --model GTR --log ${base}.phyml.log #the size is in bp so you have to play with --verbose 
#--outgroup C_270_A

#do the analysis only for a subset of individuals DWI
#--outgroup F_272,F_280,F_276
python $phyml --phyml $phyml_bin -p $base -g $geno -T $nt --windType sites  \
       --verbose --windSize 100  -Mi 50  --tmp . --model GTR --log ${base}.phyml.log \
       --individuals F_272,F_280,F_276,I_364,I_370,I_380,W_025,W_001,W_009,D_001,D_002,D_013

echo "phyml finished at `date`"
# run twisst
# group definition variable:
#gr="-g C C_254_A,C_254_B,C_256_A,C_256_B,C_270_A,C_270_B -g F F_272_A,F_272_B,F_276_A,F_276_B,F_280_A,F_280_B -g I I_364_A,I_364_B,I_370_A,I_370_B,I_380_A,I_380_B -g M M_182_A,M_182_B,M_192_A,M_192_B,M_208_A,M_208_B -g P P_001_A,P_001_B,P_016_A,P_016_B,P_017_A,P_017_B -g S S_001_A,S_001_B,S_005_A,S_005_B,S_140_A,S_140_B -g W W_001_A,W_001_B,W_009_A,W_009_B,W_025_A,W_025_B -g D D_001_A,D_001_B,D_002_A,D_002_B,D_013_A,D_013_B "
#gr="-g I I_364_A,I_364_B,I_370_A,I_370_B,I_380_A,I_380_B -g D D_001_A,D_001_B,D_002_A,D_002_B,D_013_A,D_013_B -g W W_001_A,W_001_B,W_009_A,W_009_B,W_025_A,W_025_B -g M M_182_A,M_182_B,M_192_A,M_192_B,M_208_A,M_208_B -g P P_001_A,P_001_B,P_016_A,P_016_B,P_017_A,P_017_B -g S S_001_A,S_001_B,S_005_A,S_005_B,S_140_A,S_140_B -g C C_254_A,C_254_B,C_256_A,C_256_B,C_270_A,C_270_B -g F F_272_A,F_272_B,F_276_A,F_276_B,F_280_A,F_280_B "

#for subset
gr="-g I I_364_A,I_364_B,I_370_A,I_370_B,I_380_A,I_380_B -g D D_001_A,D_001_B,D_002_A,D_002_B,D_013_A,D_013_B -g W W_001_A,W_001_B,W_009_A,W_009_B,W_025_A,W_025_B -g F F_272_A,F_272_B,F_276_A,F_276_B,F_280_A,F_280_B "

#python $twisst --threads $nt -t ${base}.trees.gz -w ${base}.weights.csv.gz --outgroup C_270_A --outputTopos ${base}.topologies.trees -g C -g D -g W -g F -g I -g M -g P -g S --method fixed --groupsFile $groups 
python $twisst --threads $nt -t ${base}.trees.gz -w ${base}.weights.csv.gz $gr --method complete  

gunzip ${base}.weights.csv.gz
### plot the tree summary
# the function source script had to be modified to plot the summary nicely.
# generate R script, to get a clue see example_plot.R:

cat << EOF > script.R
source ("plot_twisst.R")

weights_file <- "${base}.weights.csv"
window_data_file <- "${base}.data.tsv" 
twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)
ss<-function(x) {
        N=x
        best.topos<-head(order(twisst_data$weights_overall_mean,decreasing = T),n = N)
        twisst_data$weights_overall_mean<-twisst_data$weights_overall_mean[best.topos]
        twisst_data$topos<-twisst_data$topos[best.topos]
        # you have to use lapply to subset the dataset list of lists madness
        twisst_data$weights<-lapply(twisst_data$weights,function(x) x[best.topos])
        twisst_data$weights_raw<-lapply(twisst_data$weights_raw,function(x) x[best.topos])
        twisst_data$weights_mean<-lapply(twisst_data$weights_mean,function(x) x[best.topos])
        return (twisst_data)
}

svg("awerage.weight.${base}.svg",pointsize=18)
plot.twisst.summary(twisst_data, cex=0.7, y_scale=0.095,x_scale = 0.1 )
dev.off()

twisst_data_smooth <- smooth.twisst(twisst_data,span=0.05, span_bp =NULL , spacing = NULL)

svg("twisst.smooth.${base}.svg",width=20,height=10,pointsize=20)
plot.twisst(twisst_data_smooth, mode=2,show_topos=T) #mode 2 overlays polygons, mode 3 would stack them
dev.off()



EOF

Rscript script.R

#rm $fgn
#rm $sgn

mkdir $outdir

mv * $outdir/
cp -r $outdir $wd/twisst_results/  || export CLEAN_SCRATCH=false
