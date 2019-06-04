# popgenome-metacentrum

**Miscellaneous scripts dealing with various problems of population genomics.**


The majority of bash scripts contain job resource section that is specific for job scheduling system [PBS Professional](https://www.pbsworks.com/SupportGT.aspx?d=PBS-Professional,-Documentation).
On servers with this system scripts can be run with qsub (e.g `qsub ROH_bcftools.sh`).
The scripts were written for running jobs on [MetaCentrum - Czech computing infrastructure](https://metavo.metacentrum.cz/en/about/index.html).
Thus some details (especially program locations) can cause bugs when running on different type of servers. 
If you want to run the scripts locally remove the job specification section (I will try to identify it with comments).


The scripts are far away from being absolutely general although major variables can be specified either inscript or with arguments.
It means you will always need to read and understand the code a bit in order to get reliable results.
Often I use definition of arguments via qsub (e.g. `qsub -v 'vcf=blabla.vcf.gz' script.sh` this defines variable *$vcf* that is used within script.sh)


## ROH_bcftools.sh 
**Filter vcf and extract runs of homozygosity with bcftools**

`qsub -v 'dp=black.mask.bed,vcf=vcf.gz' ROH_bcftools.sh`

## plot_ROH.R 
**Plot and analyze runs of homozygosity**

R script that follows extraction of ROH by previous script

The ROH file must be parsed with in order to contain only regions with `grep RG vcf.gz.roh >ROH.rg`
 
`Rscript ROH.rg`

## msmc0_split_scaff_depth_filt.sh
**Prepare multihetsep for msmc analysis from vcf - split scaffolds and filter missingness**
This program wraps steps that are necessary to take before running msmc on data that are originally in vcf format.
Input vcf must contain both snps and invariant sites because the script uses [vcfAllSiteParser.py](https://github.com/stschiff/msmc-tools)

## msmc2_msmc.sh
**Run msmc**
[Using code developed by Stephan Schiffels](https://github.com/stschiff/msmc2)
## msmc3_plot.R
**Plot msmc results**

## shapeit.sh
**Phase vcf by shapeit**

Phasing of vcf = linking snps into haplotypes A/C ---> A|C 
It makes use of reads aligned to reference that contain two snps of known phase = PIR = **P**hase **I**nformative **R**eads
The method is implemented in program [shapeit.v2](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html) developed by O. Delaneau.
My script only wraps up the program commands and facilitates the analysis.   

Example of qsub command to run this script for multiple chromosomes defined in chrlist.txt:

`parallel -a chrlist.txt "qsub -v 'vcf=snp_gala.masked.mindp7.0.75missing.vcf.gz,scf={},outdir=shapeit_result' -N shapeit.{} ~/scripts/shapeit.sh"`

There are again inscript variables that need to be modified to match your situation. So read the script and try to understand what does it do before you run it ;).
Particularly it needs following resources:
- genetic map (variable recrate)
- list of sample names with corresponding full path to bam alignment (a2th2_01       /storage/plzen1/home/vlkofly/rgadd2/a2th2_01.final.bam)
- program binaries

## prepare_selscan.sh
**Polarize phased vcf and create genetic map**
This script takes phased vcf from shapeit and returns polarised vcf split by population that will go as an input to selscan.
Example of qsub command to run this:
`qsub -v 'vcf=shapeit_result/scaffold_1.phased.vcf.gz,scf=1' ~/scripts/prepare_selscan.sh`
In-house polarisation python script `polarize_phased.py` needs a specific polarisation key - a list of sites.
Another input is genetic map from which genetic distances are inferred for each site in vcf by [predictGMAP](https://github.com/szpiech/predictGMAP)
The vcf has to be also split by population and population map must be supplied.

## selscan.sh
**Identify sites under selection by selscan**
Just a wrapper of [selscan](https://github.com/szpiech/selscan)
Example of qsub command to run this:
qsub -v 'vcf=scaffold_2.phased.DRG.vcf' ~/scripts/selscan.sh


## run_twisst.sh

