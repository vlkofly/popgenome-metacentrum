# popgenome-metacentrum

**Miscellaneous scripts dealing with various problems of population genomics.**
The majority of bash scripts contain job resource section that is specific for job scheduling system [PBS Professional](https://www.pbsworks.com/SupportGT.aspx?d=PBS-Professional,-Documentation)
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
R script that follows extraction of ROH by previous script
