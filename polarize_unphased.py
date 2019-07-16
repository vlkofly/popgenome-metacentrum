#!/usr/bin/env python3

### Description:
# written 23.10.2018 as a part of selscan pipeline by Jakub Vlcek
# take vcf after phasing and swap ALT and REF and genotypes there where minor allel is ancestral
# the information about which site to swap is supplied in a table poltab where format is CHROM\tPOS\n
# output is repolarised vcf with a string "polarised" added at the end of original vcf file name.

# original script in vi ~/../filip_kolar/scripts/polarize_phased.py now rewritten for the purpose of repolarisation of mockingbird vcf
# the difference is that mockingbirds have unphased vcf and the poltab is a bit different as a result of est-sfs polarisation pipeline
# Chr3  64015389    64015390    G
# the base in the third column is ancestral so you need to switch alleles in case where ancestral is in ALT column
# do not separate it per scaffold this time

# generate two output vcfs, one where ancestral allele will be only added to INFO files in format AA=X
# second vcf where ALT and REF is switched for the purpose of snpeff (will snpeff be ok with ref allele being something else than what is ref seq?)

### Import:
import argparse
import re
import os
import time
import gzip

### Arguments definition:
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', type=str, metavar='vcf', required=True, help='REQUIRED: Vcf file to repolarise, must be unzipped')
parser.add_argument('-poltab', type=str, metavar='poltab', required=True, help='REQUIRED: table defining which sites needs to be polarised ')
#parser.add_argument('-scf', type=str, metavar='scf', required=True, help='REQUIRED: scaffold id ')

args = parser.parse_args()
### generate list of sites that needs to be polarised
# list of lists
poltab = open(args.poltab,'r')
#site_list = set() # put it in a list, no set will be much more efficient
site_list = {} # let's try if dictionary will be faster?

log = open ("polarisation.log","w")

for ln,l in enumerate(poltab):
    l = l.strip("\n")
    chro,pos1,pos2,aa = l.split("\t") # list chr, position, AA
    chro="Mmel_"+chro   # filter can have only one iterable as input
    idx = pos2+chro
    row = (ln,chro,pos2,aa) # tuple
    site_list[idx] = row
    #site_list.add(row) # list of rows, better create set of rows for faster lookup

def getpos (chrom,pos):
    """function to get particular site from site_list based on chrom and position"""
    #res = list(filter(lambda row: row[1] == chrom and row[2] == pos,site_list) ) # 5 seconds per site
    #res = [x for x in site_list if x[1] == chrom and x[2] == pos] # still 4 seconds per site, lets try to merge site and chrom to one variable
    #res = [x for x in site_list if x[4] == pos+chrom] # no this does not help much
    #res = list((x for x in site_list if x[4] == pos+chrom)) # generator expression ... not really faster either ~3.9
    #res = site_list[pos+chrom] # pity you can not do it this way as there are missing sites
    res = site_list.get(pos+chrom,"NA") # oh god this is the thing ~0.0004s
    if res == "NA":
        log.write("site "+chrom+" "+pos+" not found\n")
        global notinpol
        notinpol+=1
        res = [chrom,pos,pos,"NA"]
    else:
        #res = res[0]
        #site_list.discard(res) # remove item from the set, so the function gets faster
        del site_list[pos+chrom]
    return res


#print(site_list[0:3])
print("length of polarisation list: "+str(len(site_list)))

repldict={}
repldict["0/0"] = "1/1"
repldict["0/1"] = "1/0"
repldict["1/1"] = "0/0"
repldict["./."] = "./."

print (repldict) # will have to add missing genotype as well
### read through vcf
if args.vcf.endswith(".gz"):
    vcf = gzip.open(args.vcf,'r')
else:
    vcf = open(args.vcf,'r')

vcfbase = os.path.basename(args.vcf)
vcf_new_name = vcfbase.replace(".vcf",".polarised.vcf")
vcf_annotate_name = vcfbase.replace(".vcf",".ancestral.annotated.vcf")

print(vcfbase)
anovcf = open (vcf_annotate_name,"w")
with open (vcf_new_name,'w') as out:
    n = 0 # n is number of switched sites
    nt = 0
    notpol = 0
    notinpol = 0
    for l in vcf:
        t1 = time.time()
        if l.startswith("#"):
            out.write(l)
            anovcf.write(l) # you will have to add the annotation to the header
        else:
            nt += 1
            ls = l.strip("\n")
            #print(l.split("\t")[0:9])
            CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT = ls.split("\t")[0:9]
            GT=list(ls.split("\t")[9:])
                   

            polar = getpos(CHROM,POS) # for every position in vcf look up this position in polarisation  
            #print(polar)
            # generate annotated vcf
            INFOann = INFO+";AA="+polar[3]
            anovcf.write(CHROM+'\t'+POS+'\t'+ID+'\t'+REF+'\t'+ALT+'\t'+QUAL+'\t'+FILTER+'\t'+INFOann+'\t'+FORMAT)
            for g in GT:
                anovcf.write('\t'+g)
            anovcf.write('\n')


            if polar[3] == ALT: # if alternative allele is the ancestral switch it
                n+=1
                #print ("repolarising position %s, swap %s as reference\n" % (POS, ALT))
                # solve the genotypes with some replacement dict
                # recalculate AC & AF in info AC = AN - AC ;AF = 1-AF

                #print("former info: "+INFO+"\n")
                ac,af,an = INFO.split(";")[0:3]
                ac = int(ac.replace("AC=",""))
                af = float(af.replace("AF=",""))
                an = int(an.replace("AN=",""))

                ac = an - ac
                af = 1 - af
                acstr = "AC="+str(ac)
                afstr = "AF="+"{:.2f}".format(af)

                INFO=re.sub("AC=[0-9]+",acstr,INFO)
                INFO=re.sub("AF=0\.[0-9]+",afstr,INFO)
                #print("new infor: "+INFO+"\n")
                
                out.write(CHROM+'\t'+POS+'\t'+ID+'\t'+ALT+'\t'+REF+'\t'+QUAL+'\t'+FILTER+'\t'+INFO+'\t'+FORMAT)
                for g in GT:
                    genrow = g.split(":")
                    gn = repldict[genrow[0]] # replace genotype
                    newgen = genrow[1:] # remove old genotype
                    newgen.insert(0,gn) # replace with new
                    #reorder depths or, no I do not care about depths, but keep in mind that depths are switched for repolarised sites
                    newgenstr = ":".join(newgen)

                    out.write('\t'+newgenstr)
                    #print ("genotype %s swapped to %s" % (g,newgenstr))
                out.write('\n')
            elif polar[3] == "NA":
                # do not output sites that could not be polarised
                notpol+=1
            else:
                out.write(l)
                log.write("position "+CHROM+" "+POS+" not polarised\n" )
        t2 = time.time()
        #print ('time per position: %f' % (t2-t1))
print ("repolarised %d out of %d and removed %d positions " % (n,nt,notpol))                
anovcf.close()

