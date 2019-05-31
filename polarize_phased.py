#!/usr/bin/env python3

### Description:
# written 23.10.2018 as a part of selscan pipeline by Jakub Vlcek
# take vcf after phasing and swap ALT and REF and genotypes there where minor allel is ancestral
# the information about which site to swap is supplied in a table poltab where format is CHROM\tPOS\n
# output is repolarised vcf with a string "polarised" added at the end of original vcf file name. 
### Import:
import argparse

### Arguments definition:
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', type=str, metavar='vcf', required=True, help='REQUIRED: Vcf file to repolarise ')
parser.add_argument('-poltab', type=str, metavar='poltab', required=True, help='REQUIRED: table defining which sites needs to be polarised ')
parser.add_argument('-scf', type=str, metavar='scf', required=True, help='REQUIRED: scaffold id ')

args = parser.parse_args()
### generate list of sites that needs to be polarised

poltab = open(args.poltab,'r')
site_list = []
for l in poltab:
    l = l.strip("\n")
    
    scaffold,site = l.split("\t")
    #print(scaffold,site)
    if scaffold == args.scf:
        site_list.append(site)

print(len(site_list))

repldict={}
repldict["0|0"] = "1|1"
repldict["0|1"] = "1|0"
repldict["1|0"] = "0|1"
repldict["1|1"] = "0|0"
repldict[".|."] = ".|."

print (repldict) # will have to add missing genotype as well
### read through vcf 
vcf = open(args.vcf,'r')
vcf_new_name = args.vcf.replace(".vcf",".polarised.vcf")
with open (vcf_new_name,'w') as out:
    n = 0
    nt = 0
    for l in vcf:
        if l.startswith("#"):
            out.write(l)
        else:
            nt += 1
            ls = l.strip("\n")
            #print(l.split("\t")[0:9])
            CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT = ls.split("\t")[0:9]
            GT=list(ls.split("\t")[9:])
            
            if POS in site_list:
                n+=1
                print ("repolarising position %s, swap %s as reference\n" % (POS, ALT))
                # solve the genotypes with some replacement dict

                out.write(CHROM+'\t'+POS+'\t'+ID+'\t'+ALT+'\t'+REF+'\t'+QUAL+'\t'+FILTER+'\t'+INFO+'\t'+FORMAT)
                for g in GT:
                   
                    gn = repldict[g]
                    out.write('\t'+gn)
                    print ("genotype %s swapped to %s" % (g,gn))
                out.write('\n')
            else:
                out.write(l)
print ("repolarised %d positions out of %d" % (n,nt))                









  



