#! /bin/python3

# this script will take output of est-sfs, where each site was assigned probability of major allele being ancestral
# the output is a bed file for all sites with information about which allele is ancestral
# there was the issue with sites where only one allele is fixed. The probability was assigned one even thought branch b1 did not support the ancestral allele, the script contains fix of this issue.



import argparse
import operator
parser = argparse.ArgumentParser()
parser.add_argument('-sfsout', type=str, metavar='sfsout', required=True, help='REQUIRED: output from est-sfs pasted with the bed (complete.bed) ')
parser.add_argument('-output', type=str, metavar='out', required=True, help='REQUIRED: output file ')

args = parser.parse_args()
out = open(args.output,'w')
inp = open(args.sfsout,'r')
minor_ancestral=0
norepol=0
total=0
for l in inp:
    total+=1
    l = l.strip("\n").rstrip("\t")
    #print(l)
    chrom,pos1,pos2,rb,focalfr,outgfr1,outgrfr2,idx1,idx2,probanc,nodeprobA,nodeprobC,nodeprobG,nodeprobT=l.split("\t")
    # now you need to get which base is major
    frlist = focalfr.split(",")
    #print(frlist)
    d=dict() # dictionary of allele frequencies in focal species
    d["A"] = int(frlist[0])
    d["C"] = int(frlist[1])
    d["G"] = int(frlist[2])
    d["T"] = int(frlist[3])
    majorbase = max(sorted(d),key=lambda key: d[key]) # get major; in the case of tie the first one is major
    majorfreq = max(d.values())
    del d[majorbase]
    minorbase = max(d,key=lambda key: d[key])
    minorfreq = max(d.values())
    #print (majorbase,minorbase)
    #print (majorfreq,minorfreq)
    probanc = float(probanc)

    # what ajorfreq=minorfreq? if it is equal then the first base is taken as major
    # if major = 48 then there is no minor
    da = dict() # dictionary of node probabilities

    if majorfreq == 48: # fix the bug, take as ancestral base  with highest prob of branch b1
        da["A"] = float(nodeprobA)
        da["C"] = float(nodeprobC)
        da["G"] = float(nodeprobG)
        da["T"] = float(nodeprobT)
        AA= max(da,key=lambda key: da[key])
    elif probanc > 0.70:

        AA=majorbase

    
    elif probanc < 0.20:
        AA=minorbase
        minor_ancestral+=1
    else:
        AA="NA" # in the case where ancestrality could not beed assigned write NA
        norepol+=1

    #print(chrom, pos2,focalfr,AA,probanc,nodeprobA,nodeprobC,nodeprobG,nodeprobT)
    out.write(chrom+'\t'+pos1+'\t'+pos2+'\t'+AA+'\n')
print("input: "+args.sfsout+" minor ancestral: "+str(minor_ancestral)+" not polarised: "+str(norepol)+" total: "+str(total))
out.close()

