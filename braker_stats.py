#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 18:03:40 2021
@author: ce9
"""


import sys

gtf_file=sys.argv[1]
gff_file=sys.argv[2]  
ID=sys.argv[3]

# read in gtf file(augustus.final.hints.gtf)
with open(gtf_file, "r") as gtf_in:
    c=gtf_in.readlines()
    

#Get n genes
s="".join(c).replace("gene_id", "Gid")
e=s.count("gene")    

print("N genes:" + str(e))


# Get all ntron and exon lens for all transcripts
intron_lens=[]
exon_lens=[]

for h in c:    
    h=h.split()
    
    if "intron" in h:
        intronL=int(h[4]) - int(h[3])
        intron_lens.append(intronL)
    if "exon" in h:
        exonL=int(h[4]) - int(h[3])
        exon_lens.append(exonL)
        

# write out intron and exon lens 
with open(str(ID) + '_intron_lens.txt', 'w') as f:
    for item1 in intron_lens:
        f.write("%s\n" % item1)        
        
with open(str(ID) + '_exon_lens.txt', 'w') as f1:
    for item2 in exon_lens:
        f1.write("%s\n" % item2)        


#Print out some stats...
print("Min exon len:" + str(min(exon_lens)))
print("Max exon len:" + str(max(exon_lens)))
print("Min intron len:" + str(min(intron_lens)))        
print("Max intron len:" + str(max(intron_lens)))     
print("N total introns:" + str(len(intron_lens)))        
print("N total exons:" + str(len(exon_lens)))        
        
Av_exon_len=int(sum(exon_lens) / len(exon_lens))
print("average exon length:" + str(Av_exon_len))

Av_intron_len=int(sum(intron_lens) / len(intron_lens))
print("average intron length:" + str(Av_intron_len))
        

#Read in gff file
with open(gff_file, "r") as gff_in:
    d=gff_in.readlines()
    
    
# Split by gene entry     
genes2="".join(d).split("# start gene")[1:]

# Find n start codons
start_codons=str(genes2).count("start_codon")

print("N start codons:" + str(start_codons))

# Find n stop codons
stop_codons=str(genes2).count("stop_codon")

print("N stop codons:" + str(stop_codons))


#Gen n single exon genes
singles=str(genes2).count("single")

print("N single exon transcripts:" + str(singles))


#Identify full length transcripts 
full_transcripts=[]

for n in genes2:
    n=n.split()
    if "start_codon" and "stop_codon" in n:
        full_transcripts.append(n[0])

    
# Find n full length transcripts       
len_full_trans=len(full_transcripts)

print("N transcripts with start and stop:" + str(len_full_trans))


# Get gene lens
gene_lens=[]


# Perc support for each trancript on each scaffold
perc_support=[]

for i in genes2:
    i = i.split("\n")
    k = "".join(i).split("% of transcript supported by hints (any source):")[1].split("#")[0]
    j=(i[1].split())
    perc_support.append(j[0] + " " + j[8] + " " + k)
    gene_len=int(j[4]) - int(j[3])
    gene_lens.append(gene_len)    
    
    
#write out stat files     
with open(str(ID) + '_gene_lens.txt', 'w') as f2:
    for item3 in gene_lens:
        f2.write("%s\n" % item3)     
        
        
with open(str(ID) + '_perc_support.txt', 'w') as f3:
    for item4 in perc_support:
        f3.write("%s\n" % item4)            

#Get average gene length
Av_gene_len=int(sum(gene_lens) / len(gene_lens))

print("average gene len:" + str(Av_gene_len))



# For testing:
#ID="test"
#gtf_file="/Users/ce9/Desktop/Annotation_stats/Ladybird1/augustus.hints.gtf"
#with open("/Users/ce9/Desktop/Annotation_stats/Ladybird2/augustus.hints.gtf", "r") as gtf_in:
#with open("/Users/ce9/Desktop/Annotation_stats/Ladybird2/augustus.hints.gff", "r") as gff_in:
#with open("/Users/ce9/Desktop/Annotation_stats/Ladybird1/augustus.hints.final.gff", "r") as gff_in:
#gff_file="/Users/ce9/Desktop/Annotation_stats/Ladybird1/augustus.hints.gff"
  
