'''
    this script will map peptides to genome and output an peptide gff3 file with coordinates.
    which genome it maps to depends on the GTF input file provided.
    written by Yafeng Zhu @ Karolinska Institutet.  Email: yafeng.zhu@ki.se
'''

import sys
import os
import getopt
import numpy as np
import re
from Bio import SeqIO

class EXON(object):
    def __init__(self,number=0,id=None,transcript=None,gene=None,chr=None,strand=None,start=0,end=0,length=0,trans_start=0,trans_end=0):
        self.id=id
        self.gene=gene
        self.transcript=transcript
        self.number=number
        self.start=start  #chromosome start coordinate
        self.end=end  #chromosome end coordinate
        self.strand=strand
        self.chr=chr
        self.trans_start=trans_start
        self.trans_end=trans_end
    def length(self):
        self.length=self.trans_end-self.trans_start+1


def cal_trans_pos(exon_list): # calculate transcript position of exon start & end, exon_list is a list of exon objects
    strand=exon_list[0].strand
    if strand=="+":
        new_exon_list=sorted(exon_list,key=lambda x:x.start)
    else:
        new_exon_list=sorted(exon_list,key=lambda x:x.start,reverse=True)
    
    sumExonlength=0
    for exon in new_exon_list:
        exon_length=exon.end-exon.start+1
        exon.trans_start=1+sumExonlength
        exon.trans_end=exon.trans_start+exon_length-1
        sumExonlength+=exon_length
        
    return new_exon_list

def get_pep_cor(exon_object_list,n1,n2): # return peptide's chromosome start and end cor given peptide's trans_start (n1) and trans_end (n2)
    pep_chr=""
    pep_strand=""
    pep_chr_start=0
    pep_chr_end=0
    pep_start_exon=0
    pep_end_exon=0
    pep_st_exon_id=""
    pep_ed_exon_id=""
    for i in range(len(exon_object_list)):
        exon=exon_object_list[i]
        if n1<=exon.trans_end and n1>=exon.trans_start:
            pep_chr=exon.chr
            pep_strand=exon.strand
            pep_start_exon=exon.number
            pep_st_exon_id=exon_id_dic[exon.transcript+"_exon"+exon.number]
            if pep_strand=='+':
                pep_chr_start=exon.start+(n1-exon.trans_start)
            else:
                pep_chr_end=exon.end-(n1-exon.trans_start)

        if n2<=exon.trans_end and n2>=exon.trans_start:
            pep_chr=exon.chr
            pep_strand=exon.strand
            pep_end_exon=exon.number
            pep_ed_exon_id=exon_id_dic[exon.transcript+"_exon"+exon.number]
            if pep_strand=='+':
                pep_chr_end=exon.start+(n2-exon.trans_start)
            else: # chr_cor of n2 is pep_chr_start
                pep_chr_start=exon.end-(n2-exon.trans_start)

    return pep_chr,pep_strand,pep_chr_start,pep_chr_end,pep_start_exon,pep_end_exon,pep_st_exon_id,pep_ed_exon_id

def parse_gtf(infile):
    dic={}
    with open(infile,"r") as infile_object:
        for line in infile_object:
            if line[0]!="#": # skip lines commmented out
                row=line.strip().split("\t")
                if row[0] not in chr_dic:
                    continue;
                if row[2]=="CDS":
                    attri_list=row[8].split(";")
                    transID=""
                    exon=EXON(start=int(row[3]),end=int(row[4]),chr=row[0],strand=row[6])
                    for attri in attri_list:
                        if "transcript_id" in attri:
                            transID=attri.strip().replace("transcript_id ","").replace('\"',"")
                            exon.transcript=transID
                        if "exon_number " in attri:
                            exon.number=attri.strip().replace("exon_number ","").replace('\"',"")
                    if transID not in dic:
                        dic[transID]=[exon]
                    else:
                        dic[transID].append(exon)
                elif row[2]=="exon":
                    attri_list=row[8].split(";")
                    transID=""
                    exonid=""
                    exon_number=""
                    for attri in attri_list:
                        if "transcript_id" in attri:
                            transID=attri.strip().replace("transcript_id ","").replace('\"',"")
                        if "exon_number " in attri:
                            exon_number=attri.strip().replace("exon_number ","").replace('\"',"")
                        if "exon_id " in attri:
                            exonid=attri.strip().replace("exon_id ","").replace('\"',"")

                    if transID+"_exon"+exon_number not in exon_id_dic:
                        exon_id_dic[transID+"_exon"+exon_number]=exonid

    return dic

################  Comand-line arguments ################


if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print "Warning! wrong command, please read the mannual in Readme.txt."
    print "Example: python map_peptide2genome.py --input input_filename --gtf Homo_sapiens.GRCh37.75.protein_coding.gtf --fasta Homo_sapiens.GRCh37.75.pep.all.fa  --IDmap Ensembl75_IDlist.txt --output output_filename"
    print "First three columns in the input file must be :  peptide, ENSP_IDs, ENSG_IDs.
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=',
                                                         'gtf=',
                                                         'fasta=',
                                                         'IDmap=',
                                                         'output='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--gtf':gtf_file=arg
        elif opt == '--fasta': fasta_file=arg
        elif opt == '--IDmap':IDmap_file=arg
        elif opt == '--output': output_file=arg
        else:
            print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

#restrict peptide mapping in the 24 chromosomes.
chr_dic={"MT":1,"X":1,"Y":1,"M":1}
for i in range(1,23):
    chr_dic[str(i)]=1


exon_id_dic={}

print "reading GTF input file"
feature_dic=parse_gtf(gtf_file)
print "number of unique transcripts in GTF file",len(feature_dic)

print "number of unique exon id in GTF file",len(exon_id_dic)

seq_dic = SeqIO.index(fasta_file,'fasta')
print "number of unique protein sequences in fasta file",len(seq_dic)

IDlist_input=open(IDmap_file,"r")
id_dic={}
gene_dic={}
for line in IDlist_input:
    row=line.strip().split("\t")
    if len(row)==3:
        ensg=row[0]
        enst=row[1]
        ensp=row[2]
        if ensp not in id_dic:
            id_dic[ensp]=enst
        if ensg not in gene_dic:
            gene_dic[ensg]=[ensp]
        else:
            gene_dic[ensg].append(ensp)

IDlist_input.close()
print "number of unique ENSP IDs in ID table",len(id_dic)

input=open(input_file,'r') # peptide table with two columns, peptide sequence in first column, protein ID in second column

header=input.readline().strip().split("\t")
newheader=header+["chr","start","end","strand","pep_start_exon","pep_end_exon","start_exon_id","end_exon_id","NumOf.coding_transcript.From.gene","NumOf.trans.mapped.to","transcripts_id_pep_mapped_to"]
output=open(output_file,'w')
output.write("\t".join(newheader)+"\n")

non_mapped_pep=0
count_pep=0
for line in input:
    count_pep+=1
    row=line.strip().split("\t")
    ensg=row[2].split(";")
    symbol=row[3].split(";")
    if len(set(ensg))>1: #
        print "skip peptides %s that mapped to multiple genes." % (row[0])
        continue
    peptide=re.sub("[\W\d]","",row[0].strip())
    ensp=row[1].split(";")[0]
    enst=id_dic[ensp]
    try:
        exons=feature_dic[enst]
    except KeyError:
        non_mapped_pep+=1
        continue;
    
    aa_seq=str(seq_dic[ensp].seq)
    pep_index=aa_seq.index(peptide)
    
    pep_trans_start=3*pep_index+1
    pep_trans_end=pep_trans_start+3*len(peptide)-1
    
    exons=cal_trans_pos(exons)
    
    #print pep_trans_start,pep_trans_end
    pep_chr,pep_strand,pep_chr_start,pep_chr_end,pep_start_exon,pep_end_exon,start_exon_id,end_exon_id=get_pep_cor(exons,pep_trans_start,pep_trans_end)
    
    
    #handle exceptions
    if pep_chr_start>pep_chr_end:
        non_mapped_pep+=1
        #print peptide,ensp,enst
        continue;
    if pep_chr_start<=0:
        non_mapped_pep+=1
        #print peptide,ensp,enst,pep_trans_start,pep_trans_end
        continue;
    
    #check how many coding variants annotated for this gene
    gene_id=ensg[0]
    var_list=set(gene_dic[gene_id])
    map_to_var=[]
    for var in var_list:
        if peptide in str(seq_dic[var].seq):
            map_to_var.append(var)

    row[2]=";".join(set(ensg))
    row[3]=";".join(set(symbol))
    insert_columns=map(str,[pep_chr,pep_chr_start,pep_chr_end,pep_strand,pep_start_exon,pep_end_exon,start_exon_id,end_exon_id,len(var_list),len(map_to_var),";".join(map_to_var)])
    newrow=row+insert_columns

    output.write("\t".join(newrow)+"\n")


input.close()
output.close()
print "total number of unique peptides",count_pep
print "total number of unmapped peptides",non_mapped_pep

