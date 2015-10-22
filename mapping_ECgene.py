import sys
import os
import getopt
import operator
from Bio import SeqIO
from collections import OrderedDict
from mapping_EVDB import ISOFORM,EXON,PEPTIDE


###################define main Function ###############################
def main(): #the main function for mapping peptides
    for i in range(len(peparray)):
        peptide=peparray[i]
        variantID=peptide.variants[0]
        variant=variant_dic[variantID]#variant here is an object, an instance of ISOFROM()
        n=variant.seq.find(peptide.seq)
        peptide.trans_start=1+n*3
        peptide.trans_end=peptide.trans_start+3*len(peptide.seq)-1
        peptide.chr=variant.chr
        peptide.strand=variant.strand
        
        exons=variant_dic[variantID].exon
        for exon in exons:
            #print exon.number,exon.start,exon.end,exon.trans_start,exon.trans_end
            if peptide.trans_start<=exon.trans_end and peptide.trans_start>=exon.trans_start:
                peptide.exon1=exon.number
                if peptide.strand=='+':
                    peptide.start=exon.start+(peptide.trans_start-exon.trans_start)+1
                else:
                    peptide.start=exon.start-(peptide.trans_start-exon.trans_start)
            
            if peptide.trans_end<=exon.trans_end and peptide.trans_end>=exon.trans_start:
                peptide.exon2=exon.number
                if peptide.strand=='+':
                    peptide.end=exon.start+peptide.trans_end-exon.trans_start+1
                else:
                    peptide.end=exon.start-(peptide.trans_end-exon.trans_start)
    if peptide.strand=="+":
        peparray.sort(key=operator.attrgetter('start'))
    else:
        peparray.sort(key=operator.attrgetter('start'),reverse=True)
    
    
    genePSMcount=0
    uniquecluster=[]
    identified_variants=[] # a variant is identified if a peptide uniquely mapped this variant
    pep_number=0
    variant_list=set(gene_variant[ECgene]) # all predicted variants from ECgene database
    #write output for in mappingout.txt
    for peptide in peparray:
        mapped_variants=[] # splice variants contain this peptide sequences
        pep_number+=1
        genePSMcount+=peptide.PSMcount
        peptide.number=pep_number
        uniquecluster.append(peptide.cluster)
        for varid in variant_list:
            if peptide.seq in variant_dic[varid].seq:
                mapped_variants.append(varid)
        if len(variant_list)>1 and len(peptide.variants)==1:
            identified_variants+=peptide.variants
        line="%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n" % (peptide.seq,peptide.number,peptide.length,peptide.HGNC,peptide.proteinID,peptide.PSMcount,peptide.ratio,peptide.error,peptide.cluster,len(variant_list),len(peptide.variants),";".join(mapped_variants),peptide.chr,peptide.start,peptide.end,peptide.strand,peptide.trans_start,peptide.trans_end,peptide.exon1,peptide.exon2)
        
        output_handle.write(line)
    
    #write output for in genestatistics.txt
    cluster=len(set(uniquecluster))
    uniquevar=list(set(identified_variants))
    num_var=len(uniquevar)
    newline="%s\t%d\t%d\t%d\t%d\t%d\t%s\n"%(ECgene,cluster,len(variant_list),genePSMcount,
                                            len(peparray),num_var,",".join(uniquevar))
    output2_handle.write(newline)


##################map peptide sequence#############################
if __name__=='__main__':
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! wrong command, please read the mannual in Readme.txt."
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['ECgene=',
                                                             'ECprotein=',
                                                             'i=','prefix='])
    for opt, arg in options:
        if opt == '--ECgene': gene_structure_file=arg
        elif opt == '--ECprotein': protein_seq_file=arg
        elif opt == '--i':infilename=arg
        elif opt == '--prefix':prefix=arg
        else:
            print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

    input_file=open(infilename,'r')
    input_file.readline()
    gene_peparray={} # store peptide info in _pepdata.txt
    variant_dic={} # store ECvariant info in _pepdata.txt
    gene_variant={} #store ECgene info in _pepdata.txt
    print "read peptide data file"
    for line in input_file:
        row=line.strip().split("\t")
        peptide=PEPTIDE(seq=row[0],HGNC=row[1],proteinID=row[2],PSMcount=int(row[3]),ratio=row[4],error=row[5],cluster=row[6])
        peptide.length=len(peptide.seq)
        peptide.variants=row[2].split(";")
        ECgene=peptide.variants[0].split(".")[0]
        gene_variant[ECgene]=[]
        if ECgene not in gene_peparray:
            gene_peparray[ECgene]=[peptide]
        else:
            gene_peparray[ECgene].append(peptide)

    print 'there are',len(gene_peparray),'genes in total'
    print gene_variant.keys()
    input_file.close()
    print 'reading EC protein sequence file.....'
    
    protein_db=open(protein_seq_file,'r')
    for line in protein_db:
        row=line.strip().split('\t')
        geneid=row[0].split(".")[0]
        if geneid in gene_variant:
            variant=ISOFORM(id=row[0],seq=row[1])
            variant_dic[row[0]]=variant

    protein_db.close()
    print "predicted variants for this gene"
    print variant_dic.keys()
    gene_db=open(gene_structure_file,'r')
    
    print 'reading EC gene structure file.....'
    for line in gene_db:
        row=line.strip().split("\t")
        variantID=row[0]
        geneID=row[0].split(".")[0]
        if geneID in gene_variant:
            gene_variant[geneID].append(row[0])
        if variantID in variant_dic:
            chr_start_list=row[8][:-1].split(",")
            chr_end_list=row[9][:-1].split(",")
            variant_dic[variantID].chr=row[1]
            variant_dic[variantID].strand=row[2]
            variant_dic[variantID].cdsStart=int(row[5])
            variant_dic[variantID].cdsEnd=int(row[6])
            
            
            variant_dic[variantID].exon=[]
            tx=1
            for i in range(len(chr_start_list)):
                exon=EXON(variant=row[0],chr=row[1],strand=row[2])
                exon.number=i+1
                if int(chr_start_list[i])>=variant_dic[variantID].cdsStart and int(chr_end_list[i])<=variant_dic[variantID].cdsEnd:
                    exon.start=int(chr_start_list[i])
                    exon.end=int(chr_end_list[i])
                    exonlen=exon.end-exon.start
                    exon.trans_start=tx
                    exon.trans_end=tx+exonlen-1
                    tx+=exonlen
                elif int(chr_end_list[i])>=variant_dic[variantID].cdsStart and int(chr_start_list[i])<=variant_dic[variantID].cdsStart:
                    exon.start=variant_dic[variantID].cdsStart
                    exon.end=int(chr_end_list[i])
                    exonlen=exon.end-exon.start
                    exon.trans_start=tx
                    exon.trans_end=tx+exonlen-1
                    tx+=exonlen
                elif int(chr_end_list[i])>=variant_dic[variantID].cdsEnd and int(chr_start_list[i])<=variant_dic[variantID].cdsEnd:
                    exon.start=int(chr_start_list[i])
                    exon.end=variant_dic[variantID].cdsEnd
                    exonlen=exon.end-exon.start
                    exon.trans_start=tx
                    exon.trans_end=tx+exonlen-1
                    tx+=exonlen
                
                if row[2]=="+":
                    variant_dic[variantID].exon.append(exon)
                else:
                    st=exon.trans_start
                    ed=exon.trans_end
                    chr_st=exon.start
                    chr_ed=exon.end
                    exon.trans_end=3*len(variant_dic[variantID].seq)-st+1
                    exon.trans_start=3*len(variant_dic[variantID].seq)-ed+1
                    exon.start=chr_ed
                    exon.end=chr_st
                    exon.number=len(chr_start_list)-i
                    variant_dic[variantID].exon.append(exon)


    gene_db.close()
    print 'start mapping...'
    output1=prefix+'_ECmappingout.txt'
    
    output_handle=open(output1,'w')
    headline1=['peptide sequence','peptide numbering','peptide length','gene symbol','protein accession','PSM count','ratio',
               'standard_dev','PQPQ cluster','ECvariants NO.','NOMV','ECvariant ID','chr','chr_start','chr_end','strand',
               'trans_start','trans_end','exon1','exon2']
    output_handle.write('%s\n'%('\t'.join(headline1)))
   
    output2=prefix+'_ECgenestatistics.txt'
    output2_handle=open(output2,'w')
    headline2=['gene symbol','detected clusters NO.','known variants NO.',
              'PSM count','unique peptides','identified variants','variants ID']
    output2_handle.write('%s\n'%('\t'.join(headline2)))
   
   
    for ECgene in gene_peparray.keys():
        peparray=gene_peparray[ECgene]
        main()

    print 'program finished'
    print output1,'saved'
    print "Attention!Genomic coordiantes in ECmappingout.txt file are from NCBI36/hg18 assembly"
    print output2,'saved'