import sys
from Bio import SeqIO
from collections import OrderedDict
import operator


class ISOFORM(object):
    def __init__(self,id=None,chr=None,strand=None):
        self.id=id
        self.strand=strand
        self.chr=chr


class EXON(object):
    def __init__(self,number=0,gene=None,variant=None,chr=None,strand=None,start=0,end=0,length=0,trans_start=0,trans_end=0):
        self.gene=gene
        self.variant=variant
        self.number=number
        self.start=start  #chromosome start coordinate
        self.end=end  #chromosome end coordinate
        self.strand=strand
        self.chr=chr
        self.length=length
        self.trans_start=trans_start
        self.trans_end=trans_end

class PEPTIDE(object):
    def __init__(self,HGNC=None,proteinID=None,seq=None,number=0,cluster=0,chr=None,strand=None,start=0,end=0,length=0,trans_start=0,trans_end=0,exon1=0,exon2=0,PSMcount=0,ratio=None,error=None,variants=[]):
        self.HGNC=HGNC
        self.proteinID=proteinID
        self.seq=seq
        self.number=number
        self.cluster=cluster #clustering result
        self.start=start #chromosome start coordinates
        self.end=end #chromosome end coordinates
        self.strand=strand
        self.chr=chr
        self.length=length
        self.trans_start=trans_start
        self.trans_end=trans_end
        self.PSMcount=PSMcount
        self.ratio=ratio
        self.error=error
        self.variants=variants #a list of variants ID the peptide is mapped to
        self.exon1=exon1 #which exon of the first aa of peptide is mapped to
        self.exon2=exon2 #which exon of the last aa of peptide is mapped to

def main(): #the main function for mapping peptides
    for i in range(len(peparray)):
        peptide=peparray[i]
        peptide.variants=[]
        CDS=""
        klist=[]
        for i in range(0,len(variant_list)):
            pep_index=str(record_dict[variant_list[i]].seq).find(peptide.seq)
            if pep_index!=-1:
                peptide.variants.append(variant_list[i])
                k=i 
                CDR=record_dict[variant_list[i]].description.split('|')[3]
                CDS=CDR[CDR.index("=")+1:] #transcript coding start position
                if CDR[-1]!="=":
                    klist.append(i)
        
        num=len(peptide.variants)
        if len(klist)>1:
            k=klist[0] #k decides which splice variant in the list the peptide will be mapped to.
        if num!=0:# if the peptide map to at least one of splice variants
            n=str(record_dict[variant_list[k]].seq).index(peptide.seq)
            CDR=record_dict[variant_list[k]].description.split('|')[3]
            CDS=CDR[CDR.index("=")+1:] #CDS, trans start 
            if CDS=="":
                CDR=record_dict[variant_list[k]].description.split('|')[1]
                CDS=CDR[CDR.index("=")+1:]
            
            peptide.trans_start=int(CDS)+n*3
            peptide.trans_end=peptide.trans_start+3*len(peptide.seq)-1
            peptide.chr=variant_dic[variant_list[k]].chr
            peptide.strand=variant_dic[variant_list[k]].strand
            
            exons=variant_exon[variant_list[k]]
            for exon in exons:
                if peptide.trans_start<=exon.trans_end and peptide.trans_start>=exon.trans_start:
                    peptide.exon1=exon.number
                    if peptide.strand=='+':
                        peptide.start=exon.start+(peptide.trans_start-exon.trans_start)
                    else:
                        peptide.start=exon.start-(peptide.trans_start-exon.trans_start)
                
                if peptide.trans_end<=exon.trans_end and peptide.trans_end>=exon.trans_start:
                    peptide.exon2=exon.number
                    if peptide.strand=='+':
                        peptide.end=exon.start+peptide.trans_end-exon.trans_start
                    else:
                        peptide.end=exon.start-(peptide.trans_end-exon.trans_start)
    
    if peptide.strand=="+":
        peparray.sort(key=operator.attrgetter('start'))
    else:
        peparray.sort(key=operator.attrgetter('start'),reverse=True)

    pep_number=0 #number the peptide based on its chromosome coordinate
    
    genePSMcount=0
    uniquecluster=[]
    identified_variants=[]
    #write output for in mappingout.txt
    for peptide in peparray:
        pep_number+=1
        genePSMcount+=peptide.PSMcount
        peptide.number=pep_number
        uniquecluster.append(peptide.cluster)
        
        if len(variant_list)>1 and len(peptide.variants)==1:
            identified_variants+=peptide.variants
        line="%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n" % (peptide.seq,peptide.number,peptide.length,peptide.HGNC,peptide.proteinID,peptide.PSMcount,peptide.ratio,peptide.error,peptide.cluster,len(variant_list),len(peptide.variants),",".join(peptide.variants),peptide.chr,peptide.start,peptide.end,peptide.strand,peptide.trans_start,peptide.trans_end,peptide.exon1,peptide.exon2)

        output_handle.write(line)
    
    #write output for in genestatistics.txt
    cluster=len(set(uniquecluster))
    uniquevar=list(set(identified_variants))
    num_var=len(uniquevar)
    newline="%s\t%d\t%d\t%d\t%d\t%d\t%s\n"%(gene,cluster,len(variant_list),genePSMcount,
                                            len(peparray),num_var,",".join(uniquevar))
    output2_handle.write(newline)




##################map peptide sequence#############################
if __name__=='__main__':
    input_file=open(sys.argv[1],'r') #_pepdata.txt or pepcluster.txt
    input_file.readline()
    gene_peparray={}
    for line in input_file:
        row=line.strip().split("\t")
        gene=row[1]
        peptide=PEPTIDE(seq=row[0],HGNC=gene,proteinID=row[2],PSMcount=int(row[3]),ratio=row[4],error=row[5],cluster=row[6])
        peptide.length=len(peptide.seq)
        if gene not in gene_peparray:
            gene_peparray[gene]=[peptide]
        else:
            gene_peparray[gene].append(peptide)
    
    print 'there are',len(gene_peparray),'genes in total'
    input_file.close()
    print 'mapping.....'
    
    record_dict=SeqIO.index('varseq.fa','fasta')
    
    handle=open('splicingvar.txt')
    gene_variant={}
    variant_dic={}
    variant_exon={}
    for line in handle:
        row=line.strip().split("\t")
        exon=EXON(gene=row[1],chr=row[2],strand=row[3],variant=row[4],number=row[5],
                  start=int(row[6]),end=int(row[7]),trans_start=int(row[8]),trans_end=int(row[9]))
        
        gene=row[0]
        variantID=row[4]
        if variantID not in variant_exon:
            variant=ISOFORM(id=variantID,chr=row[2],strand=row[3])
            variant_dic[variantID]=variant
            variant_exon[variantID]=[exon]
        else:
            variant_exon[variantID].append(exon)
        if gene not in gene_variant:
            gene_variant[gene]=[variantID]
        else:
            gene_variant[gene].append(variantID)
    
    
    handle.close()

    handle3=open('subexon.txt')
    gene_subexon={}
    for line in handle3:
        row=line.strip().split("\t")
        subexon=EXON(gene=row[1],chr=row[2],strand=row[3],number=int(row[4]),start=int(row[6]),end=int(row[7]))
        gene=row[0]
        if gene not in gene_subexon:
            gene_subexon[gene]=[subexon]
        else:
            gene_subexon[gene].append(subexon)

    handle3.close()
    
    prefix=sys.argv[2]
    output1=prefix+'_mappingout.txt'
    output2=prefix+'_genestatistics.txt'
    
    output_handle=open(output1,'w')
    headline1=['peptide sequence','peptide numbering',
               'peptide length','gene symbol','protein accession','PSM count','ratio',
               'variation between PSMs','PQPQ cluster','known variants NO.','NOMV',
               'known variants peptide is mapped to','chr','chr_start','chr_end','strand',
               'trans_start','trans_end','exon1','exon2']
    output_handle.write('%s\n'%('\t'.join(headline1)))
    
    output2_handle=open(output2,'w')
    headline2=['gene symbol','detected clusters NO.','known variants NO.',
               'PSM count','unique peptides','identified variants','variants ID']
    output2_handle.write('%s\n'%('\t'.join(headline2)))
    
    icount=0
    for gene in gene_peparray.keys():
        icount+=1
        if gene in gene_variant:#if the gene's splice info downloaded
            peparray=gene_peparray[gene]
            variant_list=list(set(gene_variant[gene]))
            main()
        if icount%1000==0:
            print icount,"genes processed"
    
    
    print 'program finished'
    print output1,'saved'
    print output2,'saved'