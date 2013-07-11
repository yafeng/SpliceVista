import sys
from Bio import SeqIO
from collections import OrderedDict
import operator

########################Splicing Variants Structure#########
def genevariant(lis): #returns a list of variants the gene has.
    varlist=OrderedDict()
    for i in range(0,len(lis)):
        if lis[i][4] not in varlist:
            varlist[lis[i][4]]=1
    return varlist.keys();

def var_cor(lis1,lis2):   #returns an array which stores exons chr_cor from one var together 
    cordinate=[]
    for i in range(0,len(lis1)):
        cordinate.append([])
        for j in range(0,len(lis2)):
            if lis1[i]==lis2[j][4]:
                cordinate[i].append((int(lis2[j][6]),int(lis2[j][7])));

    return cordinate;

def Exon_cor(lis1,lis2):   #returns an array which stores exons aa_cor from one var together 
    cordinate=[]
    for i in range(0,len(lis1)):
        cordinate.append([])
        for j in range(0,len(lis2)):
            if lis1[i]==lis2[j][4]:
                cordinate[i].append((int(lis2[j][8]),int(lis2[j][9])));

    return cordinate;

def Exon_no(lis1,lis2):  #returns an list which stores exon position of variants
    cordinate=[]
    for i in range(0,len(lis1)):
        cordinate.append([])
        for j in range(0,len(lis2)):
            if lis1[i]==lis2[j][4]:
                cordinate[i].append(int(lis2[j][5]));

    return cordinate;

def Countexon(lis):
    exonlis=[]
    for i in range(0,len(lis)):
        exonlis.append(int(lis[i][5]))

    return max(exonlis)
        

###################Maping variants from PQPQ######################

def pepmap(x,lis): #this function search tuple x in a tuple list then return its position.
    index=0 
    for i in range(0,len(lis)):
        if x>=min(lis[i]) and x<=max(lis[i]):
            index=i
            break;

    return index


def getrecord(lis,records):
    recordlist=[]
    for record in records:
        if record.id in lis:
            recordlist.append(record)

    return recordlist

def sort_table(table, col=0):
    return sorted(table, key=operator.itemgetter(col))
###################Function part end###############################
def main():
    pep_cluster=[]
    NOMV=[]
    MTV=[]
    Exon_st=[]
    Exon_ed=[]
    Chr_start=[]
    Chr_end=[]
    pep_chr=[]
    
    var=genevariant(vararray)
    Exon_pos1=Exon_no(var,vararray) 
    varcor=var_cor(var,vararray)
    CDR_cor=Exon_cor(var,vararray)
    
    for j in range(0,len(peparray)):
        pep_cluster.append(int(peparray[j][-1]))
        seq=peparray[j][0]
        num=0 #count how many known variants the peptide mapped to
        MTV.append([])
        CDS=""
        klist=[]
        for i in range(0,len(var)):
            m=str(record_dict[var[i]].seq).find(seq)
            if m!=-1:
                num=num+1
                MTV[j].append(var[i])
                k=i 
                CDR=record_dict[var[i]].description.split('|')[3]
                CDS=CDR[CDR.index("=")+1:]
                if CDR[-1]!="=":
                    klist.append(i)
        
        NOMV.append(num)
        if len(klist)>1:
            k=klist[0]
        if num!=0:# if the peptide map to one of splicing variants
            n=str(record_dict[var[k]].seq).index(seq)
            CDR=record_dict[var[k]].description.split('|')[3]
            CDS=CDR[CDR.index("=")+1:] #CDS, trans start 
            if CDS=="":
                CDR=record_dict[var[k]].description.split('|')[1]
                CDS=CDR[CDR.index("=")+1:]
            exon_st=int(CDS)+n*3
            exon_ed=exon_st+3*len(seq)-1
            
            
            index1=pepmap(exon_st,CDR_cor[k])
            index2=pepmap(exon_ed,CDR_cor[k])
            exon_order1=Exon_pos1[k][index1]
            exon_order2=Exon_pos1[k][index2]
            #print seq,CDS,n,exon_st,exon_ed,exon_order1,exon_order2
            Exon_st.append(exon_order1)
            Exon_ed.append(exon_order2)

            if vararray[0][3]=='+':                       
                chr_st=min(varcor[k][index1])+exon_st-min(CDR_cor[k][index1])
                Chr_start.append(chr_st)
                chr_ed=min(varcor[k][index2])+exon_ed-min(CDR_cor[k][index2])
                Chr_end.append(chr_ed)
            elif vararray[0][3]=='-':
                chr_st=max(varcor[k][index1])-(exon_st-min(CDR_cor[k][index1]))
                Chr_start.append(chr_st)
                chr_ed=max(varcor[k][index2])-(exon_ed-min(CDR_cor[k][index2]))
                Chr_end.append(chr_ed)
        else:#if peptide doesn't map to any of splicing variants
            Exon_st.append('')
            Exon_ed.append('')
            Chr_start.append('')
            Chr_end.append('')
            chr_st=0

        pep_chr.append([seq,chr_st])         

    pep_chr_sorted=sort_table(pep_chr,1) #sort the peptides according to genomic coordinates

    uniqpep={}
    n=0
    for i in range(0,len(pep_chr_sorted)):
        pep=pep_chr_sorted[i][0]
        if pep not in uniqpep:
            n+=1
            uniqpep[pep]=n
    
    uniq_cluster=list(set(pep_cluster))
    if 0 in uniq_cluster:
        uniq_cluster.remove(0) # 0 means unclustered peptides
    
        
    var_identified={} #store the ID of identified splice variants for this gene
    gene_psmcount=0
    for i in range(0,len(peparray)):
        pep=peparray[i][0]
        gene_psmcount+=int(peparray[i][3])
        pepratio=peparray[i][4]
        s1=str(len(pep))
        s2=str(len(uniq_cluster))
        s3=str(len(var))
        s4=str(NOMV[i])
        s5=','.join(MTV[i])
        s6=str(Chr_start[i])
        s7=str(Chr_end[i])
        s8=str(Exon_st[i])
        s9=str(Exon_ed[i])
        if len(var)>1 and NOMV[i]==1:
            if s5 not in var_identified:
                var_identified[s5]=1

        peplabel=str(uniqpep[pep])
        line="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (pep,peplabel,s1,"\t".join(peparray[i][1:]),
                s2,s3,s4,s5,s6,s7,s8,s9)
        
        output_handle.write(line)

    n1=len(uniq_cluster)
    n2=len(var)    
    n3=gene_psmcount
    n4=len(uniqpep)
    n5=len(var_identified)
    s10=','.join(var_identified.keys())
    exonlis=Exon_st+Exon_ed
    uniqexon=list(set(exonlis))
    exoncoverage=float(len(uniqexon))/Countexon(vararray)
    newline="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(gene,str(n1),str(n2),
                                                        str(n3),str(n4),
                                                        str(exoncoverage),
                                                        str(n5),s10)
    
    output2_handle.write(newline)

    
##################map peptide sequence#############################
if __name__=='__main__':
    input_file=open(sys.argv[1],'r')
    input_file.readline()
    gene_peparray={}
    for line in input_file:
        row=line[:-1].split("\t")
        gene=row[1]
        if gene not in gene_peparray:
            gene_peparray[gene]=[row]
        else:
            gene_peparray[gene].append(row)
    
    print 'there are',len(gene_peparray),'genes in total'
    input_file.close()
    print 'mapping.....'

    record_dict=SeqIO.index('varseq.fa','fasta')
    
    handle=open('splicingvar.txt')
    gene_vararray={}
    for line in handle:
        row=line[:-1].split("\t")
        gene=row[1]
        if gene not in gene_vararray:
            gene_vararray[gene]=[row]
        else:
            gene_vararray[gene].append(row)
    
    handle.close()
    
    prefix=sys.argv[2]
    output1=prefix+'_mappingout.txt'
    output2=prefix+'_genestatistics.txt'
    
    output_handle=open(output1,'w')
    headline1=['peptide sequence','peptide numbering',
               'peptide length','gene symbol','protein accession','PSM count','ratio',
               'standard_dev','PQPQ cluster',
               'detected clusters NO.','known variants NO.','NOMV',
               'known variants peptide is mapped to','chr_start','chr_end',
               'exon_start','exon_end']
    output_handle.write('%s\n'%('\t'.join(headline1)))
    output2_handle=open(output2,'w')
    headline2=['gene symbol','detected clusters NO.','known variants NO.',
               'PSM count','unique peptides','exon coverage','identified variants','variants ID']
    output2_handle.write('%s\n'%('\t'.join(headline2)))

    icount=0
    for gene in gene_peparray.keys():
        icount+=1
        if gene in gene_vararray:#if the gene's splice info downloaded
            peparray=gene_peparray[gene]
            vararray=gene_vararray[gene]
            main()
        if icount%1000==0:
            print icount,"genes processed"
            
                
    print 'program finished'
    print output1,'saved'
    print output2,'saved'
