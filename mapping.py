import sys
from Bio import SeqIO
from collections import OrderedDict
import operator

def getgene(infile): #get gene symbol list out of pqpq file 
    lis=OrderedDict()
    infile.readline()
    for line in infile:
        gene=line.split('\t')[0]
        if gene not in lis:
            lis[gene]=1

    return lis.keys();

def extractvar(gene,infile):#extract splicevariant info for one gene
    array=[]
    for line in infile:
        if line.split('\t')[0]==gene:
            array.append(line[:-2].split('\t')); #the last field in splicingvar file is '0\r\n'

    return array

def extractpep(gene,infile):#extract all pep for one gene
    array=[]
    for line in infile:
        if line.split('\t')[0]==gene:
            array.append(line[:-1].split('\t')); 
    return array  

#######################Subexons structure############


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

def normalize(lis):
    normalized=[]
    for i in range(0,len(lis)):
        normvalue=2*float(lis[i])/(float(lis[0])+float(lis[1]))
        normalized.append(normvalue)

    return normalized

def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if length % 2==0:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return sorts[length / 2]

def sort_table(table, col=0):
    return sorted(table, key=operator.itemgetter(col))
###################Function part end###############################
def main():
    pqpq_cluster=[]
    NOMV=[]
    MTV=[]
    Exon_st=[]
    Exon_ed=[]
    Chr_start=[]
    Chr_end=[]
    pep_chr=[]
    plex7=[]
    plex8=[]
    
    handle1=open('splicingvar.txt')                           
    array1=extractvar(gene,handle1)# the array stores splicing variants file the gene has
    var=genevariant(array1)
    Exon_pos1=Exon_no(var,array1) 
    varcor=var_cor(var,array1)
    CDR_cor=Exon_cor(var,array1)
    
    handle1.close()

    handle2=open(infilename)
    array2=extractpep(gene,handle2)
    handle2.close()

    for j in range(0,len(array2)):
        pqpq_cluster.append(int(array2[j][5]))
        seq=array2[j][1]
        num=0 #count how many known variants the peptide mapped to
        MTV.append([])
        for i in range(0,len(var)):
            m=str(record_dict[var[i]].seq).find(seq)
            if m!=-1:
                num=num+1
                MTV[j].append(var[i])
                k=i
                n=m;

        #print var[i],k,seq,n
        NOMV.append(num)
        if num!=0:# if the peptide map to one of splicing variants
            CDR=record_dict[var[k]].description.split('\t')[1]
            if CDR.find('<')!=-1: #usually in 
                CDS=int(CDR[CDR.index('<')+1:CDR.index(':')])
            else:
                CDS=int(CDR[CDR.index('[')+1:CDR.index(':')])

            exon_st=CDS+(n-1)*3
            exon_ed=exon_st+3*len(seq)
            
            
            index1=pepmap(exon_st,CDR_cor[k])
            index2=pepmap(exon_ed,CDR_cor[k])
            exon_order1=Exon_pos1[k][index1]
            exon_order2=Exon_pos1[k][index2]
            #print seq,CDS,n,exon_st,exon_ed,exon_order1,exon_order2
            Exon_st.append(exon_order1)
            Exon_ed.append(exon_order2)
            

            if array1[0][3]=='+':                       
                chr_st=min(varcor[k][index1])+exon_st-min(CDR_cor[k][index1])
                Chr_start.append(chr_st)
                chr_ed=min(varcor[k][index2])+exon_ed-min(CDR_cor[k][index2])
                Chr_end.append(chr_ed)
            elif array1[0][3]=='-':
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

    pep_chr_sorted=sort_table(pep_chr,1)

    uniqpep={}
    n=0
    for i in range(0,len(pep_chr_sorted)):
        pep=pep_chr_sorted[i][0]
        if pep not in uniqpep:
            n+=1
            uniqpep[pep]=n
            
    uniq_cluster=list(set(pqpq_cluster))
    if 0 in uniq_cluster:
        uniq_cluster.remove(0) # 0 means unclustered peptides
    
        
    var_identified={} #store the ID of identified splice variants for this gene
    psmsum=0
    for i in range(0,len(array2)):
        pep=array2[i][1]
        psmsum+=int(array2[i][2])
        pepratio=array2[i][3]
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
        line="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene,pep,
                peplabel,s1,array2[i][2],pepratio,array2[i][5],array2[i][4],
                s2,s3,s4,s5,s6,s7,s8,s9)

       
        peplabel=str(uniqpep[pep])

        
        output_handle.write(line)

    n1=len(uniq_cluster)
    n2=len(var)    
    n3=psmsum
    n4=len(uniqpep)
    n5=len(var_identified)
    s10=','.join(var_identified.keys())
    exonlis=Exon_st+Exon_ed
    uniqexon=list(set(exonlis))
    exoncoverage=float(len(uniqexon))/Countexon(array1)
    newline="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(gene,str(n1),str(n2),
                                                        str(n3),str(n4),
                                                        str(exoncoverage),
                                                        str(n5),s10)
    
    output2_handle.write(newline)

    
##################map peptide sequence#############################
if __name__=='__main__':
    prefix=sys.argv[1]
    infilename=prefix+'_cluster.txt'
    input_file=open(infilename,'r')
    genelist=getgene(input_file)
    print 'there are',len(genelist),'genes in the file'
    input_file.close()

    record_dict=SeqIO.index('varseq.fa','fasta')

    output1=prefix+'_mappingout.txt'
    output2=prefix+'_genestatistics.txt'
    
    output_handle=open(output1,'w')
    headline1=['gene symbol','peptide sequence','peptide label',
               'peptide length','PSM count','foldchange','PQPQ cluster',
               'standard_dev',
               'detected clusters NO.','known variants NO.','NOMV',
               'known variants peptide is mapped to','chr_start','chr_end',
               'exon_start','exon_end']
    output_handle.write('%s\n'%('\t'.join(headline1)))
    output2_handle=open(output2,'w')
    headline2=['gene symbol','detected clusters NO.','known variants NO.',
               'PSM count','unique peptides','exon coverage','identified variants','variants ID']
    output2_handle.write('%s\n'%('\t'.join(headline2)))

    file3=open('gene_notfound.txt','r')
    notfounddict={}
    
    for line in file3:
        if line[:-1] not in notfounddict:
            notfounddict[line[:-1]]=1 #get a genelist that are not found in EVDB

    print 'mapping.....'  
    for i in range(0,20):
        gene=genelist[i]
        if gene in notfounddict:
            continue;
        else:
            main()
            print gene,i
                
    print 'program finished'
    print output1,'saved'
    print output2,'saved'
