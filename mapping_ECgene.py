import sys
import os
import getopt
from Bio import SeqIO
from collections import OrderedDict
import operator

########################Splicing Variants Structure#########
def Transcript(start,end,list1,list2):#from a list of eoxn chr start coordiantes and end coordinates generates two list of tuples for each exon [(chr_start,chr_end)],[(Tran_start,Tran_end)]
    exon_chr=[]
    exon_tx=[]
    tx=1
    for i in range(len(list1)):
        if int(list2[i])>start and int(list1[i])<start:
            exon_chr.append((start,int(list2[i])))
            exonlen=int(list2[i])-start
            exon_tx.append((tx,exonlen+tx-1))
            tx+=exonlen
        elif int(list1[i])>start and int(list2[i])<end:
            exon_chr.append((int(list1[i]),int(list2[i])))
            exonlen=int(list2[i])-int(list1[i])
            exon_tx.append((tx,tx+exonlen-1))
            tx+=exonlen
        elif int(list2[i])>end and int(list1[i])<end:
            exon_chr.append((int(list1[i]),end))
            exonlen=end-int(list1[i])
            exon_tx.append((tx,exonlen+tx-1))
    return exon_chr,exon_tx

###################Maping variants from PQPQ######################

def pepmap(x,lis): #this function search tuple x in a tuple list then return its position.
    index=0 
    for i in range(0,len(lis)):
        if x>=min(lis[i]) and x<=max(lis[i]):
            index=i
            break;
    
    return index

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

    var=gene_var[gene]

    for j in range(0,len(peparray)):
        pep_cluster.append(int(peparray[j][-1]))
        seq=peparray[j][0]
        protein_acc=peparray[j][2]
        chr=var_dic[protein_acc][1]
        strand=var_dic[protein_acc][2]
        num=0 #count how many known variants the peptide mapped to
        MTV.append([])
        CDS=""
        klist=[]
        for i in range(0,len(var)):
            m=protein_dict[var[i]].find(seq)
            if m!=-1:
                num=num+1
                MTV[j].append(var[i])
                k=i 
        
        NOMV.append(num)
        if num!=0:# if the peptide map to one of splice variants
            peptide_index=protein_dict[protein_acc].index(seq)
            var_len=3*len(protein_dict[protein_acc])
            pep_st=1+peptide_index*3 #start position of peptide at the transcript
            pep_ed=pep_st+3*len(seq)-1

            coding_start=int(var_dic[protein_acc][5])
            coding_end=int(var_dic[protein_acc][6])
            exon_starts=var_dic[protein_acc][8][:-1].split(",")
            exon_ends=var_dic[protein_acc][9][:-1].split(",")
            exon_chr,exon_tx=Transcript(coding_start,coding_end,exon_starts,exon_ends)
            if seq=="GLSDFLGVISDTFAPSPDKTIDCDVITLMGTPSGTAEPYDGTK":
                print protein_dict[protein_acc],pep_st,pep_ed,coding_start,coding_end,exon_chr,exon_tx
            
            if strand=='+':
                index1=pepmap(pep_st,exon_tx)
                index2=pepmap(pep_ed,exon_tx)
                Exon_st.append(index1+1)
                Exon_ed.append(index2+1)
                chr_st=min(exon_chr[index1])+pep_st-exon_tx[index1][0]+1
                Chr_start.append(chr_st)
                chr_ed=min(exon_chr[index2])+pep_ed-exon_tx[index2][0]+1
                Chr_end.append(chr_ed)
            elif strand=='-':
                index1=pepmap(var_len-pep_st+1,exon_tx)
                index2=pepmap(var_len-pep_ed+1,exon_tx)
                Exon_st.append(len(exon_tx)-index1)
                Exon_ed.append(len(exon_tx)-index2)
                chr_st=max(exon_chr[index2])-(exon_tx[index2][1]-(var_len-pep_ed+1))
                Chr_start.append(chr_st)
                chr_ed=max(exon_chr[index1])-(exon_tx[index1][1]-(var_len-pep_st+1))
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
        line="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (pep,peplabel,s1,"\t".join(peparray[i][1:]),s2,s3,s4,s5,chr,s6,s7,strand,s8,s9)
        
        output_handle.write(line)
    
    n1=len(uniq_cluster)
    n2=len(var)    
    n3=gene_psmcount
    n4=len(uniqpep)
    n5=len(var_identified)
    s10=','.join(var_identified.keys())
    exonlis=Exon_st+Exon_ed
    uniqexon=list(set(exonlis))
    exoncoverage=float(len(uniqexon))/len(exon_tx)
    newline="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(gene,str(n1),str(n2),
                                                str(n3),str(n4),
                                                str(exoncoverage),
                                                str(n5),s10)
    
    output2_handle.write(newline)


##################map peptide sequence#############################
if __name__=='__main__':
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! wrong command, please read the mannual in Readme.txt."
        print "Example: python converter.py --i heavy_msout.txt --prefix heavy --database uniprot"
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
    for line in input_file:
        row=line[:-1].split("\t")
        gene=row[2].split(".")[0] #use ECgene cluster as gene id
        if gene not in gene_peparray:
            gene_peparray[gene]=[row]
        else:
            gene_peparray[gene].append(row)
    
    print 'there are',len(gene_peparray),'genes in total'
    input_file.close()
    print 'mapping.....'
    
    protein_dict={}
    protein_db=open(protein_seq_file,'r')
    for line in protein_db:
        row=line.strip().split('\t')
        protein_dict[row[0]]=row[1]
    
    protein_db.close()
    
    gene_db=open(gene_structure_file,'r')
    var_dic={} #store splice variant chromosome position
    gene_var={} #store var ids for each gene 
    for line in gene_db:
        row=line[:-1].split("\t")
        var=row[0]
        gene=var.split(".")[0]
        var_dic[var]=row
        if gene not in gene_var:
            gene_var[gene]=[var]
        else:
            gene_var[gene].append(var)
        
    
    output1=prefix+'_mappingout.txt'
    output2=prefix+'_genestatistics.txt'
    
    output_handle=open(output1,'w')
    headline1=['peptide sequence','peptide numbering',
               'peptide length','gene symbol','protein accession','PSM count','peptide ratio',
               'standard_dev','cluster',
               'detected clusters NO.','variants NO.','NOMV',
               'splice variants peptide is mapped to','chromosome','chr_start','chr_end','strand',
               'exon_start','exon_end']
    output_handle.write('%s\n'%('\t'.join(headline1)))
    output2_handle=open(output2,'w')
    headline2=['gene symbol','detected clusters NO.','variants NO.',
               'PSM count','unique peptides','exon coverage','identified variants','variants ID']
    output2_handle.write('%s\n'%('\t'.join(headline2)))
    
    icount=0
    for gene in gene_peparray.keys():
        icount+=1
        if gene in gene_var:#if the gene's splice info downloaded
            peparray=gene_peparray[gene]
            main()
        if icount%1000==0:
            print icount,"genes processed"
    
    
    print 'program finished'
    print output1,'saved'
    print output2,'saved'