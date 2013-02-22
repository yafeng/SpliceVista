import sys
import numpy

def getpep(array): #this function is to find unique peptides in the psmarray
    dic={}
    for i in range(0,len(array)):
        pep=array[i][2].upper()
        if pep not in dic:
            dic[pep]=1
    return dic.keys()


def main():
    peplist=getpep(genearray)
    for i in range(0,len(peplist)):
        pep=peplist[i]
        pepratio=[]
        for j in range(0,len(genearray)):
            if genearray[j][2].upper()==pep:
                pepratio.append(genearray[j][3:])

        a=numpy.array(pepratio,dtype=float)
        mean=numpy.mean(a,axis=0)
        mean_round=[ '%.2f' % elem for elem in mean ]
        stdv=numpy.std(a,axis=0)
        stdv_round=[ '%.2f' % elem for elem in stdv ]

        newline="%s\t%s\t%i\t%s\t%s\t%s\n" % (gene,pep,len(a),','.join(mean_round),
                                      ','.join(stdv_round),'0')
        output.write(newline)
        
if __name__=='__main__':
    prefix=sys.argv[1]
    infilename=prefix+'_psmdata.txt'
    handle=open(infilename,'r')
    handle.readline()
    newheader=['gene','pep','PSM count','ratio','standard_dev','cluster']
    firstline='\t'.join(newheader)+'\n'
    
    outfilename=prefix+'_pepdata.txt'
    output=open(outfilename,'w')
    output.write(firstline)
    
    gene_psmarray={}
    for line in handle:
        row=line[:-1].split("\t")
        gene=row[1]
        if ";" in gene: # this will discard PSMs with ambiguous identity
            continue;
        else:
            if gene not in gene_psmarray:
                gene_psmarray[gene]=[row]
            else:
                gene_psmarray[gene].append(row)
    
    print 'there are',len(gene_psmarray),'genes in total'
    handle.close()
    
    i=0
    for gene in gene_psmarray.keys():
        i+=1
        genearray=gene_psmarray[gene]
        main()
        if i%1000==0:
            print i,'genes processed'
    print len(gene_psmarray),'genes processed'
    print 'program finished'    
