import sys
import numpy

def getgene(infile): #get gene list out of infile 
    dic={}
    for line in infile:
        gene=line.split('\t')[1]
        if gene not in dic:
            dic[gene]=1

    return dic.keys();

def extractpsm(gene,infile):#extract all psm for one gene
    array=[]
    for line in infile:
        if line.split('\t')[1]==gene:
            array.append(line[:-1].split('\t')); 
    return array

def getpep(array):
    dic={}
    for i in range(0,len(array)):
        pep=array[i][2].upper()
        if pep not in dic:
            dic[pep]=1
    return dic.keys()

def normalize(lis):
    normalized=[]
    for i in range(0,len(lis)):
        normvalue=2*float(lis[i])/(float(lis[0])+float(lis[1]))
        normalized.append(normvalue)

    return normalized

def main():
    infile=open(sys.argv[1],'r')
    genearray=extractpsm(gene,infile)
    infile.close()
    peplist=getpep(genearray)
    for i in range(0,len(peplist)):
        pep=peplist[i]
        foldchange=[]
        for j in range(0,len(genearray)):
            if genearray[j][2].upper()==pep:
                pepsignal=genearray[j][3:11]
                pepratio=normalize(pepsignal)
                foldchange.append(pepratio)

        a=numpy.array(foldchange,dtype=float)
        mean=numpy.mean(a,axis=0)
        mean_f=[ '%.2f' % elem for elem in mean ]
        stdv=numpy.std(a,axis=0)
        stdv_f=[ '%.2f' % elem for elem in stdv ]

        newline="%s\t%s\t%i\t%s\t%s\t%s\n" % (gene,pep,len(a),','.join(mean_f),
                                      ','.join(stdv_f),'0')
        output.write(newline)
        
if __name__=='__main__':
    handle=open(sys.argv[1],'r')# take sample_pepdata.txt as input file
    handle.readline()
    newheader=['gene','pep','PSM count','foldchange','standard_dev','cluster']
    firstline='\t'.join(newheader)+'\n'
    genelist=getgene(handle)
    print 'there are',len(genelist),'genes in total'
    handle.close()
    
    output=open(sys.argv[2],'w')
    output.write(firstline)
    for i in range(0,len(genelist)):
        gene=genelist[i]
        main()
        if i%100==0:
            print i,'genes processed'
    print len(genelist),'genes processed'
    print 'program finished'    
