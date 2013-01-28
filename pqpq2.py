import sys,os
import getopt
import numpy
import scipy.cluster.hierarchy as sch

def getgene(infile): #get gene list out of infile 
    dic={}
    for line in infile:
        gene=line.split('\t')[0]
        if gene not in dic:
            dic[gene]=1

    return dic.keys();
def extractpep(gene,infile):#extract all pep for one gene
    array=[]
    for line in infile:
        row=line[:-1].split('\t')
        if row[0]==gene:
            array.append(row[:-1]); #remove last column from pepdata.txt
    return array

def main(): #clustering and write output
    infile=open(infilename,'r')
    pep_array=extractpep(gene,infile)
    infile.close()
    if len(pep_array)>1:
        matrix=[]
        for i in range(0,len(pep_array)):
            matrix.append(pep_array[i][3].split(','))

        dataMatrix=numpy.array(matrix,dtype=float)
        d = sch.distance.pdist(dataMatrix,metric)# vector of pairwise distances
        D = numpy.clip(d,0,2)
        L = sch.linkage(D, method,metric)
        ind = sch.fcluster(L,distance,'distance')#distance is dissmilarity(1-correlation)
        p=numpy.array(pep_array)
        p=numpy.column_stack([p,ind])
        formatoutput(p)
    else:
        p=numpy.array(pep_array)
        p=numpy.column_stack([p,[0]])
        formatoutput(p)
        
def formatoutput(array): #format output
    for i in range(0,len(array)):
        newline='\t'.join(array[i]) +'\n'
        output.write(newline)
    
    
if __name__=='__main__':
    ################  Default  ################
    method = 'average'
    distance= 0.2
    metric = 'euclidean'
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! wrong command, please read the mannual in Readme.txt."
        print "Example: python pqpq2.py --i heavy_pepdata.txt --o heavy_pqpqout.txt"
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['metric=',
                                                             'method=','d=',
                                                             'i=','o='])
        print options,remainder
        for opt, arg in options:
            if opt == '--metric': metric=arg
            elif opt == '--method': method=arg
            elif opt == '--d': distance=arg
            elif opt == '--i':infilename=arg
            elif opt == '--o':outfilename=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    
    handle=open(infilename,'r')
    handle.readline()
    newheader=['gene','pep','PSM count','foldchange','standard_dev','cluster']
    firstline='\t'.join(newheader)+'\n'
    genelist=getgene(handle)
    print metric,"metric is used to measure the distance of elements"
    print "distance cutoff to form a new cluster is",distance
    print 'there are',len(genelist),'genes in total'
    handle.close()
    
    output=open(outfilename,'w')
    output.write(firstline)
    for i in range(0,len(genelist)):
        gene=genelist[i]
        main()
        if i%1000==0:
            print i,'genes processed'
    print len(genelist),'genes processed'
    print 'program finished'    
