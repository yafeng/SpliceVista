import sys
import os
import getopt
import numpy as np

def normalize1(lis):
    newlis=[]
    for i in range(0,len(lis)):
        normvalue=(float(lis[i])/factor[i])/(float(lis[0])/factor[0])
        newlis.append(normvalue)
    
    return newlis

def normalize2(lis):
    newlis=[]
    for i in range(0,len(lis)):
        normvalue=2*(float(lis[i])/factor[i])/(float(lis[0])/factor[0]+float(lis[1])/factor[1])
        newlis.append(normvalue)
    
    return newlis

def main():
    N=0
    for protein_acc in acclist.keys():
        proteinList=protein_acc.split(";")
        genename=[]
        for i in range(0,len(proteinList)):
            proteinID=proteinList[i]
            if database=="ECgene":
                ECgene=proteinID.split(".")[0]
                genename.append(db.get(ECgene,ECgene))
            elif database=="ensembl+splice":
                if proteinID[:3]=="SJP":
                    ensg=proteinID.split("_")[1]
                    genename.append(db.get(ensg,"None"))
                else:
                    genename.append(db.get(proteinID,"None"))
            else:
                genename.append(db.get(proteinID,"None"))
        
        if len(set(genename))==1: #if the proIDs correspond to one gene
            acclist[protein_acc]=genename[0]   #return one gene name
        else: #if not, return gene name for each proID
            acclist[protein_acc]=";".join(genename)
        
        N=N+1
        if N%1000==0:
            print N,'protein groups processed'

if __name__=='__main__':
    ################  Default  ################
    database = 'ensembl'
    n='0' #normalization option
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! wrong command, please read the mannual in Readme.txt."
        print "Example: python converter.py --i heavy_msout.txt --prefix heavy --database uniprot"
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['database=',
                                                             'n=',
                                                             'i=','prefix='])
        for opt, arg in options:
            if opt == '--database': database=arg
            elif opt == '--n': n=arg
            elif opt == '--i':infilename=arg
            elif opt == '--prefix':prefix=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    
    if database=="ensembl":
        handle=open("ENSP_Gene.txt",'r')
    elif database=="uniprot":
        handle=open("Uniprotid_Gene.txt",'r')
    elif database=="IPI":
        handle=open("IPIaccesion_gene.txt",'r')
    elif database=="ECgene":
        handle=open("ECgene2HUGO.txt",'r')
    elif database=="ensembl+splice":
        handle=open("ENSP_ENSG_Gene.txt",'r')
    else:
        print database,"is not supported database"
        sys.exit()
    
    db={}
    for line in handle:
        row=line.strip().split('\t')
        db[row[0]]=row[1]
    ############read the input file##############
    input_file=open(infilename,'r')
    acclist={}
    totalPSM=0
    for line in input_file:
        totalPSM+=1
        protein_acc=line.split('\t')[1]
        if protein_acc not in acclist:
            acclist[protein_acc]=1
    
    print totalPSM,'PSMs in this  file'
    input_file.close()
    
    main() #protein_acc in acclist now has the gene symbol as value
    
    ########write the output file################
    print "writing output"
    infile2=open(infilename,'r')
    outfilename=prefix+'_psmdata.txt'
    outfile=open(outfilename,'w')
    
    line=infile2.readline()
    cols=line.split('\t')
    samplesize=len(cols[2:])
    cols.insert(1,'gene symbol')
    outfile.write("\t".join(cols))
    
    intensity_array=[]
    for line in infile2:
        cols=line.strip().split('\t')
        if len(cols[2:])==samplesize:
            try:
                newcols=map(float,cols[2:])
                intensity_array.append(newcols)
            except ValueError:
                continue;
    
    a=np.array(intensity_array,dtype=float)
    median=np.median(a,axis=0)
    factor=median/min(median)
    print median
    print factor
    infile2.close()

    infile2=open(infilename,'r')
    infile2.readline()
    if n=="0":
        for line in infile2:
            cols=line.split('\t')
            proteinID=cols[1]
            cols.insert(1,acclist[proteinID])
            outfile.write("\t".join(cols))
    
    elif n=="1":
        for line in infile2:
            cols=line.strip().split('\t')
            proteinID=cols[1]
            intensity=cols[2:]
            if len(intensity)==samplesize:
                try:
                    ratio=normalize1(intensity)
                    ratio_round=[ '%.3f' % elem for elem in ratio ] #round the value, keep 3 digit after decimal
                
                    cols.insert(1,acclist[proteinID])
                    newcols=cols[:3]+ratio_round
                    outfile.write("\t".join(newcols)+'\n')
                except ValueError:
                    #print "non numeric value found, skip that PSM"
                    continue;
                except IndexError:
                    continue;

    
    elif n=="2":
        for line in infile2:
            cols=line.strip().split('\t')
            proteinID=cols[1]
            intensity=cols[2:]
            if len(intensity)==samplesize:
                try:
                    ratio=normalize2(intensity)
                    ratio_round=[ '%.3f' % elem for elem in ratio ] #round the value, keep 3 digit after decimal
                
                    cols.insert(1,acclist[proteinID])
                    newcols=cols[:3]+ratio_round
                    outfile.write("\t".join(newcols)+'\n')
                except ValueError:
                    #print "non numeric value found, skip that PSM"
                    continue;
                except IndexError:
                    continue;
                    

    
    print 'file including gene symbol saved'    
    outfile.close()