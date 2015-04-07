import sys
import os
import getopt
import numpy as np

def getpep(array): #this function is to find unique peptides in the psmarray
    dic={}
    for i in range(0,len(array)):
        pep=array[i][0].upper()
        if pep not in dic:
            dic[pep]=1
    return dic

def MAD(data,axis=0):
    return np.median(np.absolute(data-np.median(data,axis)),axis)

def main():
    pepdic=getpep(genearray)
    for pep in pepdic.keys():
        pepratio=[]
        for j in range(0,len(genearray)):
            if genearray[j][0].upper()==pep:
                pepratio.append(genearray[j][3:])
    
        try:
            a=np.array(pepratio,dtype=float)
        except ValueError:
            continue;
        if method=="mean":
            mean=np.mean(a,axis=0)
            stdv=np.std(a,axis=0)
            
            mean_round=[ '%.3f' % elem for elem in mean ]
            stdv_round=[ '%.3f' % elem for elem in stdv ]
            
            ensp=pep_ensp[pep]
            newline="%s\t%s\t%s\t%i\t%s\t%s\t%s\n" % (pep,gene,ensp,len(a),','.join(mean_round),
                                                      ','.join(stdv_round),'0')
            output.write(newline)
        elif method=="median":
            median=np.median(a,axis=0)
            mad=MAD(a)
            
            median_round=[ '%.3f' % elem for elem in median ]
            mad_round=[ '%.3f' % elem for elem in mad ]
            
            ensp=pep_ensp[pep]
            newline="%s\t%s\t%s\t%i\t%s\t%s\t%s\n" % (pep,gene,ensp,len(a),','.join(median_round),
                                                      ','.join(mad_round),'0')
            output.write(newline)

if __name__=='__main__':
    ################  Default  ################
    method="median"
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! wrong command, please read the mannual in Readme.txt."
        print "Example: python mergepsm.py --prefix heavy --method median"
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['method=','prefix='])
        for opt, arg in options:
            if opt == '--method': method=arg
            elif opt == '--prefix': prefix=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    infilename=prefix+'_psmdata.txt'
    handle=open(infilename,'r')
    header=handle.readline().split('\t')
    samplesize=len(header[3:])
    if method=="mean":
        newheader=['pep','gene symbol','protein accession','PSM count','mean_ratio','standard deviation','cluster']
    elif method=="median":
        newheader=['pep','gene symbol','protein accession','PSM count','median_ratio','median absolute deviation','cluster']
    firstline='\t'.join(newheader)+'\n'

    outfilename=prefix+'_pepdata.txt'
    output=open(outfilename,'w')
    output.write(firstline)

    gene_psmarray={}
    pep_ensp={}
    for line in handle:
        row=line.strip().replace(",",".").split("\t")
        try:
            pep=row[0].upper()
            gene=row[1]
            ensp=row[2]
            pep_ensp[pep]=ensp
        except IndexError:
            print row
            break;
        if gene=="None":#this will discard PSMs with unknown identity
            continue;
        if ";" in gene: # this will discard PSMs with multiple gene identities
            continue;
        if len(row[3:])!=samplesize:
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
    print "peptides quantities were calculated as %s of PSMs quantities" % method
