import sys

def normalize(lis):
    newlis=[]
    for i in range(0,len(lis)):
        normvalue=float(lis[i])/float(lis[0])
        newlis.append(normvalue)
    
    return newlis

handle=open("Uniproid_Gene.txt",'r')
uniproid_gene={}
for line in handle:
    row=line[:-1].split('\t')
    uniproid_gene[row[0]]=row[1]

msoutput=sys.argv[1]
input_file=open(msoutput,'r')
acclist={}
totalPSM=0
for line in input_file:
    totalPSM+=1
    protein_acc=line.split('\t')[0]
    if protein_acc not in acclist:
        acclist[protein_acc]=1

print totalPSM,'PSMs in this  file'
print len(acclist),'proteins'
input_file.close()

N=0
for protein_acc in acclist.keys():
    proteinList=protein_acc.split(";")
    genename=[]
    for i in range(0,len(proteinList)):
        proteinID=proteinList[i]
        try:
            if '-' not in proteinID:
                genename.append(uniproid_gene[proteinID])
            else:
                index=proteinID.index('-')
                genename.append(uniproid_gene[proteinID[:index]])
        
        except KeyError:
            genename.append('None')
    
    if len(set(genename))==1: #if the proIDs correspond to one gene
        acclist[protein_acc]=genename[0]   #return one gene name
    else: #if not, return gene name for each proID
        acclist[protein_acc]=";".join(genename)
    
    N=N+1
    if N%1000==0:
        print N,'proteins processed'

print N,'proteins processed'
print 'writing output...'


infile2=open(msoutput,'r')

prefix=sys.argv[2]
outfilename=prefix+'_psmdata.txt'
outfile=open(outfilename,'w')

line=infile2.readline()
cols=line.split('\t')
cols.insert(1,'gene symbol')
outfile.write("\t".join(cols))

for line in infile2:
    cols=line.split('\t')
    intensity=cols[2:]
    ratio=normalize(intensity)
    ratio_round=[ '%.2f' % elem for elem in ratio ] #round the value, keep 2 digit after decimal
    proteinID=cols[0]
    cols.insert(1,acclist[proteinID])
    newcols=cols[:3]+ratio_round
    outfile.write("\t".join(newcols)+'\n')

print 'file including gene symbol saved'    
outfile.close()