import sys
from Bio import SeqIO

def normalize(lis):
    newlis=[]
    for i in range(0,len(lis)):
        normvalue=2*float(lis[i])/(float(lis[0])+float(lis[1]))
        newlis.append(normvalue)
    
    return newlis


print "reading input..."

protein_symbol={}
handle=open('ENSP_ENSG_genesymbol.txt','r')

for line in handle:
    row=line[:-1].split('\t')
    symbol=row[2]
    proid=row[0]
    if proid not in protein_symbol:
        protein_symbol[proid]=symbol
    
msoutput=sys.argv[1]

infile=open(msoutput,'r')

infile.readline()
ENSPlist={}
totalPSM=0

for line in infile:
    totalPSM+=1
    cols=line.split('\t')
    if cols[0] not in ENSPlist:
        ENSPlist[cols[0]]=1

    		
print totalPSM,'PSMs in this  file'
print len(ENSPlist),'proteins'
infile.close()

N=0
for protein_acc in ENSPlist.keys():
    proteinList=protein_acc.split(";")
    genename=[]
    for i in range(0,len(proteinList)):
        proteinID=proteinList[i]
        try:
            genename.append(protein_symbol[proteinID])

        except KeyError:
            genename.append('Not Found')

    if len(set(genename))==1: #if the proIDs correspond to one gene
        ENSPlist[protein_acc]=genename[0]   #return one gene name
    else: #if not, return gene name for each proID
        ENSPlist[protein_acc]=";".join(genename)

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
    cols=line[:-1].split('\t')
    intensity=cols[2:]
    ratio=normalize(intensity)
    ratio_round=[ '%.2f' % elem for elem in ratio ] #round the value, keep 2 digit after decimal
    proteinID=cols[0]
    cols.insert(1,ENSPlist[proteinID])
    newcols=cols[:3]+ratio_round
    outfile.write("\t".join(newcols)+'\n')

print 'file including gene symbol saved'    
outfile.close()
