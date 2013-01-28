import sys

def getprotein(infile): #get protein id list out of ms file 
    lis={}
    for line in infile:
        protein_acc=line.split('\t')[0]
        lis[protein_acc]=1
    return lis;

handle=open("Uniproid_Gene.txt",'r')
uniproid_gene={}
for line in handle:
    row=line[:-1].split('\t')
    uniproid_gene[row[0]]=row[1]

msoutput=sys.argv[1]
input_file=open(msoutput,'r')
acclist=getprotein(input_file)
input_file.close()
print len(acclist)
N=0
for protein_acc in acclist.keys():
    proteinList=protein_acc.split(";")
    genename=[]
    for i in range(0,len(proteinList)):
        proteinID=proteinList[i]
        try:
            genename.append(uniproid_gene[proteinID])
        
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
    proteinID=cols[0]
    cols.insert(1,acclist[proteinID])
    outfile.write("\t".join(cols))

print 'file including gene symbol saved'    
outfile.close()