import sys
from Bio import SeqIO
from cogent.db.ensembl import Genome

print "reading input..."

human=Genome('human',Release=63)

handle = open("Homo_sapiens.GRCh37.63.pep.all.fa", "rU")#fetch geneID by proID
record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()

gene_symbol={}
handle2=open('proteincoding_genesymbol.txt','r')

for line in handle2:
    row=line[:-1].split('\t')
    symbol=row[1]
    geneid=row[0]
    if geneid not in gene_symbol:
        gene_symbol[geneid]=symbol
    

pqpqoutput=sys.argv[1]


infile=open(pqpqoutput,'r')

infile.readline()
ENSPlist={}
totalPSM=0

for line in infile:
    cols=line.split('\t')
    if cols[0] not in ENSPlist:
        ENSPlist[cols[0]]=1

    totalPSM+=1

print totalPSM,'PSMs in this  file'
print len(ENSPlist),'proteins'
infile.close()

N=0
for protein_acc in ENSPlist.keys():
    proteinList=protein_acc.split(";")
    genename=[]
    for i in range(0,len(proteinList)):
        proteinID=proteinList[i]
        descrip=record_dict[proteinID].description
        ind=descrip.index('gene:')
        geneID=descrip[ind+5:ind+20]
        try:
            genename.append(gene_symbol[geneID])
        except KeyError:
            gene = human.getGeneByStableId(StableId=geneID)
            genename.append(gene.Symbol)
        except AttributeError:
            genename.append('Not Found')

    if len(set(genename))==1: #if the proIDs correspond to one gene
        ENSPlist[protein_acc]=genename[0]   #return one gene name
    else:                   #if not, return gene name for each proID
        ENSPlist[protein_acc]=";".join(genename)

    N=N+1
    if N%1000==0:
        print N,'proteins processed'

print N,'proteins processed'
print 'writing output...'
infile2=open(pqpqoutput,'r')

prefix=sys.argv[2]
outfilename=prefix+'_pqpqout.txt'
outfile=open(outfilename,'w')

line=infile2.readline()
cols=line.split('\t')
cols.insert(1,'gene name')
outfile.write("\t".join(cols))

for line in infile2:
    cols=line.split('\t')
    proteinID=cols[0]
    cols.insert(1,ENSPlist[proteinID])
    outfile.write("\t".join(cols))

print 'file including gene symbol saved'    
outfile.close()
