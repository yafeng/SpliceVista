import sys
import urllib,urllib2
from Bio import Entrez,SeqIO
from collections import OrderedDict

def getgene(infile): #get gene symbol list out of pepdata file 
    dic={}
    infile.readline()
    for line in infile:
        gene=line.split('\t')[1]
        if ";" not in gene and gene not in dic: #abandon PSMs mapped to multiple genes
            dic[gene]=1

    return dic;

def getdownloadgene(infile): #get downloaded gene list out of splicingvar file 
    dic={}
    for line in infile:
        gene=line.split('\t')[0]
        if gene not in dic:
            dic[gene]=1

    return dic;

def getvariant(infile): #get variant acc list out of splicevariant file
    dic={}
    for line in infile:
        array=line.split('\t')
        acc=array[4]
        if acc not in dic:
            dic[acc]=1
    return dic

###################Function part end###############################


#####################fetch splicing variants from EVDB for all genes########
prefix=sys.argv[1]
infilename=prefix+'_pepdata.txt'
pepdatafile=open(infilename,'r')
genelist=getgene(pepdatafile).keys()
print 'there are',len(genelist),'genes in total' #after discard PSMs mapped to multiple genes

pepdatafile.close()

url='http://projects.insilico.us.com/SpliceMiner/Batch?';

variantfile='splicingvar.txt'
output_file1=open(variantfile,'a')

subexonfile='subexon.txt'
file1=open(subexonfile,'r')
downloadgene=getdownloadgene(file1)
print len(downloadgene),'genes exon structure stored locally'


output_file2=open(subexonfile,'a')

file3=open('gene_notfound.txt','r')
notfounddict={}
for line in file3:
    if line!='' and line[:-1] not in notfounddict:
        notfounddict[line[:-1]]=1 #get a genelist that are not found in EVDB

file3.close()

newgene=0

for i in range(0,len(genelist)):
    if genelist[i] not in downloadgene:
        values={'queryBy':'Gene',
                'organism':'9606',
                'inputType':'Text',
                'text':genelist[i],
                }

        data=urllib.urlencode(values)

        response=urllib.urlopen(url+data)
        response.readline()
        s=response.read()
        if s[0:10]=='No results':
            if genelist[i] not in notfounddict:
                notfounddict[genelist[i]]=1        
        else:
            newgene+=1
            print genelist[i]
            output_file1.write(s)

            values['resultType']='subexon'

            newdata=urllib.urlencode(values)
            newresponse=urllib.urlopen(url+newdata)
            newresponse.readline()
            output_file2.write(newresponse.read())
    if i%1000==0:
        print genelist[i], i,'\tgenes processed'

print len(genelist),'genes processed'
print newgene,'new genes exon structure downloaded'

output_file2.close()
print subexonfile,'saved'
output_file1.close()
print variantfile,'saved'
output_file3=open('gene_notfound.txt','w')
for key in notfounddict.keys():
    output_file3.write(key+'\n')

output_file3.close()
print 'gene_notfound.txt saved'

############fetch translated sequence from NCBI for all splicing variants#####
output_file1=open(variantfile,'r')
varlist=getvariant(output_file1).keys()
print 'there are',len(varlist),'splice variants for those genes in your input file'
output_file1.close()
Entrez.email=sys.argv[2]

record_dict = SeqIO.index("varseq.fa", "fasta")
print len(record_dict),'splice variants sequence have been downloaded in the varseq.fa'
var_download={}
for key in record_dict.keys():
	var_download[key]=1
	
file4=open('var_notfound.txt','a')
output_handle=open('varseq.fa','a')

newvar=0

for i in range(0,len(varlist)):
    if varlist[i] not in var_download:
        try:
            handle=Entrez.efetch(db='nucleotide',
                                 id=varlist[i],
                                 rettype='gb',
                                 retmode='text',
                                 tool='biopython')

            record=SeqIO.read(handle,'gb')
            seqtype=[]
            for seq_feature in record.features:
                seqtype.append(seq_feature.type)
                
            if 'CDS' in seqtype:
                x=seqtype.index('CDS')
                seq_feature=record.features[x]
                try:
                    sequence=seq_feature.qualifiers['translation'][0]
                    output_handle.write(">%s\t%s\t%s\n%s\n" % (
                    varlist[i],seq_feature.location,record.description,sequence))
                    newvar+=1
                except KeyError: #if there is no annotated translation(often occurs when this is a pseudogene), then generate artificial translation.
                    sequence=str(seq_feature.extract(record.seq).translate())
                    output_handle.write(">%s\t%s\t%s\n%s\n" % (
                    varlist[i],seq_feature.location,record.description,sequence))
            else:
                x=seqtype.index('source')
                seq_feature=record.features[x]
                sequence=str(seq_feature.extract(record.seq).translate())
                output_handle.write(">%s\t%s\t%s\n%s\n" % (
                    varlist[i],seq_feature.location,record.description,sequence))
	except urllib2.HTTPError:
            file4.write('%s\n'%(varlist[i]))
            continue ## if HTTPError occurs, continue fetching from next one.
        except httplib.BadStatusLine:
            time.sleep(1)
            file4.write('%s\n'%(varlist[i]))
            continue	            
    if i>0 and i%1000==0:
        print i,'\tsplice variants processed'

print len(varlist),'splice variants processed'
print newvar,'new splice variants sequence downloaded'
output_handle.close()
file4.close()
print 'varseq.fa\tsaved'
