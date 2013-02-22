import sys
import urllib,urllib2
from Bio import Entrez,SeqIO

def getgene(infile): #get gene symbol list out of file 
    dic={}
    infile.readline()
    for line in infile:
        gene=line.split('\t')[0]
        if gene not in dic:
            dic[gene]=1

    return dic;


def getvariant(infile): #get variant acc list out of splicingvar.txt file
    dic={}
    for line in infile:
        array=line.split('\t')
        var_acc=array[4]
        if var_acc not in dic:
            dic[var_acc]=1

    return dic

###################Function part end###############################

Entrez.email=raw_input('please provide a valid email address in order to retrieve data from NCBI:(finish by ENTER)\n')

#####################fetch splice variants from EVDB for all genes########
prefix=sys.argv[1]
infilename=prefix+'_pepdata.txt'
pepdatafile=open(infilename,'r')
genedic=getgene(pepdatafile)
print 'there are',len(genedic),'genes in total'

pepdatafile.close()

url='http://projects.insilico.us.com/SpliceMiner/Batch?';

variantfile='splicingvar.txt'
output_file1=open(variantfile,'a') #open the file in appending mode

subexonfile='subexon.txt'
file1=open(subexonfile,'r')
downloadgene=getgene(file1)
print len(downloadgene),'genes exon composition stored locally'


output_file2=open(subexonfile,'a')


output_file3=open('gene_notfound.txt','a')

newgene=0
i=0
for gene in genedic.keys():
    i+=1
    if gene not in downloadgene:
        values={'queryBy':'Gene',
                'organism':'9606',
                'inputType':'Text',
                'text':gene,
                }

        data=urllib.urlencode(values)
        response=urllib.urlopen(url+data)
        response.readline()
        s=response.read()
        if s[0:10]=='No results':
            output_file3.write(gene+'\n')       
        else:
            newgene+=1
            output_file1.write(s)

            values['resultType']='subexon'

            newdata=urllib.urlencode(values)
            newresponse=urllib.urlopen(url+newdata)
            newresponse.readline()
            output_file2.write(newresponse.read())
    if i%1000==0:
        print i,'\tgenes processed'

print len(genedic),'genes processed'
print newgene,'new genes exon composition downloaded'

output_file2.close()
print subexonfile,'saved'
output_file1.close()
print variantfile,'saved'
output_file3.close()
print 'gene_notfound.txt saved'

############fetch translated sequence from NCBI for all splicing variants#####
handle_var=open(variantfile,'r')
vardic=getvariant(handle_var)
handle_var.close()

print len(vardic),"splice variants to be downloaded"

record_dict = SeqIO.index("varseq.fa", "fasta")
var_download={}
for key in record_dict.keys():
	var_download[key]=1

file4=open('var_notfound.txt','a')
output_handle=open('varseq.fa','a')

newvar=0
icount=0
for var in vardic.keys():
    icount+=1
    if icount%1000==0:
        print icount,"known splice varaints in this set of genes "
    if var not in var_download:
        try:
            handle=Entrez.efetch(db='nucleotide',
                                 id=var,
                                 rettype='gb',
                                 retmode='text',
                                 tool='biopython')

            record=SeqIO.read(handle,'gb')
            seqtype=[] #gather sequence feature type
            for seq_feature in record.features:
                seqtype.append(seq_feature.type)
                
            if 'CDS' in seqtype:
                x=seqtype.index('CDS')
                seq_feature=record.features[x]
                try:
                    sequence=seq_feature.qualifiers['translation'][0]
                    output_handle.write(">%s\t%s\t%s\n%s\n" % (
                    var,seq_feature.location,record.description,sequence))
                    newvar+=1
                except KeyError: #if there is no annotated translation(often occurs when this is a pseudogene), then translate the nucleotide sequence into amino acids.
                    sequence=str(seq_feature.extract(record.seq).translate())
                    output_handle.write(">%s\t%s\t%s\n%s\n" % (
                    var,seq_feature.location,record.description,sequence))
            else:
                x=seqtype.index('source')
                seq_feature=record.features[x]
                sequence=str(seq_feature.extract(record.seq).translate())
                output_handle.write(">%s\t%s\t%s\n%s\n" % (
                    var,seq_feature.location,record.description,sequence))
        except urllib2.HTTPError:
            file4.write('%s\n'%(var))
            continue ## if HTTPError occurs, continue fetching from next one.
        except httplib.BadStatusLine:
            time.sleep(1)
            file4.write('%s\n'%(var))
            continue	            

print len(vardic),"splice variants processed"
print newvar,'new splice variants sequence downloaded'
output_handle.close()
file4.close()
print 'varseq.fa\tsaved'
