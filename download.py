import sys
import urllib,urllib2
from Bio import Entrez,SeqIO

def getuniqacc(infile,i): #get a unique accession id list out of file, i is the ith coloum(start from 0) i the file which contain accession id
    dic={}
    infile.readline()
    for line in infile:
        acc=line.split('\t')[i]
        if acc not in dic:
            dic[acc]=1
    
    return dic;
###################Function part end###############################

Entrez.email=raw_input('please provide a valid email address in order to retrieve data from NCBI:(finish by ENTER)\n')

#####################fetch splice variants from EVDB for all genes########
prefix=sys.argv[1]
infilename=prefix+'_pepdata.txt'
pepdatafile=open(infilename,'r')
genedic=getuniqacc(pepdatafile,0)
print 'there are',len(genedic),'genes in total'

pepdatafile.close()

url='http://projects.insilico.us.com/SpliceMiner/Batch?';

variantfile='splicingvar.txt'
output_file1=open(variantfile,'a') #open the file in appending mode

subexonfile='subexon.txt'
file1=open(subexonfile,'r')
downloadgene=getuniqacc(file1,0)
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
vardic=getuniqacc(handle_var,4) #vardic contains all splice variants to be downloaded
handle_var.close()

print len(vardic),"splice variants to be downloaded"

record_dict = SeqIO.index("varseq.fa", "fasta")
var_download={} #var_download contains all splice variants that have been downloaded previously
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
                codon_start=int(seq_feature.qualifiers['codon_start'][0])+seq_feature.location.nofuzzy_start #nofuzzy_start returns an integer instead of a string
                codon_end=seq_feature.location.nofuzzy_end
                try:
                    sequence=seq_feature.qualifiers['translation'][0]
                    
                    output_handle.write(">%s %s|codon_start=%s|codon_end=%s\n%s\n" % (var,record.description,str(codon_start),str(codon_end),sequence))
                    newvar+=1
                except KeyError: #if there is no annotated translation(often occurs when this is a pseudogene), then translate the nucleotide sequence into amino acids.
                    sequence=str(seq_feature.extract(record.seq).translate()) #extract CDS region and translate into amino acids
                    output_handle.write(">%s %s|codon_start=%s|codon_end=%s\n%s\n" % (var,record.description,str(codon_start),str(codon_end),sequence))
                    newvar+=1
            else:
                x=seqtype.index('source')
                seq_feature=record.features[x]
                codon_start=seq_feature.location.nofuzzy_start+1 #traslation starts from the first nucleotide 
                codon_end=seq_feature.location.nofuzzy_end
                sequence=str(seq_feature.extract(record.seq).translate())
                output_handle.write(">%s %s|codon_start=%s|codon_end=%s\n%s\n" % (var,record.description,str(codon_start),str(codon_end),sequence))
                newvar+=1
        except ValueError:
            print var,"not downloaded"
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
