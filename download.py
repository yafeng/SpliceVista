import sys
import os
import getopt
import urllib,urllib2
import httplib
from Bio import Entrez,SeqIO

def getuniqacc(infile,i): #get a unique accession id list out of file, i is the ith coloum(start from 0) i the file which contain accession id
    dic={}
    infile.readline()
    for line in infile:
        acc=line.split('\t')[i]
        if acc not in dic:
            dic[acc]=1
    
    return dic;

def getTrans_Start(var):
    start_trans=""
    query=">"+var+"\n"+sequence
    values={'queryBy':'Protein',
        'organism':organism,
        'inputType':'Text',
        'text':query,
    }
    
    data=urllib.urlencode(values)
    response=urllib.urlopen(url+data)
    response.readline()
    s=response.read()
    if s[0:10]=='No results':
        pass       
    else:
        lines=s.split("\n")
        for i in range(0,len(lines)):
            row=lines[i].split("\t")
            try:
                if row[4]==var:
                    start_trans=row[12]
                    break
            except IndexError:
                print var,"Start(Trans) empty"
    
    return start_trans

###################Function part end###############################

################  Default  ################
organism = '9606' #H.sapiens

################  Comand-line arguments ################
if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print "Warning! wrong command, please read the mannual in Readme.txt."
    print "Example: python download.py --organism 9606 --prefix heavy --email youremailaddress@XXX"
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['organism=',
                                                         'email=',
                                                         'prefix='])
    for opt, arg in options:
        if opt == '--organism': organism=arg
        elif opt == '--email': Entrez.email=arg
        elif opt == '--prefix':prefix=arg
        else:
            print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()


#####################fetch splice variants from EVDB for all genes########
infilename=prefix+'_pepdata.txt'
pepdatafile=open(infilename,'r')
genedic=getuniqacc(pepdatafile,1)
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

file3name=prefix+'_gene_notfound.txt'
output_file3=open(file3name,'w')

newgene=0
i=0
for gene in genedic.keys():
    i+=1
    if gene not in downloadgene:
        values={'queryBy':'Gene',
                'organism':organism,
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

print "Downloading sequences of splice variant from NCBI GenBank"

record_dict = SeqIO.index("varseq.fa", "fasta")
var_download={} #the dictionary contains all splice variants that have been downloaded previously
for key in record_dict.keys():
	var_download[key]=1

output_handle=open('varseq.fa','a')

newvar=0
icount=0
for var in vardic.keys():
    icount+=1
    if icount%1000==0:
        print icount,"splice variants processed"
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
                    trans_st=getTrans_Start(var)
                    output_handle.write(">%s %s|codon_start=%s|codon_end=%s|Start(Trans)=%s\n%s\n" % (var,record.description,str(codon_start),str(codon_end),trans_st,sequence))
                    newvar+=1
                except KeyError: #if there is no annotated translation(often occurs when this is a pseudogene), then translate the nucleotide sequence into amino acids.
                    sequence=str(seq_feature.extract(record.seq).translate()) #extract CDS region and translate into amino acids
                    trans_st=getTrans_Start(var)
                    output_handle.write(">%s %s|codon_start=%s|codon_end=%s|Start(Trans)=%s\n%s\n" % (var,record.description,str(codon_start),str(codon_end),trans_st,sequence))
                    newvar+=1
            else:
                x=seqtype.index('source')
                seq_feature=record.features[x]
                codon_start=seq_feature.location.nofuzzy_start+1 #traslation starts from the first nucleotide 
                codon_end=seq_feature.location.nofuzzy_end
                sequence=str(seq_feature.extract(record.seq).translate())
                trans_st=getTrans_Start(var)
                output_handle.write(">%s %s|codon_start=%s|codon_end=%s|Start(Trans)=%s\n%s\n" % (var,record.description,str(codon_start),str(codon_end),trans_st,sequence))
                newvar+=1
        except IndexError:
            var,"Start(Trans) empty"
        except ValueError:
            print var,"not downloaded"
        except urllib2.HTTPError:
            print var,"not downloaded"
        except httplib.BadStatusLine:
            time.sleep(1)
            print var,"not downloaded"
                        

print len(vardic),"splice variants sequences stored locally"
print newvar,'new splice variants downloaded'
output_handle.close()
print 'varseq.fa\tsaved'
