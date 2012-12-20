import sys
from Bio import SeqIO
from cogent.db.ensembl import Genome

handle = open("Homo_sapiens.GRCh37.63.pep.all.fa", "rU")#fetch geneID by proID
record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()

pro_gene={}
for proID in record_dict.keys():
	descrip=record_dict[proID].description
	ind=descrip.index("gene:")
	geneID=descrip[ind+5:ind+20]
	pro_gene[proID]=geneID

infile=open(sys.argv[1],'r')

gene_symbol={}
for line in infile:
	row=line[:-1].split('\t')
	geneid=row[0]
	gene_symbol[geneid]=row[1]

outfile=open(sys.argv[2],'w')

human=Genome('human',Release=63)
for pro in pro_gene.keys():
	gene=pro_gene[pro]
	if gene in gene_symbol:
		line="%s\t%s\t%s\n" % (pro,gene,gene_symbol[gene])
		outfile.write(line)
	else:
		entry = human.getGeneByStableId(StableId=gene)
		line="%s\t%s\t%s\n" % (pro,gene,entry.Symbol)

infile.close()
outfile.close()
	
