import sys
import os
import getopt
import operator
from PIL import Image, ImageDraw, ImageFont
from mapping_EVDB import EXON,ISOFORM,PEPTIDE
import numpy

def drawxy(ymax):#ymax is the maximum peptide ratio
    draw.text((200,axis_start),'ratio',font=font,fill='black')
    draw.rectangle([350,axis_start,360,axis_start+100*ymax],fill='black',outline='black')
    draw.rectangle([360,axis_start+100*ymax,setwidth-200,axis_start+100*ymax-10],fill='black',outline='black')
    for i in range(1,ymax+1):# draw grid on y axis
        y3=axis_start+100*ymax-100*i
        draw.text((300,y3),str(i),font=font,fill='black')
        draw.rectangle([350,y3,390,y3+10],fill='black',outline='black')

def histogram(peptide):#take in peptide object
    a=[peptide.ratio.split(","),peptide.error.split(",")]
    b=numpy.array(a,dtype=float)
    pattern=numpy.transpose(b)
    color=colorlist[int(peptide.cluster)]
    for i in range(0,len(pattern)):
        draw.rectangle([barstart+25*i,y4,barstart+25*(i+1),y4-100*pattern[i][0]],fill=color,outline='black')
        draw.rectangle([barstart+10+25*i,y4-100*(pattern[i][0]+pattern[i][1]),barstart+15+25*i,y4-100*(pattern[i][0]-pattern[i][1])],fill='black',outline='black')


###################Function part end###############################
format="png" #default image output format

################  Comand-line arguments ################
if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print "Warning! wrong command, please read the mannual in Readme.txt."
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['f=','prefix=','gene=',])
    for opt, arg in options:
        if opt == '--prefix': sample=arg
        elif opt == '--f': format=arg
        elif opt == '--gene': gene=arg
        else:
            print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
gene=gene.upper()

inputfilename=sample+'_mappingout.txt'
handle1=open(inputfilename)
peparray=[]
#read peptide coordinates from mappingout.txt
samplesize=0
for line in handle1:
    row=line[:-1].split("\t")
    if row[3].upper()==gene:
        peptide=PEPTIDE(seq=row[0],number=row[1],length=row[2],HGNC=row[3],proteinID=row[4],PSMcount=int(row[5]),ratio=row[6],error=row[7],cluster=row[8],variants=row[11].split(","),chr=row[12],start=int(row[13]),
                        end=int(row[14]),strand=row[15],trans_start=int(row[16]),trans_end=int(row[17]),exon1=row[18],exon2=row[19])
        peparray.append(peptide)
        samplesize=len(peptide.ratio.split(","))

handle1.close()

if len(peparray)==0:
    print "%s not in the mappingout.txt file. Exiting..." % gene; sys.exit()

print len(peparray),"peptides identified from gene",gene

ratio=[]
cluster=[]
for peptide in peparray:
    ratio.append(max(peptide.ratio.split(',')))
    cluster.append(peptide.cluster)

uniq_cluster=list(set(cluster))
print "maximum peptide ratio",max(ratio)
print len(uniq_cluster),"unique peptide cluster discovered for this gene"

ymax=int(float(max(ratio)))+1
yaxis=100+len(uniq_cluster)*120*ymax

#read splice variant coordinates from file
handle2=open('splicingvar.txt')
variant_exon={}
for line in handle2:
    row=line[:-1].split("\t")
    if row[0].upper()==gene:
        exon=EXON(gene=row[0].upper(),chr=row[2],strand=row[3],variant=row[4],number=int(row[5]),
                  start=int(row[6]),end=int(row[7]),trans_start=int(row[8]),trans_end=int(row[9]))
        
        variantID=row[4]
        if variantID not in variant_exon:
            variant_exon[variantID]=[exon]
        else:
            variant_exon[variantID].append(exon)

handle2.close()

if len(variant_exon)==0:
    print "%s doesn't have splice variant structure downloaded. Exiting..." % gene; sys.exit()

print len(variant_exon),"known splice variants for gene",gene

#read gene subexons coordinates from file
handle3=open('subexon.txt')
gene_subexon=[]
exon_chr={}
for line in handle3:
    row=line[:-1].split("\t")
    if row[0].upper()==gene:
        subexon=EXON(gene=row[0].upper(),chr=row[2],strand=row[3],number=int(row[4]),start=int(row[6]),end=int(row[7]))
        gene_subexon.append(subexon)
        if subexon.number not in exon_chr:
            exon_chr[subexon.number]=[subexon.start,subexon.end]
        else:
            exon_chr[subexon.number][1]=subexon.end

handle3.close()

if len(gene_subexon)==0:
    print "%s doesn't have subexon structure downloaded. Exiting..." % gene; sys.exit()

#set size for the image
max_transcript=0 #maximum transcript length
for exon in exon_chr:
    max_transcript+=abs(exon_chr[exon][0]-exon_chr[exon][1])

uniqcount=[] #count the unique peptides for each cluster
for i in range(0,len(uniq_cluster)):
    count=cluster.count(uniq_cluster[i])
    uniqcount.append(count)

start=400 #start position for first subexon
height=50
background='#eee' #grey color
setwidth1=start+max_transcript+20*len(exon_chr)
setwidth2=400+(samplesize*25+35)*max(uniqcount)

setwidth=max(setwidth1,setwidth2)
setheight=200+60*len(variant_exon)+60*len(uniq_cluster)+yaxis
im=Image.new('RGB',(setwidth,setheight),background)    
draw=ImageDraw.Draw(im)
font=ImageFont.truetype("Noxchi_Arial.ttf",30)
font2=ImageFont.truetype("Noxchi_Arial.ttf",20)


################draw image of subexons###########
draw.text((10,50),gene,font=font,fill='black')
draw.text((10,100),'Subexons: ',font=font,fill='black')

width=0
stx=0
stx2=0
y=100 #y axis of the subexon

exon_plt_start={} #subexon start position in the plot
exon_plt_end={} #subexon end position in the plot
intron_color=(125,125,125)
for i in range(0,len(gene_subexon)):
    subexon=gene_subexon[i]
    stx1=(subexon.number-1)*20
    stx=stx1+stx2
    width=abs(subexon.start-subexon.end) #width of subexon
    y1=y+height/2-5 #y axix of the intron line
    
    exon_plt_end[subexon.number]=start+stx+width
    if subexon.number not in exon_plt_start:
        exon_plt_start[subexon.number]=start+stx
    
    draw.rectangle([start+stx,y,start+stx+width,y+height],fill='white',outline='black');
    stx2+=width
    if i!=len(gene_subexon)-1:
        draw.rectangle([start+stx+width,y1,start+stx+width+20,y1+10],fill=intron_color,outline=intron_color);


#label the subexons
for exon in exon_plt_start:
    x_cor=(exon_plt_start[exon]+exon_plt_end[exon])/2
    y_cor=50
    draw.text((x_cor,y_cor),str(exon),font=font,fill='blue');
############draw image of known splicing variants#####
j=0
variant_plt_start={} #splice variant's all exons' start position in the plot
variant_plt_end={} #splice variant's all exons' end position in the plot
for variant in variant_exon.keys():
    variant_plt_start[variant]={}
    variant_plt_end[variant]={}
    exonarray=variant_exon[variant]
    exonarray.sort(key=operator.attrgetter('trans_start'))
    color=(245,245,220) #splice variant exon color
    y=160+60*j
    y1=y+height/2-5
    draw.text((10,y),variant,font=font,fill='black')
    for i in range(0,len(exonarray)):
        exon=exonarray[i]
        if exon.trans_start>0:
            stx=abs(exon.start-exon_chr[exon.number][0])+exon_plt_start[exon.number]
            width=abs(exon.start-exon.end)
            draw.rectangle([stx,y,stx+width,y+height],fill=color,outline='black');
            
            variant_plt_end[variant][exon.number]=stx+width
            if exon.number not in variant_plt_start[variant]:
                variant_plt_start[variant][exon.number]=stx
            
            if i!=len(exonarray)-1:
                exon2=exonarray[i+1]
                stx2=abs(exon2.start-exon_chr[exon2.number][0])+exon_plt_start[exon2.number]
                width2=abs(exon2.start-exon2.end)
                draw.rectangle([stx+width,y1,stx2,y1+10],fill=intron_color,outline=intron_color);
    j+=1


######draw peptides in different clusters##############
colorlist=[(255,250,250),(0,0,205),(65,105,225),(135,206,250),(224,255,255),(0,100,0),
           (34,139,34),(46,139,87),(144,238,144),(32,178,170),(50,205,50)]

panel2=y+60

uniq_cluster.sort()
for i in range(0,len(uniq_cluster)):
    y1=panel2+60*int(uniq_cluster[i])
    clus=int(uniq_cluster[i])
    count=cluster.count(uniq_cluster[i]) #count how many peptide in each cluster
    draw.text((150,y1+10),'%s(%s)'%('cluster'+uniq_cluster[i],str(count)),font=font,fill='black')
    draw.rectangle([10,y1,150,y1+50],fill=colorlist[clus],outline="black")


for peptide in peparray:  #draw peptides for each cluster
    if peptide.start==0:
        print "Omit peptide %s which doesn't map to any of splice variants to this gene" % peptide.seq
        continue;
        
    color=colorlist[int(peptide.cluster)]
    exon1=int(peptide.exon1)
    exon2=int(peptide.exon2)
    stx1=abs(peptide.start-exon_chr[exon1][0])+exon_plt_start[exon1]
    stx2=abs(peptide.end-exon_chr[exon2][0])+exon_plt_start[exon2]
    y2=panel2+60*int(peptide.cluster)
    draw.text((stx1,y2+10),peptide.number,font=font2,fill='black')
    if exon1==exon2:
        draw.rectangle([stx1,y2,stx2,y2+10],fill=color,outline=color)
    else:
        if peptide.variants[0]=="": #for novel splice junction peptides which can't map to any of known variants
            draw.rectangle([stx1,y2,stx1+peptide.trans_start,y2+10],fill=color,outline="black")
            draw.rectangle([stx1+peptide.trans_start,y2+4,stx2-peptide.trans_end,y2+6],fill=intron_color,outline=intron_color)
            draw.rectangle([stx2-peptide.trans_end,y2,stx2,y2+10],fill=color,outline="black")
            print peptide.seq,peptide.number,stx1,stx1+peptide.trans_start,stx2-peptide.trans_end,stx2
        else:
            variant_exon2_start=variant_plt_start[peptide.variants[0]][exon2]
            variant_exon1_end=variant_plt_end[peptide.variants[0]][exon1]
            draw.rectangle([stx1,y2,variant_exon1_end,y2+10],fill=color,outline="black")
            draw.rectangle([variant_exon1_end,y2+4,variant_exon2_start,y2+6],fill=intron_color,outline=intron_color)
            draw.rectangle([variant_exon2_start,y2,stx2,y2+10],fill=color,outline="black")


############draw quantitative pattern of the peptides###############
for i in range(0,len(uniq_cluster)):
    axis_start=y1+100+(100*ymax+20)*i
    y4=axis_start+100*ymax
    drawxy(ymax)
    draw.text((150,(y4+axis_start)/2),'cluster'+str(uniq_cluster[i]),font=font,fill='black')
    j=0
    for peptide in peparray:
        if peptide.cluster==uniq_cluster[i]:
            barstart=400+j*(samplesize*25+30)
            mid=samplesize*25/4
            histogram(peptide)
            xlabel='peptide '+peptide.number+'('+str(peptide.PSMcount)+')'
            draw.text((barstart+mid,y4),xlabel,font=font2,fill='black')
            j+=1


###draw mapped peptide on the transcript######
im2=Image.new('RGB',(setwidth,setheight),'#eee')    
draw2=ImageDraw.Draw(im2)

k=0
for variant in variant_exon.keys():
    y=160+60*k
    for peptide in peparray:
        if peptide.start==0:
            continue;
        
        color=(132,112,255)
        MTV=peptide.variants
        exon1=int(peptide.exon1)
        exon2=int(peptide.exon2)
        if variant in MTV:
            #print array2[i][4],var[k],MTV
            stx1=abs(peptide.start-exon_chr[exon1][0])+exon_plt_start[exon1]
            stx2=abs(peptide.end-exon_chr[exon2][0])+exon_plt_start[exon2]
            if exon1==exon2:
                draw2.rectangle([stx1,y,stx2,y+height],fill=color)
            else:
                variant_exon2_start=variant_plt_start[variant][exon2]
                variant_exon1_end=variant_plt_end[variant][exon1]
                draw2.rectangle([stx1,y,variant_exon1_end,y+height],fill=color,outline=color)
                draw2.rectangle([variant_exon2_start,y,stx2,y+height],fill=color,outline=color)

    k+=1


newimage=Image.blend(im,im2,0.3)

imagename=gene+'_pattern_'+sample+'.'+format
newimage.save(imagename,dpi=(300,300))
print imagename+' saved'
