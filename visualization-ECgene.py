import sys
import os
import getopt
import operator
import numpy
from PIL import Image, ImageDraw, ImageFont
from mapping_EVDB import EXON,ISOFORM,PEPTIDE


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

scale=5 #default, scale 3 base into 1 pixel
format="png" #default image output format
################  Comand-line arguments ################
if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print "Warning! wrong command, please read the mannual in Readme.txt."
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['sample=',
                                                         'scale=',
                                                         'gene-structure-file=',
                                                         'id=','f='])
    for opt, arg in options:
        if opt == '--sample': sample=arg
        elif opt == '--scale': scale=int(arg)
        elif opt == '--gene-structure-file':ECfilename=arg
        elif opt == '--id':id=arg
        elif opt == '--f': format=arg
        else:
            print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()


inputfilename=sample+'_mappingout.txt'
handle1=open(inputfilename)
peparray=[]
#read peptide coordinates from mappingout.txt
samplesize=0
handle1.readline()
for line in handle1:
    row=line[:-1].split("\t")
    gene=row[4].split('.')[0]
    if gene==id:
        peptide=PEPTIDE(seq=row[0],number=row[1],length=row[2],HGNC=row[3],proteinID=row[4],PSMcount=int(row[5]),ratio=row[6],error=row[7],cluster=row[8],variants=row[11].split(","),chr=row[12],start=int(row[13]),end=int(row[14]),strand=row[15],trans_start=int(row[16]),trans_end=int(row[17]),exon1=row[18],exon2=row[19])
        peparray.append(peptide)
        samplesize=len(peptide.ratio.split(","))

handle1.close()

ratio=[]
cluster=[]
for peptide in peparray:
    ratio.append(max(peptide.ratio.split(',')))
    cluster.append(peptide.cluster)

uniq_cluster=list(set(cluster))
print max(ratio),len(uniq_cluster)
ymax=int(float(max(ratio)))+1
yaxis=100+len(uniq_cluster)*120*ymax

variant_objects=[]
variant_exons={} #variant ID as key, value is a dictionary which stores exon object as value,exon number as key
handle2=open(ECfilename ,'r')
strand="+"
for line in handle2:
    row=line.strip().split("\t")
    gene=row[0].split('.')[0]
    if gene==id:
        transcript_length=0
        variant_exons[row[0]]={}
        variant=ISOFORM(id=row[0],chr=row[1],strand=row[2],chr_start=int(row[3]),chr_end=int(row[4]),cds_start=int(row[5]),cds_end=int(row[6]),exon=int(row[7]))
        variant.chr_span=int(row[4])-int(row[3])
        strand=row[2]
        
        startlist=row[8].split(',')
        endlist=row[9].split(',')
        for i in range(len(endlist[:-1])):
            exon=EXON(variant=row[0],chr=row[1],strand=row[2])
            
            if variant.strand=="+":
                exon.number=i+1
                exon.start=int(startlist[i])
                exon.end=int(endlist[i])
            else:
                exon.number=len(endlist[:-1])-i
                exon.start=int(endlist[i])
                exon.end=int(startlist[i])
            
            exon.length=abs(exon.end-exon.start)        
            transcript_length+=exon.length
            variant_exons[row[0]][exon.number]=exon
        
        variant.transcript_length=transcript_length
        variant_objects.append(variant)


############set image parameters#####
min_chr=min(var.chr_start for var in variant_objects)
max_chr=max(var.chr_end for var in variant_objects)
print "maximum chromosome span,bp",(max_chr-min_chr)
setwidth=400+(max_chr-min_chr)/scale
setheight=200+60*len(variant_exons)+60*len(uniq_cluster)+yaxis
print setwidth,setheight
background='#eee' #grey color
im=Image.new('RGB',(setwidth,setheight),background)    
draw=ImageDraw.Draw(im)
font=ImageFont.truetype("Noxchi_Arial.ttf",30)
font2=ImageFont.truetype("Noxchi_Arial.ttf",20)

############draw image of EC splice variants#####
exonheight=50
start=250
introncolor=(125,125,125)
exoncolor=(245,245,220)
stcodon_color='green'
edcodon_color='red'
pep_shadow=(132,112,255)

print len(variant_exons)
print len(variant_objects)
if strand=="+":
    variant_objects.sort(key=operator.attrgetter('chr_start'))
else:
    variant_objects.sort(key=operator.attrgetter('chr_end'),reverse=True)

for j in range(len(variant_objects)):
    var=variant_objects[j]
    print var.id
    y=160+60*j
    y1=y+exonheight/2-2
    draw.text((10,y1),var.id,font=font,fill='black')
    #splice variant exon color
    if j==0:
        indent=0
    else:
        if strand=="+":
            indent=abs(var.chr_start-variant_objects[0].chr_start)/scale
        else:
            indent=abs(var.chr_end-variant_objects[0].chr_end)/scale
        print var.id,indent
    
    stx_chr=variant_exons[var.id][1].start
    print "stx_chr",stx_chr
    for i in range(var.exon):
        exon1=variant_exons[var.id][i+1]
        stx=indent+abs(exon1.start-stx_chr)/scale
        exon1size=exon1.length/scale
        draw.rectangle([start+stx,y,start+stx+exon1size,y+exonheight],fill=exoncolor,outline='black');
        if i+1<var.exon:
            exon2=variant_exons[var.id][i+2]
            stx2=indent+abs(exon2.start-stx_chr)/scale
            draw.rectangle([start+stx+exon1size,y1,start+stx2,y1+4],fill=introncolor,outline=introncolor)
        
        exon1.plt_st=start+stx
        exon1.plt_ed=start+stx+exon1size
        exon1.plt_ycor=y
        
        print exon1.start,start+stx,start+stx+exon1size,i+1,exon1.number,exon1.length,exon1size
    
    
    cds_stx=variant_exons[var.id][1].plt_st+abs(var.cds_start-stx_chr)/scale
    cds_ed=variant_exons[var.id][1].plt_st+abs(var.cds_end-stx_chr)/scale
    if strand=="+":
        draw.rectangle([cds_stx,y,cds_stx+2,y+exonheight],fill=stcodon_color,outline=stcodon_color);
        draw.rectangle([cds_ed,y,cds_ed+2,y+exonheight],fill=edcodon_color,outline=edcodon_color);
    else:
        draw.rectangle([cds_stx,y,cds_stx+2,y+exonheight],fill=edcodon_color,outline=edcodon_color);
        draw.rectangle([cds_ed,y,cds_ed+2,y+exonheight],fill=stcodon_color,outline=stcodon_color);





######draw peptides in different clusters##############
colorlist=[(0,0,205),(65,105,225),(135,206,250),(224,255,255),(0,100,0),
           (34,139,34),(46,139,87),(144,238,144),(32,178,170),(50,205,50)]

panel2=y+80

uniq_cluster.sort()
for i in range(0,len(uniq_cluster)):
    y1=panel2+60*int(uniq_cluster[i])
    clus=int(uniq_cluster[i])
    count=cluster.count(uniq_cluster[i]) #count how many peptide in each cluster
    draw.text((150,y1+10),'%s(%s)'%('cluster'+uniq_cluster[i],str(count)),font=font,fill='black')
    draw.rectangle([10,y1,150,y1+50],fill=colorlist[clus],outline="black")


for peptide in peparray:  #draw peptides for each cluster
    if peptide.start==None:
        continue;
    
    color=colorlist[int(peptide.cluster)]
    for ECvar in peptide.proteinID.split(','):
        exon1=variant_exons[ECvar][int(peptide.exon1)] #exon1 and exon2 are EXON object
        exon2=variant_exons[ECvar][int(peptide.exon2)]
        
        stx1=abs(peptide.start-exon1.start)/scale+exon1.plt_st
        stx2=abs(peptide.end-exon2.start)/scale+exon2.plt_st
        
        y2=panel2+60*int(peptide.cluster)
        draw.text((stx1,y2+10),peptide.number,font=font2,fill='black')
        width=float(peptide.length)/scale
        
        varid=peptide.proteinID
        tx_y2=variant_exons[varid][1].plt_ycor
        if exon1==exon2:
            draw.rectangle([stx1,y2,stx1+width,y2+10],fill=color,outline=color) #draw peptide below splice variants
            draw.rectangle([stx1,tx_y2,stx1+width,tx_y2+exonheight],fill=pep_shadow,outline=pep_shadow)#draw peptide on transcript
        else:
            draw.rectangle([stx1,y2,exon1.plt_ed,y2+10],fill=color,outline="black")
            draw.rectangle([exon1.plt_ed,y2+4,exon2.plt_st,y2+6],fill=introncolor,outline=introncolor)
            draw.rectangle([exon2.plt_st,y2,stx2,y2+10],fill=color,outline="black")
            
            draw.rectangle([stx1,tx_y2,exon1.plt_ed,tx_y2+exonheight],fill=pep_shadow,outline=pep_shadow)
            draw.rectangle([exon2.plt_st,tx_y2,stx2,tx_y2+exonheight],fill=pep_shadow,outline=pep_shadow)
            
            
            print peptide.exon1,peptide.exon2,peptide.start,peptide.end,exon1.start,exon2.start
            print stx1,exon1.plt_ed,exon2.plt_st,stx2



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


imagename=id+'_pattern_'+sample+'.'+format
im.save(imagename,dpi=(300,300))
print imagename+' saved'
