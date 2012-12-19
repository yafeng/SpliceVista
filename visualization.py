#! /usr/bin/env python
import sys
import Image, ImageDraw, ImageFont
from collections import OrderedDict
from math import sqrt
import operator

def extract(gene,infile):#extract splicevariant, subexon info for one gene
    array=[]
    for line in infile:
        if line.split('\t')[0]==gene:
            array.append(line[:-2].split('\t'));

    return array

def extractmap(gene,infile):#extract mappingout info for one gene
    array=[]
    for line in infile:
        if line.split('\t')[1]==gene:
            array.append(line[:-1].split('\t'));

    return array  
#######################Subexons structure############
def Exoncordinate(lis): #returns a list of cordinates of all exons
    cor=[]
    for i in range(0,len(lis)):
            cor.append((int(lis[i][6]),int(lis[i][7])))

    return cor;


def Exonorder(lis):# returns a list of exonorder of all exons
    order=[]
    for i in range(0,len(lis)):
        order.append(int(lis[i][4]))

    return order

def Exonstart(lis): #returns a list of genomic cordinates of starting site of each exon
    cor=[]
    cor.append(int(lis[0][6]))
    for i in range(1,len(lis)):
        if lis[i-1][4]!=lis[i][4]:
            cor.append(int(lis[i][6]));

    return cor;

def Exonsend(lis): #returns a list of genomic cordinates of starting site of each exon
    cor=[]
    cor.append(int(lis[0][7]))
    for i in range(1,len(lis)):
        if lis[i-1][4]!=lis[i][4]:
            cor.append(int(lis[i][7]));

    return cor;

def sumlen(index,lis): #sum all elements before the element lis[index]
    sumlen=0
    for i in range(0,index):
        sumlen+=max(lis[i])-min(lis[i]);

    return sumlen;

########################Splicing Variants Structure#########
def collectvariant(lis): #returns a list of variants the gene has.
    varlist=[]
    for i in range(0,len(lis)):
        if lis[i][4] not in varlist:
            varlist.append((lis[i][4]));
    return varlist;

def var_cor(lis1,lis2):   #returns an array which stores exons chr_cor from one var together 
    cordinate=[]
    for i in range(0,len(lis1)):
        cordinate.append([])
        for j in range(0,len(lis2)):
            if lis1[i]==lis2[j][4]:
                cordinate[i].append((int(lis2[j][6]),int(lis2[j][7])));

    return cordinate;

def Exon_cor(lis1,lis2):   #returns an array which stores exons aa_cor from one var together 
    cordinate=[]
    for i in range(0,len(lis1)):
        cordinate.append([])
        for j in range(0,len(lis2)):
            if lis1[i]==lis2[j][4]:
                cordinate[i].append((int(lis2[j][8]),int(lis2[j][9])));

    return cordinate;

def Exon_no(lis1,lis2):  #returns an list which stores exon position of variants
    cordinate=[]
    for i in range(0,len(lis1)):
        cordinate.append([])
        for j in range(0,len(lis2)):
            if lis1[i]==lis2[j][4]:
                cordinate[i].append(int(lis2[j][5]));

    return cordinate;

def overlap(t1,t2):# check if two tuples overlap
    if min(t2)<min(t1)<max(t2) or min(t2)<max(t1)<max(t2):
        return True;
    else:
        return False;
    
def getcolor(string,lis1,colorlist):
    i=lis1.index(string)
    return colorlist[i];

def drawxy(ymax):#ymax is the maximum peptide ratio
    draw.text((100,axis_start),'Foldchange',font=font,fill='black')
    draw.rectangle([350,axis_start,360,axis_start+100*ymax],fill='black',outline='black')
    draw.rectangle([360,axis_start+100*ymax,setwidth-200,axis_start+100*ymax-10],fill='black',outline='black')
    for i in range(1,ymax+1):# draw grid on y axis
        y3=axis_start+100*ymax-100*i
        draw.text((300,y3),str(i),font=font,fill='black')
        draw.rectangle([350,y3,390,y3+10],fill='black',outline='black')

def meanstdv(x):
    n,mean,std=len(x),0,0
    for a in x:
        mean=mean+a
    mean=mean/float(n)
    if n>1:
        for a in x:
            std=std+(a-mean)**2
        std=sqrt(std/float(n-1))

    return mean,std

def histogram(pep,cluster):
    h0,h2,h6,h24=[],[],[],[]
    h0r,h2r,h6r,h24r=[],[],[],[]
    pattern=[]
    for i in range(0,len(array2)):
        if array2[i][2].upper()==pep and int(array2[i][6])==cluster:
            pepratio=array2[i][7].split(',')
            h0.append(float(pepratio[0]))
            h0r.append(float(pepratio[1]))
            h2.append(float(pepratio[2]))
            h2r.append(float(pepratio[3]))
            h6.append(float(pepratio[4]))
            h6r.append(float(pepratio[5]))
            h24.append(float(pepratio[6]))
            h24r.append(float(pepratio[7]))

    #pattern=[(h0mean,h0std),(h0rmean,h0rstd),...]
    pattern.append(meanstdv(h0)) 
    pattern.append(meanstdv(h0r))
    pattern.append(meanstdv(h2))
    pattern.append(meanstdv(h2r))
    pattern.append(meanstdv(h6))
    pattern.append(meanstdv(h6r))
    pattern.append(meanstdv(h24))
    pattern.append(meanstdv(h24r))
        
    if cluster!=0:
        color=getcolor(cluster,uniq_cluster,colorlist)
    else:
        color='white'
    for i in range(0,len(pattern)):
        draw.rectangle([barstart+25*i,y4,barstart+25*(i+1),y4-100*pattern[i][0]],fill=color,outline='black')
        draw.rectangle([barstart+10+25*i,y4-100*(pattern[i][0]+pattern[i][1]),barstart+15+25*i,y4-100*(pattern[i][0]-pattern[i][1])],fill='black',outline='black')
    

###################Function part end###############################
gene=sys.argv[2]
handle=open('subexon.txt')
array=extract(gene,handle)

Exon_pos=Exonorder(array)
Subexons=Exoncordinate(array) #stores genomic cordinates of all exons
total=sumlen(len(Subexons),Subexons)# It will be used to set the size of picture

startcor=Exonstart(array) #stores genomic cordinates of starting site of each unique exon 
handle.close()

handle1=open('splicingvar.txt')
array1=extract(gene,handle1)   #handle1 should be splicing variants file

var=collectvariant(array1) #get the unique accession of all splicing variants                                 
Exon_pos1=Exon_no(var,array1) #
varcor=var_cor(var,array1)

handle1.close()

sample=sys.argv[1]
inputfilename=sample+'_mappingout.txt'
handle2=open(inputfilename) #open  mappingoutput
array2=extractmap(gene,handle2)   #read mappingoutput
handle2.close()


var_cluster=[]
quanti=[]

for i in range(0,len(array2)):
    var_cluster.append(int(array2[i][6])); #get a list of clusterred result 
    quanti.append(max(array2[i][7].split(',')))


start=400
height=50
uniq_cluster=list(set(var_cluster)) #get the unique cluster, used for color setting
ymax=int(float(max(quanti)))+1
print len(uniq_cluster),uniq_cluster,ymax
yaxis=100+len(uniq_cluster)*120*ymax
setwidth1=start+total+20*len(Subexons)
setheight=200+60*len(var)+60*len(uniq_cluster)+yaxis

uniqcount=[] #count the unique peptides for each cluster
for i in range(0,len(uniq_cluster)):
    clus=uniq_cluster[i]
    psm=[]
    for i in range(0,len(array2)):
        if int(array2[i][6])==clus:
            psm.append(array2[i][2].upper())

    countuniqpep=len(list(set(psm)))
    #print countuniqpep
    uniqcount.append(countuniqpep)

print uniqcount
maxuniq=max(uniqcount)
setwidth2=400+230*maxuniq
print setwidth1,setwidth2
setwidth=max(setwidth1,setwidth2)

im=Image.new('RGB',(setwidth,setheight),'#ddd')    
draw=ImageDraw.Draw(im)
fontPath='/usr/share/fonts/truetype/ttf-dejavu/DejaVuSansCondensed-Bold.ttf'
font=ImageFont.truetype(fontPath,30)
fontPath2='/usr/share/fonts/truetype/ttf-dejavu/DejaVuSans-Bold.ttf'
font2=ImageFont.truetype(fontPath2,16)


exoncor=[] #this list will store starting cordinates of exons 
exoncor1=[] #this list will store ending cordinates of exons
exonstart=[] #this list will store starting cordinates of unique exons
exonend=[]   #this list will store starting cordinates of unique exons


################draw image of subexons###########
draw.text((10,50),gene,font=font,fill='black')
draw.text((10,100),'Subexons: ',font=font,fill='black')

for i in range(0,len(Subexons)):
    stx=(Exon_pos[i]-1)*20+sumlen(i,Subexons)
    exoncor.append(stx)
    y=100
    width=max(Subexons[i])-min(Subexons[i])
    exoncor1.append(stx+width)
    y1=y+height/2-5
    draw.rectangle([start+stx,y,start+stx+width,y+height],fill='white',outline='black');
    if i!=len(Subexons)-1:
        draw.rectangle([start+stx+width,y1,start+stx+width+20,y1+10],fill='blue',outline='blue');


for i in range(1,max(Exon_pos)+1):
    index=Exon_pos.index(i)
    exonstart.append(exoncor[index])
    exonend.append(exoncor1[index]);

for i in range(0,len(exonstart)):
    x_cor=exonstart[i]+(exonend[i]-exonstart[i])/2
    y_cor=50
    draw.text((start+x_cor,y_cor),str(i+1),font=font,fill='blue');

############draw image of known splicing variants#####
for j in range(0,len(var)):
    y=160+60*j
    y1=y+height/2-5
    stx=0
    draw.text((10,y),var[j],font=font,fill='black')
    minnum=min(Exon_pos1[j])# the num of the first exon of this variant
    if array[0][3]=='+': #check seq direction is + or -
        stx1=start+min(min(varcor[j]))-startcor[minnum-1]+exonstart[minnum-1]
    else:
        stx1=start+startcor[minnum-1]-max(max(varcor[j]))+exonstart[minnum-1]
    for i in range(0,len(varcor[j])):
        cor=varcor[j][i]
        num=Exon_pos1[j][i]#corresponding exon of the cordinates
        width=max(varcor[j][i])-min(varcor[j][i])
        if array[0][3]=='+':
            stx=min(cor)-startcor[num-1]+exonstart[num-1]
            draw.rectangle([start+stx,y,start+stx+width,y+height],fill='green',outline='black');
        else:
            stx=startcor[num-1]-max(cor)+exonstart[num-1]
            draw.rectangle([start+stx,y,start+stx+width,y+height],fill='green',outline='black');
            #print startcor[num-1],max(cor),exonstart[num-1],stx
    if minnum!=num:        
        draw.rectangle([stx1,y1,start+stx,y1+10],fill='green',outline='green')
    else:
        if array[0][3]=='-':
            stx=startcor[num-1]-max(min(varcor[j]))+exonstart[num-1]
            draw.rectangle([stx1,y1,start+stx,y1+10],fill='green',outline='green')
        else:
            stx=min(max(varcor[j]))-startcor[num-1]+exonstart[num-1]
            draw.rectangle([stx1,y1,start+stx,y1+10],fill='green',outline='green')
    
######Map clustered Peptide to known splicing variants##############
colorlist=['blue','lime','red','magenta','yellow','cyan',
           'purple','navy','orange','maroon','brown','teal','violet']


###draw mapped peptide######

for k in range(0,len(var)):
    y=160+60*k
    for i in range(0,len(array2)):
        color='violet'
        MTV=array2[i][11].split(',')
        if array2[i][11]=='':
            continue;
        tuple1=(int(array2[i][12]),int(array2[i][13]))
        tuple2=(min(min(varcor[k])),max(max(varcor[k])))
        exon1=int(array2[i][14])
        exon2=int(array2[i][15])
        if var[k] in MTV and overlap(tuple1,tuple2): #double check
            #print array2[i][4],var[k],MTV
            stx1=abs(int(array2[i][12])-startcor[exon1-1])+exonstart[exon1-1]
            stx2=abs(int(array2[i][13])-startcor[exon2-1])+exonstart[exon2-1]
            if exon1==exon2:
                draw.rectangle([start+stx1,y,start+stx2,y+height],fill=color,outline=color)
            else:
                
                draw.rectangle([start+stx1,y,start+exonstart[exon2-1]-20,y+height],fill=color,outline=color)
                draw.rectangle([start+exonstart[exon2-1],y,start+stx2,y+height],fill=color,outline=color)
                    
                                            
bottom=y+60

for i in range(0,len(uniq_cluster)):
    y1=bottom+60*i
    count=var_cluster.count(uniq_cluster[i])#count how many psm in each cluster

    draw.text((150,y1),'%s(%s)'%('cluster'+str(uniq_cluster[i]),str(count)),font=font,fill='black')
    if uniq_cluster[i]!=0:
        draw.rectangle([150,y1,0,y1+50],fill=colorlist[i],outline=colorlist[i])
        
    else:
        draw.rectangle([150,y1,0,y1+50],fill='white',outline='white')


for i in range(0,len(array2)):  #draw peptides for each cluster
    if array2[i][11]=='':
        continue;
    if int(array2[i][6])!=0:
        color=getcolor(int(array2[i][6]),uniq_cluster,colorlist)
    else:
        color='white'
        
    exon1=int(array2[i][14])
    exon2=int(array2[i][15])
    stx1=abs(int(array2[i][12])-startcor[exon1-1])+exonstart[exon1-1]
    stx2=abs(int(array2[i][13])-startcor[exon2-1])+exonstart[exon2-1]
    k1=uniq_cluster.index(int(array2[i][6]))
    y2=bottom+60*k1
    draw.rectangle([start+stx1,y2,start+stx2,y2+10],fill=color,outline=color)
    pepNO=array2[i][4]
    draw.text((start+stx1,y2+10),pepNO,font=font2,fill='black')
    if int(array2[i][10])==1 and int(array2[i][9])>1:
        k2=var.index(array2[i][11])
        Y=160+k2*60
        draw.text((10,Y),var[k2]+'(*)',font=font,fill='black')


############draw quantitative pattern of the peptides###############
    
for i in range(0,len(uniq_cluster)):
    cluster=uniq_cluster[i]
    axis_start=y1+100+(100*ymax+20)*i
    y4=axis_start+100*ymax
    drawxy(ymax)
    draw.text((100,(y4+axis_start)/2),'cluster'+str(cluster),font=font,fill='black')
    cluster_psm={}
    for i in range(0,len(array2)):
        if int(array2[i][6])==cluster:
            pep=array2[i][2].upper()
            cluster_psm[pep]=int(array2[i][4])

    uniqpep=sorted(cluster_psm.iteritems(), key=operator.itemgetter(1))
    for i in range(0,len(uniqpep)):
        pep=uniqpep[i][0]
        barstart=400+i*230
        histogram(pep,cluster)
        pepNO=cluster_psm[pep]
        
        draw.text((barstart+50,y4),'peptide '+str(pepNO),font=font2,fill='black')
        
image2=gene+'_pattern_'+sample+'.jpg'
im.save(image2)
print image2+' saved'
