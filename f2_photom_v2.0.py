#!/usr/bin/env python


#################################################################
#
#  It does aperture photometry.
#  Use:
#  python f2_photom.py imagen.fits  Mag_star_from_catalog M_err_from_catalog
# 
#  DO NOT FORGET TO PROPER SETUP YOUR PYTHON PATH!!!
#
#  in a terminal: copy the output of which python to the first line of this program adding #!
#                 if the first line does not work
#
#  
#
#  Image must have been  reduced using f2pipeline i.e image.fits[sci]
#  
# 0. Software asumes that the star is centred within a region of 200x200 pix, therefore
#    it uses an area of 200 x 200 pix   
#    See area = "[924:1124,924:1124]"  (2048/2 +-100)
#
#   The script does the following:
#    
#
#
# 1. catch the important header keywords from the imagen.fits 
#    GAIN
#    RDNOISE
#    EXPTIME
#    AIRMASS
# 
#  2. Uses imstatistic to calculate the STDDEV of the Background: 
#     Iterates 10 times : low  = mode - 3 x STDEEV
#                         high = mode + 3 x STDEEV
#
#     It converges most of the time.
#
#  3. Use gemseeing to calculate the FWHM of the stars and to get the coordinates (which have an offset of 923 pixels) 
#
#  4. Run PHOT
#     Several apertures 
#     with a sky radius of 50 pix and 10 pix width. This can be edited deppending of the crowding level of the image. 
#   
#  
#  5. Run txdump
#
#  6. Creation of the final file where it is added 923+XC and 923+YC using awk
#     If you change the area, you need to change this as well !!!
#     PLEASE NOTE THAT is + 923 instead of + 924, CCD coords starts a 1 
#  
# 
#  7. DS9 must be open at the final stage
#
#
#    V.1.1 Changes done: Corrected the X and Y position
#                        Added display and tvmark of the objects with photometry. (point) 
#    V2.0
#
#
#
#
#
#    BUGGS=
#
# 
#
#
#
#
#
#
#    M.Mora 20.Aug.2014 mmora@astro.puc.cl
#
#
#
####################################################################




#############################
#
# importing systems and iraf
#
############################

import sys,os,string
from subprocess import call
import pyraf
from pyraf import iraf
from iraf import  gemini, gnirs, niri,f2,gemtools,images,tv
from iraf import imutil, imheader as imhead
from iraf import noao,digiphot,daophot
from iraf import beep
import math


# iraf language module cannot be imported, 
# so each task should be imported from iraf
#
from iraf import unlearn,lparam, system,sleep, concatenate as concat
from iraf import nsheaders, niri,niflat
from iraf import images,immatch, geomap,imutil

import math
import numpy as np
from numpy import *


image=sys.argv[1]
STAR_MAG=float(sys.argv[2])   #Measured magnitude from catalog
STAR_MERR=float(sys.argv[3])  #Measured error from catalog
RA1=float(sys.argv[4]) 
RA2=float(sys.argv[5]) 
RA3=float(sys.argv[6]) 
DEC1=float(sys.argv[7]) 
DEC2=float(sys.argv[8]) 
DEC3=float(sys.argv[9]) 


########
#
# area: central area of the image: 200 x 200 pix  EDIT THIS IF YOU WANT TO INCREASE THE AREA!!!
#
# change to the total area of the image [1:2048,1:2048]
#
#
#########

area="[924:1124,924:1124]"

################
#
# EDIT THESE PARAMETERS:
#
#
###################

ANNULUS="40"
DANNULUS="10"
NAPERTURES="10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39"
ZERO_POINT="25"




############################
#
# unlearn all parameters 
#
#########################


iraf.gemini.gnirs.nsheaders.instrument ='f2'
iraf.gemini.gnirs.nsheaders.instrument ='f2'
iraf.gemini.gnirs.nsheaders.instrument ="f2"
iraf.gemini.gnirs.nsheaders.instrument ="f2"
iraf.gemini.gnirs.nsheaders.logfile ="f2.logfile.log"
iraf.gemini.gnirs.nsheaders ('f2',logfile="f2.logfile.log")
iraf.gemini.gnirs.nsheaders ('f2',logfile="f2.logfile.log")
iraf.gemini.gnirs.nsheaders ('f2',logfile="f2.logfile.log")


print "\n Unlearning the following tasks:\n"

print "hselect\n" 
unlearn("hselect")
print "gemseeing\n"  
unlearn("gemseeing")

print "phot\n" 
unlearn("phot")
print "datapars"  
unlearn("datapars")
print "centerpars\n"
unlearn("centerpars")
print "fitskypars\n"
unlearn("fitskypars")
print "photpars\n"
unlearn("photpars")




#iraf.hselect(images=image+'[sci]',fields="$I,GAIN,RDNOISE",expr='yes')


#Write the hsel outpu to hsel.tmp
class out_to_a_file:                                                                               # iraf:  
    def __init__(self,imagen): 
        temp=sys.stdout                                                                            # Store original stdout object for later
        sys.stdout =open("hsel.tmp",'w')                                                           # redirect all output to the file
        self.test=  iraf.hselect(images=imagen+'[0]',fields="$I,GAIN,RDNOISE,EXPTIME,AIRMASS",expr='yes')
        sys.stdout.close()                          # closing the file
        sys.stdout=temp                             # back to the normal print command
        
out1=out_to_a_file(image)


#opening the file hsel.tmp

FILE=open('hsel.tmp','r')
lines=FILE.readlines()
FILE.close()

#creating empty arrays
name=[]
GAINt=[]
RDNOISEt=[]
EXPTIMEt=[]
AIRMASSt=[]
for line in lines:
    p=line.split()
    name.append(p[0])
    GAINt.append(float(p[1]))
    RDNOISEt.append(float(p[2]))
    EXPTIMEt.append(float(p[3]))
    AIRMASSt.append(float(p[4]))


GAINarr=np.array(GAINt)
RDNOISEarr=np.array(RDNOISEt)
EXPTIMEarr=np.array(EXPTIMEt)
AIRMASSarr=np.array(AIRMASSt)

GAIN=GAINarr[0]
RDNOISE=RDNOISEarr[0]
EXPTIME=EXPTIMEarr[0]
AIRMASS=AIRMASSarr[0]

print"\n"
print "GAIN=",GAIN
print "RDNOSIE=",RDNOISE
print "EXPTIME=",EXPTIME
print "AIRMASS=",AIRMASS

#iraf.imstatistics(images=image+"[sci]",fields="image,npix,mean,stddev,min,max",lsigma=1,usigma=1)




#############################################
#
#
# Calculationg the STDDEV of the background in a loop of 10 iterations 
#
#
##############################################

N=0
niter=10    # N of iterations, it can be changed
up='INDEF'
low='INDEF'

while (N < niter):


####################################
#
#   Put iraf output into a file
#       
#
#
####################################


    class iraf_to_a_file:                                # iraf:  sections "lista" > "output_file"
        def __init__(self,imagen,outfile):
            temp=sys.stdout                             # Store original stdout object for later
            sys.stdout =open(outfile,'w')               # redirect all output to the file
            self.test=iraf.imstatistics(images=imagen,fields="stddev,min,max,mode",lower=low,upper=up,lsigma=1,usigma=1)                     
            sys.stdout.close()                          # closing the file
            sys.stdout=temp                             # back to the normal print command

    #hsel_test=iraf.hselect(images=image+'[sci]',fields="$I,GAIN,RDNOISE",expr='yes')

    iraf_to_a_file(image+area+"[sci]",'imstat.tmp')

    #opening the file imstat.tmp

    FILE=open('imstat.tmp','r')
    lines=FILE.readlines()
    FILE.close()

    #creating empty arrays
    STD=[]
    MINt=[]
    MAXt=[]
    MODEt=[]

    for line in lines:
        if not line.startswith("#"):
            p=line.split()
            STD.append(float(p[0]))
            MINt.append(float(p[1]))
            MAXt.append(float(p[2]))
            MODEt.append(float(p[3]))
    STDarr=np.array(STD)
    MINarr=np.array(MINt)
    MAXarr=np.array(MAXt)
    MODEarr=np.array(MODEt)
    #Saving values for daofind

    STDDEV=STDarr[0]
    MIN=MINarr[0]
    MAX=MAXarr[0]
    MODE=MODEarr[0]
    print"\n"
    print "iter ",N
    print "STDDEV=",STDDEV
    print "MIN=",MIN
    print "MAX=",MAX
    print "MODE=",MODE
    up=MODE+3*STDDEV
    low=MODE-3*STDDEV
    
    N=N+1

#iraf.imstatistics(images=image+"[sci]",fields="stddev,min,max",lsigma=1,usigma=1)
#iraf.imstatistics(images=image+area+"[sci]",fields="stddev,min,max",lsigma=1,usigma=1)


#############
#
#
# gemseeing
#
#
#################


print "\n"
print "\n"
print "running gemseeing\n"


iraf.gemseeing.image=image
iraf.gemseeing.coords=image+".coo"
iraf.gemseeing.fl_update='no'
iraf.gemseeing.fl_keep='yes'
iraf.gemseeing.fl_overwrite='yes'
iraf.gemseeing.ron=RDNOISE
iraf.gemseeing.gain=GAIN
iraf.gemseeing.pixscale=0.18





#Write the hsel outpu to hsel.tmp
class out_to_a_file2:                               # iraf:  
    """write bla"""
    def __init__(self,imagen): 
        temp=sys.stdout                             # Store original stdout object for later
        sys.stdout =open("gemseeing.tmp",'w')       # redirect all output to the file
        self.test=  iraf.gemseeing(image=image+area)
        sys.stdout.close()                          # closing the file
        sys.stdout=temp                             # back to the normal print command
        
out1=out_to_a_file2(image)




os.system("tail -n 4 gemseeing.tmp  >delme.tmp")


FILE=open("delme.tmp",'r')
lines=FILE.readlines()
FILE.close()



#creating empty arrays
name=[]
Npsf=[]
FWHM=[]

for line in lines:
    if any(str(image) in line for image in lines):
        p=line.split()
#        print image,"here"
        #name.append(float(p[0]))
        Npsf.append(float(p[1]))
        FWHM.append(float(p[2]))
    break

FWHMpsf=FWHM[0]/0.18
#print FWHMpsf

SEEING=FWHM[0]

print "FINAL Value of STDDEV =", STDDEV
print "FINAL Value of FWHM =", FWHMpsf

print "MEASURED SEEING [arcsec]=",SEEING

unlearn("datapars") 

#######################################
#
#
#  Preparing PSET Parameters
#
#
########################################

iraf.daophot.datapars.fwhmpsf=FWHMpsf
iraf.daophot.datapars.sigma=STDDEV
iraf.daophot.datapars.readnoise=RDNOISE
iraf.daophot.datapars.epadu=GAIN
iraf.daophot.datapars.itime=EXPTIME
iraf.daophot.datapars.xairmass=AIRMASS



#########################################################
#
#
# APERTURES, ANNULUS, DANNULUS AND ZMAG
#
#
########################################################

#iraf.daophot.fitskypars.annulus=40
#iraf.daophot.fitskypars.dannulus=10
#iraf.daophot.photpars.apertures="10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39"
#iraf.daophot.photpars.zmag=25.0


########
#
# Aperture photometry 3*Seeing 
#
#######

PHOT_APER=FWHMpsf*3.

TOT_APER=str(PHOT_APER)+","+NAPERTURES
print TOT_APER

iraf.daophot.fitskypars.annulus=ANNULUS   
iraf.daophot.fitskypars.dannulus=DANNULUS  
#iraf.daophot.photpars.apertures=NAPERTURES
iraf.daophot.photpars.apertures=TOT_APER
iraf.daophot.photpars.zmag=ZERO_POINT



iraf.daophot.phot.coords=image+".coo"
iraf.daophot.phot.output=image+".phot.1"

iraf.phot(image=image+area+"[sci]")

###################
#
#
#  txdump, 
#
#
#################

print "running txdump",image+".phot.1", "XC,YC,MAG,MERR" 



#Write the hsel output of txdump
class out_to_a_file3:                               # iraf:  
    def __init__(self,imagen): 
        temp=sys.stdout                             # Store original stdout object for later
        sys.stdout =open(imagen+".phot.1.dum",'w')       # redirect all output to the file
        self.test=  iraf.daophot.txdump(textfile=imagen+".phot.1",field="XC,YC,MAG,MER",expr='yes')
        sys.stdout.close()                          # closing the file
        sys.stdout=temp                             # back to the normal print command
        
out3=out_to_a_file3(image)


###################
#
#
# Archive is in the coordinates of the small area, therefore it needs to add the central pixel coord: 1024 and 1024
#
#
###################


final=image+".phot.1.dum"
#print final
ARCHIVO=image+".phot.1.END.dum"


#################################
#
#
# The following 2 lines should be deleted if all image is being used
# in such a case the final archive is final i.e. image+.phot.1.dum
#
####################################


#[924:1124,924:1124]"


string="awk '{printf(\"%f %f \", $1+923, $2+923); for (i=3;i<=NF;i++) printf(\"%s \",$i); print \" \"}' "+final+" > "+ARCHIVO
os.system(string)



print "\n"
print "\n"
print "\n"
print "\n"
print "Apertures\n"
print  NAPERTURES
#"10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39\n"
print "Final Archive:",ARCHIVO


#Showing the image with the detected star 

unlearn("tvmark")

iraf.images.tv.display(image+'[sci]',1)
iraf.tv.tvmark(frame=1,coords=ARCHIVO,mark='point',number='no')


#########
#
# Reading the archive
#
#######

FILE=open(ARCHIVO,'r')
lines=FILE.readlines() 
FILE.close()
#creating empty arrays

XCt=[]
YCt=[]
STAR_MARt=[]
STAR_ERRt=[]
j=0
for line in lines:
    p2=line.replace("INDEF","99.999") #reemplazo INDEF POR 99.999
    p=p2.split()
    XCt.append(float(p[0]))
    YCt.append(float(p[1]))
    STAR_MARt.append(float(p[2]))
    STAR_ERRt.append(float(p[33]))
    j=j+1
XC=np.array(XCt)
YC=np.array(YCt)
STAR_MAR=np.array(STAR_MARt)
STAR_ERR=np.array(STAR_ERRt)

print XC[0],j,len(XC)

from iraf import imgtools, xy2rd

n=0
#sys.stdout =open(imagen+".phot.1.dum",'w')       # redirect all output to the file
temp=sys.stdout                             # Store original stdout object for later
sys.stdout =open("radec.tmp",'w')       # redirect all output to the file
while (n < len(XC)):
    Zp=25-(STAR_MAR[n]-STAR_MAG)
    Zperr=sqrt(STAR_ERR[n]*STAR_ERR[n] + STAR_MERR*STAR_MERR)
    print Zp,Zperr," ",
    iraf.xy2rd(infile=image+"[sci]",x=XC[n],y=YC[n],hms='yes')
    n=n+1
sys.stdout.close()                          # closing the file
sys.stdout=temp                             # back to the normal print command




FILE=open("radec.tmp",'r')
lines=FILE.readlines() 
FILE.close()
#creating empty arrays

Mt=[]
ERRt=[]
RAt=[]
DECt=[]

tmp=[]
ra_s=[]

dec_m=[]
dec_s=[]




#        sys.stdout =open(imagen+".phot.1.dum",'w')       # redirect all output to the file
#        self.test=  iraf.daophot.txdump(textfile=imagen+".phot.1",field="XC,YC,MAG,MER",expr='yes')
#        sys.stdout.close()                          # closing the file
#        sys.stdout=temp                             # back to the normal print command





temp=sys.stdout                             # Store original stdout object for later
sys.stdout=open("radec2.tmp",'w')

#creating empty arrays
import json

j=0
for line in lines:
    p2=line.replace(",","") #reemplazo la , de xy2rd por nada
    p=p2.split()
    #Mt.append(float(p[0]))
    #ERRt.append(float(p[1]))
    M=str(p[0])#.append(str(p[1]))
    Me=str(p[1])
    #ra_h[j],ra_m[j],ra_s[j]=p[4].split(":")
    ra_h=p[4].split(":")
    RAt.append(str(p[4]))
    dec_g=p[7].split(":")
    #dec_g[j],dec_m[j],dec_s[j]=p[7].split(":")
    DECt.append(str(p[7]))
    final=ra_h+dec_g
    final2=final.append(M)
    final3=final.append(Me)
    tmp=final
    #print final#ra_h,dec_g,M,Me#,ra_m,ra_s,dec_g,dec_m,dec_s 
    print tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7]
    j=j+1

sys.stdout.close()                          # closing the file
sys.stdout=temp                             # back to the normal print command



FILE4=open('radec2.tmp','r')
lines=FILE4.readlines()
FILE4.close()

#creating empty arrays
RA1p1=[]
RA2p1=[]
RA3p1=[]
DEC1p1=[]
DEC2p1=[]
DEC3p1=[]
Magp1=[]
Merrp1=[]

for line in lines:
    p=line.split()
    RA1p1.append(float(p[0]))
    RA2p1.append(float(p[1]))
    RA3p1.append(float(p[2]))
    DEC1p1.append(float(p[3]))
    DEC2p1.append(float(p[4]))
    DEC3p1.append(float(p[5]))
    Magp1.append(float(p[6]))
    Merrp1.append(float(p[7]))



RA1p=np.array(RA1p1) 
RA2p=np.array(RA2p1)
RA3p=np.array(RA3p1)
DEC1p=np.array(DEC1p1)
DEC2p=np.array(DEC2p1)
DEC3p=np.array(DEC3p1)
Magp=np.array(Magp1)
Merrp=np.array(Merrp1)

nn=len(RA1p)

DIST_RA=[]
DIST_DEC=[]
DIST_TOT=[]

DIST_RA=zeros(nn)
DIST_DEC=zeros(nn)
DIST_TOT=zeros(nn)

RA_DEG=[]
DEC_DEG=[]
RA_DEG=zeros(nn)
DEC_DEG=zeros(nn)

PI=math.pi
rad=PI/180.

print RA1,len(RA1p)

for l in range (len(RA1p)):
    #DIST_RA[l]= abs(RA1p[l]-RA1)/3600 + abs(RA2p[l]-RA2)/60 +abs(RA3p[l]-RA3)
    

    RA_DEG[l]= ((RA1p[l] + RA2p[l]*60. +RA3p[l]*3600))# -(RA1/3600. + RA2/60. +RA3))
    DEC_DEG[l]= DEC1p[l]*15 +DEC2p[l]/60*15+ DEC3p[l]/3600*15
                #DIST_DEC[l]=abs(DEC1p[l]-DEC1)/15 + abs(DEC2p[l]-DEC2)/60 +abs(DEC3p[l]-DEC3)
    #DIST_DEC[l]= ((abs(DEC1p[l]/15.) + DEC2p[l]/60. +DEC3p[l]) -(abs(DEC1/15.) + DEC2/60. +DEC3)) 
    RA_STR=RA1+RA2*60+RA3*3600
    DEC_STR=DEC1/15+DEC2/15*60+DEC3/15*3600
    
    DIST_TOT[l]= (math.sin((DEC_DEG[l])*rad)*math.sin((DEC_STR)*rad) + math.cos((DEC_DEG[l])*rad)*math.cos((DEC_STR)*rad)*math.cos((RA_DEG[l]-RA_STR)*rad))
    print DIST_TOT[l],RA1p[l],RA2p[l],RA3p[l],DEC1p[l],DEC2p[l],DEC3p[l],RA1,RA2,RA3,DEC1,DEC2,DEC3#,Magp[l],Merrp[l]

    #DIST_TOT[l]=sqrt(DIST_RA[l]*DIST_RA[l] + DIST_DEC[l]*DIST_DEC[l])

print min(DIST_TOT),argmin(DIST_TOT)    

f=argmin(DIST_TOT)    

print RA1p[f],RA2p[f],RA3p[f],DEC1p[f],DEC2p[f],DEC3p[f],RA1,RA2,RA3,DEC1,DEC2,DEC3,Magp[f],Merrp[f]




sys.exit("paro")
RA=np.array(RAt)
DEC=np.array(DECt)
M=np.array(Mt)
ERR=np.array(ERRt)
jmax=j
#print RA,DEC,jmax
k=0
kk=0
RA_SEC=[]
DEC_SEC=[]
while (k<jmax):
    Ra2=str(RA[k])#.split(':',1)); 
    Dec2=str(DEC[k])#.split(':',1));
   
    coordRA= Ra2.split(':',2)
    coordDEC=Dec2.split(':',2)
    New_ra=coordRA[2].replace(",","")
 #   print coordRA,coordDEC,New_ra


    RA_SEC=(float(New_ra[kk]))
    DEC_SEC=(float(coordDEC[kk]))
    k=k+1
    kk=kk+1

sys.exit("AQUI PARO")

#s=0
#while(s < jmax):
#    print New_ra[s],coordDEC[s]
#    s=s+1



#Write the hsel output of txdump
#class out_to_a_file3:                               # iraf:  
#    def __init__(self,imagen): 
#        temp=sys.stdout                             # Store original stdout object for later
#        sys.stdout =open(imagen+".phot.1.dum",'w')       # redirect all output to the file
#        self.test=  iraf.daophot.txdump(textfile=imagen+".phot.1",field="XC,YC,MAG,MER",expr='yes')
#        sys.stdout.close()                          # closing the file
#        sys.stdout=temp                             # back to the normal print command
#        
#out3=out_to_a_file3(image)




#str = "Line1-abcdef \nLine2-abc \nLine4-abcd";
#print str.split( );
#print str.split(' ', 1 );





#os.system("awk '{printf(\"%f %f \", $1+1024, $2+1024); for (i=3;i<=NF;i++) printf(\"%s \",$i); print \" \"}\'  'final' >delme.plis")
#awk '{printf("%f %f ", $1+1024, $2+1024); for (i=3;i<=NF;i++) printf("%s ",$i); print " "}' P9181_1.6_5x5sH.fits.phot.1.dum >lala
#os.sys( "awk '{print ($1+1024) " " ($2+1024) $3,$4,$5}' P9181_1.6_5x5sH.fits.phot.1.dum
#$6,$7,$8
#awk '{ for (x=3; x<63; x++)  print ($1+1024) " " ($2+1024)} $x' P9181_1.6_5x5sH.fits.phot.1.dum
#awk '{print $1+1024, $2+1024; for (i=3;i<=NF;i++) {print ${i}}'
# echo 20:31:20.0576, |awk '{split($0,a,":"); print a[1],a[2],a[3]}'
