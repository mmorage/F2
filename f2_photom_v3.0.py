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
#     with a sky radius of 40 pix and 10 pix width. This can be edited deppending of the crowding level of the image. 
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
#    V2.0 Search for the closest star from the input coordinate
#    V3.0 Ask the user to identify the Number of the detected star. The user should write this number in the terminal
#         Search for the closest star deleted
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
#    M.Mora 25.Aug.2014 mmora@astro.puc.cl
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
#RA1=float(sys.argv[4]) 
#RA2=float(sys.argv[5]) 
#RA3=float(sys.argv[6]) 
#DEC1=float(sys.argv[7]) 
#DEC2=float(sys.argv[8]) 
#DEC3=float(sys.argv[9]) 


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
        self.test=  iraf.hselect(images=imagen+'[0]',fields="$I,GAIN,RDNOISE,EXPTIME,AIRMASS,OBSID,FILTER",expr='yes')
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
OBSIDt=[]
FILTERt=[]

for line in lines:
    p=line.split()
    name.append(p[0])
    GAINt.append(float(p[1]))
    RDNOISEt.append(float(p[2]))
    EXPTIMEt.append(float(p[3]))
    AIRMASSt.append(float(p[4]))
    OBSIDt.append((p[5]))
    FILTERt.append((p[6]))



GAINarr=np.array(GAINt)
RDNOISEarr=np.array(RDNOISEt)
EXPTIMEarr=np.array(EXPTIMEt)
AIRMASSarr=np.array(AIRMASSt)
OBSIDarr=np.array(OBSIDt)
FILTERarr=np.array(FILTERt)


GAIN=GAINarr[0]
RDNOISE=RDNOISEarr[0]
EXPTIME=EXPTIMEarr[0]
AIRMASS=AIRMASSarr[0]
OBSID=OBSIDarr[0]
#OBJECT=OBJECTarr[0]
FILTER=FILTERarr[0]



print"\n"
print "GAIN=",GAIN
print "RDNOSIE=",RDNOISE
print "EXPTIME=",EXPTIME
print "AIRMASS=",AIRMASS
print "OBSID=",OBSID
#print "OBJECT=",OBJECT
print "FILTER=",FILTER


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
        self.test=  iraf.daophot.txdump(textfile=imagen+".phot.1",field="XC,YC,ID,MAG,MER",expr='yes')
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
iraf.tv.tvmark(frame=1,coords=ARCHIVO,mark='point',number='yes')






FILE=FILE=open(ARCHIVO,'r')
lines=FILE.readlines()
FILE.close()

#creating empty arrays
XTt=[]
YTt=[]
IDTt=[]


for line in lines:
    p=line.split()
    XTt.append(p[0])
    YTt.append(p[1])
    IDTt.append(int(p[2]))

Xcoo=np.array(XTt)
Ycoo=np.array(YTt)
ID=np.array(IDTt)

mmax=len(ID)-1

print "Detected stars are plotted in the display"
print "If the star ID is not visible, type EXIT"
numero2=raw_input ("Type the number corresponding to the STD star: ")

print ID[mmax]

found='no'
if (numero2=='EXIT'):
    sys.exit("STD not detected, program aborted")
else:
    numero=int(numero2)
    while (found =='no'):
        for m in range (len(ID)):
            print ID[m]
            if (numero==ID[m]):
#                print  Xcoo[m],Ycoo[m],ID[m]
                found='yes'
                fila_std=m
                print m
                m=m+1
            #if (numero != ID[m]):
            else:
                m=m+1                
            if int(numero) > ID[mmax]:
                #m=m+1
                print "\nERROR, Please open the file named:",ARCHIVO,"\n"
                print ARCHIVO," has a total of",mmax+1,"objetcs\n"
                sys.exit("NUMBER NOT FOUND, please re-run this task.")




#print fila_std
#import json

XCt=[]
YCt=[]
STAR_MARt=[]
STAR_ERRt=[]

#FILE_OUT=open(OBSID+FILTER,"w")
FP=FILE=open(ARCHIVO,'r')
for i, line in enumerate(FP):
    if i==fila_std:
        p2=line.replace("INDEF","99.999") #reemplazo INDEF POR 99.999
        p=p2.split()
        XCt.append(float(p[0]))
        YCt.append(float(p[1]))
        STAR_MARt.append(float(p[3]))
        STAR_ERRt.append(float(p[34]))
        Zp=25-(STAR_MARt[0]-STAR_MAG)
        Zperr=sqrt(STAR_ERRt[0]*STAR_ERRt[0] + STAR_MERR*STAR_MERR)
        print Zp,Zperr," "
        temp=sys.stdout                             # Store original stdout object for later
        sys.stdout =open(OBSID+FILTER+".txt",'w')       # redirect all output to the file
        print "#01 XC"
        print "#02 YC"
        print "#03 Zpoint=25-(Mag_photom - MAg Catalog)"
        print "#04 Error=SQRT(ERRO_MAG**2 + ERRO_Catalog**2)"
        print "#05 Photomertry Zp=25 3*FWHM"
        print "#06 ERROR Zp=25 3*FWHN"
        print "#07 XC remaining columns correspond to file NNNN.END.dum"
        print "#08 YC"
        print "#09 ID"
        print "#10 Photometry 3*FWHM = 05"
        print "#11 Potomerty 10 pix aperture"
        print "#12 ... etc
        print "#X  Y ZP ZPerr",FILTER,FILTER+"_err  Apertures from 10 to 39 and error"
        print p[0],p[1],Zp,Zperr,p[3],p[34],line
        sys.stdout.close()                          # closing the file
        sys.stdout=temp                             # back to the normal print command        
        print p[0],p[1],p[3],p[34],line
FP.close()

print "\n\nFINAL ARCHIVE:",OBSID+FILTER+".txt\n\n"



sys.exit("Success !!!")
