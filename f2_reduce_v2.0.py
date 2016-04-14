#!/astro/iraf/UR/Ureka_osx-6_64_1.0/variants/internal/bin/python


##!/usr/bin/env python
#
# Test pipeline to reduce imaging F2 data on the fly.
# Based on the F2exmaples  procedure
# MM
#
# FIRST Edit the python PATH in the first line of this file :
# on a terminal type: which python and copy the path to python
# Do not forget to add #! at the beginning 
#
#
#
#
#  Use
#  python f2_reduce.py  obj.lis sky.lis darks.lis flats.lis flatdarks.lis shortdarks.lis  createcal(yes or no) interact (yes or no)
#  
#  E.g. First time running the script:
#
#  python f2_reduce.py  objects1K.txt  skyobject1K.txt darks_forK.lis flats.lis flatdarks.lis shortdarks.lis  yes
#
#  second time: Used the flats and shortflats already processed by the first running: 
#
#  python f2_reduce.py  objects2K.txt  skyobject2K.txt darks_forK.lis  flats.lis  flatdarks.lis shortdarks.lis  no no
#
#
# Output: A file like "GS-2014A-SV-NNN-NNN.fits" "GS-2014A-SV-NNN-NNN_sky.fits" and "GS-2014A-SV-NNN-NNN.txt  
#
#
#  
#
#
#  where: - obj.lis is the list of science images   
#         - sky.lis is the list of sky frames
#         - darks.lis is the list of darks frames
#         - flats.lis is the list of flats frames
#         - flatdarks.lis is  long exposure darks list. The name of the file MUST BE flatdarks.lis
#         - shortdarks.lis is the short exposure darks list. The name of the file MUST BE shortdarks.lis
#         - Yes 'yes' no capital letter. Used for the processing of the flats and darks if is typed no. No darks
#           and flats are processed. Useful for more than one set of data reduction. I.e. several sci images in one band.  
#           The first time runing this script, createcal must be yes, if not, it will not create calibration files such 
#           as flats and longdarks  
# 
#  IMPORTANT: 
#         - ds9 needs to be open  
#         - list files needs to have one f2 image without the .fits
#         - one f2 image per line
#         - last line MUST be an empty line:
#
#################### 
#f2image1
#f2image2
#f2iamge3
#
####################
# 
#
#  BUGS: Sometimes f2prepare fails after 48 files, uncomment the for loop of the line 230, 231 and 232  
#        I do not know why.  
#
#       The final alignements sometimes fails due to the P2 noise. In this case the script needs to be run
#       interactively i.e 
# 
#       python f2_reduce.py  objects2K.txt  skyobject2K.txt darks_forK.lis  flats.lis  flatdarks.lis shortdarks.lis  no yes
#
#
#
#
#
#
#
#
#
#
#
#
#
#







#############################
#
# importing systems and iraf
#
############################

import sys,os,string
from subprocess import call
import numpy as np
import pyraf
from pyraf import iraf
from iraf import  gemini, gnirs, niri,f2,gemtools,images,tv
from iraf import imutil, imheader as imhead
from iraf import beep

# iraf language module cannot be imported, 
# so each task should be imported from iraf
#
from iraf import unlearn,lparam, system,sleep, concatenate as concat
from iraf import nsheaders, niri,niflat
from iraf import images,immatch, geomap

##########################
#
# Reading the comamnd line 
#
#########################

obj=sys.argv[1]
sky=sys.argv[2]
darks=sys.argv[3]
flats=sys.argv[4]
flatdarks=sys.argv[5]
shortdarks=sys.argv[6]
createcal=sys.argv[7]
interact=sys.argv[8]




############################
#
# unlearn all parameters 
#
#########################


print "\n Unlearning the following tasks:\n"

print "gemini\n" 
unlearn ("gemini") 
print "gnirs\n" 
unlearn ("gnirs")
print "niri\n"  
unlearn ("niri")
print "f2\n" 
unlearn ("f2")
print "gemtools\n" 
unlearn ("gemtools")

print "geomap\n" 
unlearn ("geomap")




#deleting logfile if exist

iraf.delete("f2_image_logfile.log",verify='no')

# creating logfile

f2.logfile= "f2_image_logfile.log"


#loading headers keywords for F2

print "loading headers keywords for F2"


#
#Two times, just to be sure that it loads
# CRITICAL STEP!!!
#

iraf.gemini.gnirs.nsheaders.instrument ='f2'
iraf.gemini.gnirs.nsheaders.instrument ='f2'
iraf.gemini.gnirs.nsheaders.logfile ="f2.logfile.log"
iraf.gemini.gnirs.nsheaders ('f2',logfile="f2.logfile.log")


############# 
#
# concatenating files:  flats + dakrs + flatdarks +  shortdarks in calib.lis
# concatenating calib.lis + objs and sky in all.lis
#
##############
concat(input_files=flats+','+darks+','+flatdarks+','+shortdarks,output_file="calib.lis") 
concat(input_files="calib.lis"+','+obj+','+sky,output_file="all.lis") 



###########################
#
# Displaying all images.  ds9 Needs to be open!! 
#
###########################


FILE=open('all.lis','r')  # Opening the file, only read
lines=FILE.readlines()    # Reading lines
FILE.close()              # Closing the file


##comment this for loop if you do not need to display the images

#for line in lines:
#    seeim=line.strip()  # <-- cuts the "\n" at the end of each line in the arch#ive
#    images.tv.display (seeim+'.fits[1]', 1)
#    sleep(3)



 
#################################
#
# F2prepare the calibration data
#
#################################

print "Running feprepare..."

f2.f2prepare("@calib.lis",fl_vardq='yes',fl_correct='yes',fl_saturated='yes',fl_nonlinear='yes')

###############################
#
#
# Sometimes f2prepare has problems if the list is too long. In such a case uncomment:
#
#
##############################

#for line in lines:
#    seeim=line.strip()  # <-- cuts the "\n" at the end of each line in the arch#ive
#    f2.f2prepare(seeim,fl_vardq='yes',fl_correct='yes',fl_saturated='yes',fl_nonlinear='yes')



#################
#
# Normalized flats field and BPM
#
##################

#
# Values taken from f2examples.  I have no idea how it will affects
# others filters
#
iraf.niri.niflat.thresh_flo = 0.70
iraf.niri.niflat.thresh_fup = 1.20
iraf.niri.niflat.thresh_dlo = -50.
iraf.niri.niflat.thresh_dup = 600.

#############################
# Change this section if there is PWFS2 shadowing.
# avoid such area
#
# larger possible area "[300:1748,300:1748]" I.e OIWFS
##############################

iraf.niri.niflat.statsec = "[300:1300,300:1500]"


#
#Deleting old flat.fits and f2_bpm.pl and list files fflats.lis,fflatdarks.lis and fshortdarks.lis
#

if createcal=='yes':
    iraf.imdelete ("flat.fits,f2_bpm.pl", verify='no')
    iraf.delete ("fflats.lis,fflatdarks.lis,fshortdarks.lis", verify='no')
else:
    print "Using calibrations already created\n"
    print "Using flat.fits"


##
#
#  redirect the outpot of the screen to a file!
#
# PLEASE KEEP temp=sys.stdout  , and at the end of the screen-to-file, sys.stdout=temp  
# this small class does in iraf: sections "f@sky.lis" > "fsky.lis" 
#
##


class out_to_a_file:                                # iraf:  sections "lista" > "output_file"
    def __init__(self,lista,output_file):
        temp=sys.stdout                             # Store original stdout object for later
        sys.stdout =open(lista,'w')                 # redirect all output to the file
        self.test =  iraf.sections(output_file) 
        sys.stdout.close()                          # closing the file
        sys.stdout=temp                             # back to the normal print command
        


out1=out_to_a_file('fflats.lis',"f@flats.lis")
out2=out_to_a_file('fflatdarks.lis',"f@flatdarks.lis")
out3=out_to_a_file('fshortdarks.lis',"f@shortdarks.lis")


######################
#
## NO INTERACTIVO !!!
#
####################
if createcal=='yes':
    iraf.niri.niflat("@fflats.lis", flatfile="flat.fits", lampsoff="@fflatdarks.lis", \
                     darks="@fshortdarks.lis", bpmfile="f2_bpm.pl", fl_inter='no')
else:
    print "using flats.fits"
#####
#
#Displaying flat and bad pixel mask
#
#####

beep()
beep()
iraf.display("flat.fits[sci]",1)
sleep(5)                                 # 5 sec 

beep()
beep()
iraf.display("f2_bpm.pl",2)
sleep(5) 



########
#
# F2prepare  science
#
########

#copy name of archive containing the sci list to a generic name list obj.lis

from subprocess import call
call(["cp",obj,"obj.lis"])

#copy name of archive containing the sky list to a generic name list sky.lis

call(["cp",sky,"sky.lis"])


iraf.imdelete ("f@obj.lis", verify='no')
f2.f2prepare ("@obj.lis",  bpm="f2_bpm.pl", fl_vardq='yes', \
           fl_correct='yes', fl_saturated='yes', fl_nonlinear='yes')

iraf.imdelete ("f@sky.lis", verify='no')
f2.f2prepare ("@sky.lis", bpm="f2_bpm.pl", fl_vardq='yes', \
           fl_correct='yes', fl_saturated='yes', fl_nonlinear='yes')

###############
#
# Sky Frames
#
###############

#
# From F2 example:
#
# Construct a sky frame by identifying objects in each image, removing them, 
# and averaging the remaining good pixels. The object masks created by nisky 
# are saved so that they can be checked to ensure that the task is masking 
# objects in the images appropriately.

#####
#
# Create the dark frame for the sky images
#
######

#copy name of archive containing the dark list to a generic name list sky.lis
#copiar aca archivo de nombre_individual_darks a darks.lis

call(["cp",darks,"darks.lis"])


#interact=sys.argv[8]
if interact=='no':
    iraf.imdelete ("dark.fits", verify='no')
    iraf.delete ("fdarks.lis", verify='no')
    out4=out_to_a_file('fdarks.lis',"f@darks.lis ")
    #########
    #
    # gemcombine fdarks
    #
    ##########

    iraf.gemcombine ("@fdarks.lis", "dark.fits", combine="average", fl_vardq='yes', logfile='f2.logfile', fl_dqprop='yes')
    iraf.display("dark.fits[sci]",3)
    sleep(5)

    # Dark subtract the prepared sky images
    # nsheaders sets nireduce.statsec = "[300:1748,300:1748]"


    iraf.delete ("fsky.lis", verify='no')

    ##sections "f@sky.lis" > "fsky.lis" #

    out4=out_to_a_file('fsky.lis',"f@sky.lis ")

    iraf.imdelete ("df@sky.lis", verify='no')
    iraf.niri.nireduce ("@fsky.lis", outprefix="d", fl_sky='no', fl_autosky='no', \
               fl_scalesky='no', fl_dark='yes', darkimage="dark.fits", fl_flat='no')


    # Create the sky frame using the dark subtracted sky images
    # nsheaders sets nisky.statsec = "[300:1748,300:1748]"

    iraf.imdelete ("sky.fits", verify='no')
    iraf.delete ("dfskymsk.lis", verify='no')

    #sections "df@sky.lis//msk.pl" > "dfskymsk.lis"
    out5=out_to_a_file("dfskymsk.lis",'df@sky.lis//msk.pl')


    iraf.imdelete("@dfskymsk.lis", verify='no')
    iraf.delete ("dfsky.lis", verify='no')



    #sections "df@sky.lis" > "dfsky.lis"
    out6=out_to_a_file('dfsky.lis',"df@sky.lis")

    #agregado a mano el minmax
    #iraf.niri.nisky ("@dfsky.lis", outimage="sky.fits", combtype="median",rejtype="minmax", fl_keepmasks='yes')

    #original from f2example
    iraf.niri.nisky ("@dfsky.lis", outimage="sky.fits", combtype="median", fl_keepmasks='yes')

    from iraf import beep
    iraf.beep()
    iraf.beep()
    iraf.display("sky.fits[sci]",4)
    iraf.beep()
    iraf.beep() 
    sleep(5)

    #####################################################################################
    # Flat divide the sky frame
    #################################################################################

    iraf.imdelete ("fsky.fits", verify='no')
    iraf.nireduce ("sky.fits", outprefix="f", fl_sky='no', fl_autosky='no', \
                   fl_scalesky='no', fl_dark='no', fl_flat='yes', flatimage="flat.fits")


    iraf.beep()
    iraf.beep()
    iraf.display("fsky.fits[sci]",1) 
    sleep(5)
else:
    print"Running in interactive mode"
    print"Using prevoous files fsky and sky"

########
#
# SCIENCE DATA REDUCTION....
#
#######

# Dark subtract the prepared science images
# nsheaders sets nireduce.statsec = "[300:1748,300:1748]"

iraf.delete ("fobj.lis", verify='no')

#sections "f@obj.lis" > "fobj.lis"
out7=out_to_a_file('fobj.lis',"f@obj.lis")


iraf.imdelete ("df@obj.lis", verify='no')
iraf.nireduce ("@fobj.lis", outprefix="d", fl_sky='no', fl_autosky='no', \
          fl_scalesky='no', fl_dark='yes', darkimage="dark.fits", fl_flat='no')




iraf.delete ("dfobj.lis", verify='no')


#sections "df@obj.lis" > "dfobj.lis"
out8=out_to_a_file('dfobj.lis',"df@obj.lis")

iraf.imdelete ("fdf@obj.lis", verify='no')

# Flat divide the dark subtracted science images
iraf.nireduce ("@dfobj.lis", outprefix="f", fl_sky='no', fl_autosky='no', \
    fl_scalesky='no', fl_dark='no', fl_flat='yes', flatimage="flat.fits")

iraf.delete ("fdfobj.lis", verify='no')

#sections "fdf@obj.lis" > "fdfobj.lis"
out9=out_to_a_file('fdfobj.lis',"fdf@obj.lis")


iraf.imdelete ("rfdf@obj.lis", verify='no')

# Sky subtract the flat divided, dark subtracted science images
iraf.nireduce ("@fdfobj.lis", outprefix="r", fl_sky='yes', skyimage="fsky.fits", \
               fl_autosky='yes', fl_scalesky='yes', fl_dark='no', fl_flat='no')


########################
#
# Combining the science data
#
########################


iraf.imdelete ("obj_comb.fits", verify='no')
iraf.delete ("rfdfobj.lis", verify='no')


#sections "rfdf@obj.lis" > "rfdfobj.lis"
out10=out_to_a_file('rfdfobj.lis',"rfdf@obj.lis")


#####
#
# Del ejemplo de F2 example. Usa "[300:1748,300:1748]", 
#
#####
#iraf.imcoadd ("@rfdfobj.lis", outimage="obj_comb.fits", rotate='no', \
#    geofitgeom="shift", niter=1, statsec="[300:1748,300:1748]", \
#    badpix="f2_bpm.pl", fl_fixpix='yes', fl_find='yes', fl_map='yes', fl_trn='yes',\
#    fl_med='yes', fl_add='yes', fl_avg='yes', fl_scale='yes', fl_overwrite='yes',\
#    logfile=f2.logfile)


# Problemas con la imagen final de la astrometria. Evitando la zona de la sombra del PWFS2...
#"[300:1300,300:1500]"
# 
# Troubleshooting: If stars looks faint in the final images, please have a look at the 
# rfdfSYYYYNNNNNNNNN_trn.fits images. Check that stars are in the same position if not, then:
# 
#PWFS2 shadow confuses imcoadd. First approach to do it interactively.  
# added: datamin=0,datamax=20000,fl_inter='yes'.\
# 
#
#
#



##########################################
#imcoadd is fl_inter yes. 
#
#
#
##########################################
iraf.immatch.geomap.reject='2'
iraf.immatch.geomap.maxiter='4'

if interact=='yes':
    iraf.imcoadd ("@rfdfobj.lis", outimage="obj_comb.fits", rotate='no', \
                  datamin=0,datamax=20000,fl_inter='yes',\
                  geofitgeom="shift", niter=3, statsec="[1:1100,*]", \
                  badpix="f2_bpm.pl", fl_fixpix='yes', fl_find='yes', fl_map='yes', fl_trn='yes',\
                  fl_med='yes', fl_add='yes', fl_avg='yes', fl_scale='yes', fl_overwrite='yes',\
                  logfile=f2.logfile)
else:
    iraf.imcoadd ("@rfdfobj.lis", outimage="obj_comb.fits", rotate='no', \
                  datamin=0,datamax=20000,fl_inter='no',\
                  geofitgeom="shift", niter=3, statsec="[1:1100,*]", \
                  badpix="f2_bpm.pl", fl_fixpix='yes', fl_find='yes', fl_map='yes', fl_trn='yes',\
                  fl_med='yes', fl_add='yes', fl_avg='yes', fl_scale='yes', fl_overwrite='yes',\
                  logfile=f2.logfile)

####################
#
# Creating a standarized name:
# reading the headers
#
#
###################




class out_to_a_file:                                                                               # iraf:  
    def __init__(self,imagen): 
        temp=sys.stdout                                                                            # Store original stdout object for later
        sys.stdout =open("hsel.tmp",'w')                                                           # redirect all output to the file
        self.test=  iraf.hselect(images="obj_comb.fits[0]",fields="$I,OBSID,FILTER,GAIN,RDNOISE,EXPTIME,AIRMASS",expr='yes')
        sys.stdout.close()                          # closing the file
        sys.stdout=temp                             # back to the normal print command
        
out1=out_to_a_file("obj_comb.fits")


#opening the file hsel.tmp

FILE=open('hsel.tmp','r')
lines=FILE.readlines()
FILE.close()

#creating empty arrays
namet=[]
OBSIDt=[]
OBJECTt=[]
FILTERt=[]
GAINt=[]
RDNOISEt=[]
EXPTIMEt=[]
AIRMASSt=[]


for line in lines:
    p=line.split()
    namet.append(p[0])
    OBSIDt.append((p[1]))
#    OBJECTt.append((p[2]))
    FILTERt.append((p[2]))
    GAINt.append(float(p[3]))
    RDNOISEt.append(float(p[4]))
    EXPTIMEt.append(float(p[5]))
    AIRMASSt.append(float(p[6]))


namearr=np.array(namet)
OBSIDarr=np.array(OBSIDt)
#OBJECTarr=np.array(OBJECTt)
FILTERarr=np.array(FILTERt)
GAINarr=np.array(GAINt)
RDNOISEarr=np.array(RDNOISEt)
EXPTIMEarr=np.array(EXPTIMEt)
AIRMASSarr=np.array(AIRMASSt)

name=namearr[0]
OBSID=OBSIDarr[0]
#OBJECT=OBJECTarr[0]
FILTER=FILTERarr[0]
GAIN=GAINarr[0]
RDNOISE=RDNOISEarr[0]
EXPTIME=EXPTIMEarr[0]
AIRMASS=AIRMASSarr[0]

print"\n"
print "name=",name
print "OBSID=",OBSID
#print "OBJECT=",OBJECT
print "FILTER=",FILTER
print "GAIN=",GAIN
print "RDNOSIE=",RDNOISE
print "EXPTIME=",EXPTIME
print "AIRMASS=",AIRMASS


call(["cp","obj_comb.fits",OBSID+".fits"])
call(["cp","sky.fits",OBSID+"_sky.fits"])


temp=sys.stdout                             # Store original stdout object for later
sys.stdout=open(OBSID+".txt",'w')
print "File Name   OBSID  FILTER GAIN RDNOISE EXPTIME AIRMASS"
print name,OBSID,FILTER,GAIN,RDNOISE,EXPTIME,AIRMASS
sys.stdout.close()                          # closing the file
sys.stdout=temp                             # back to the normal print command


####################
#
# Clearning files 
#
###################
 

#iraf.delete(files="fS*.fits,df*,rdf*,rfdfS*,*.tmp",verify='no')

 
