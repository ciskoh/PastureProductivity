#!/usr/bin/python
#Script to combine different sentinel Granules and transform them into an array for further analysis

#------------------------------SETTINGS-------------------------
# Prepare the environment
import sys
import os
import qgis
from qgis.core import *
from PyQt4.QtGui import *
from PyQt4.QtCore import *

QgsApplication.setPrefixPath("/usr", False)
app = QgsApplication([], False)
app.initQgis()
sys.path.append('/usr/share/qgis/python/plugins')
sys.path.append('/home/jkm2/GIT/QGIS-scripts/pyModule')

from processing.core.Processing import Processing
Processing.initialize()
import processing as p
import shutil
import sqlite3
from decimal import *
import numpy as np
import gdal
import subprocess as sb

#-----------------------------PARAMETERS-------------------------
#### ---folders
# folder containing preprocessed NDVI raster images
imgFold="/home/jkm2/GIS/Sentinel_preprocess/test/ndvi_output22"
# path to landscape map /study area map
lsc="/home/jkm2/GIS/bound/watershedUTM.tif"
# working directory for temporary and output folders
parent="/home/jkm2/GIS/Sentinel_preprocess/"
orWodName="database"

###--settings

#starting Date as YYYYMDD

startD=20160101

#-----------------------------FUNCTIONS--------------------------

def mkWod(parent, orWodName):
# 0.1 create or define working folder
#input: wod (STRING)     - path to parent directory
#       wodName (STRING) - name of folder
#output:global variables wod - main folder
#                        tempWod - working folder ~/workingFiles
#                        finWod  - output folder ~/output

# create folder names
    global wod
    global tempWod 
    global finWod
    wod=os.path.join(parent, orWodName)
    tempWod=os.path.join(wod,"workingFiles")
    finWod=os.path.join(wod, "output")
    # remove previous folders
    if not os.path.exists(wod):
        os.makedirs(wod)
    #TODO REMOVE comment to line below
    if os.path.exists(tempWod):
        shutil.rmtree(tempWod)

    #create new folders
    os.makedirs(tempWod)
    if not os.path.exists(finWod):
        os.makedirs(finWod)
    # check
    if not os.path.exists(wod):
         print "problem creating working directory"
    else:
        print "created working folder  at %s\n" %(wod)

def makeList(fold, ext):
# make dictionary of file names present in the folder
# input: fold - STRING - path to folder to scan for files
#        ext  - STRING - extendion of files to scan for e.g. ".tif"

# output: Nested list with dates and paths
    imgPath=[]
    dateList=[]
    for r,d,files in os.walk(fold):
        for f in files:
            if f.endswith(".tif"):
                imgD=f[11:19]
                if imgD > startD:
                    imgPath.append(os.path.join(r, f))
                    dateList.append(imgD)
    if len(imgPath)<1:
        print "couldn't find any image in folder"
        return []
    else:
        print "found %s  files in folder %s" %(len(imgPath), fold)
        return [dateList, imgPath]

def checkLst(orList, outList):
    # function to compare different lists and return only the ones that are not present in the second
    # input: orList  - LIST - list of dates to be processed
    #        outList - LIST - list of already processed images
    #comparing image lists
    RlistNP=[y for y in orList if y not in outList]
    #if sublist is of length 0 exit right now
    if len(RlistNP) == 0:
        print  "no image to correct.\n\n -->quitting script on trigger"
        raise SystemExit
    else:
        logstr2=  "\nfound %d unprocessed images" %(len(RlistNP))
        print logstr2
        print "\nthe following images will be processed"
        print "output of function checkList:\n\n"
        return RlistNP

def sameDate(imgList, date):
# group images with same date
# input: imgList - LIST - list of paths to images
#        date - STRING - date in YYYYMMDD format
# output: list of images with the same date 
    smallList=[]
    for inner in imgList:
        if os.path.basename(inner)[11:19] == str(date):
            smallList.append(inner)
    if len(smallList) < 4:
        print "for this date there are only %s images:\n" %(str(len(smallList)))
        print smallList
    return smallList

#TODO: discarded function?
def getExtOld(rast, format):
    # function to extract extension parameters of raster and export it in qgis and gdal format
    #input: rast - STRING - path to raster
    #        format - LOGICAL - format of extension (0: QGIS; 1: GDAL)
    #output: STRING with new extension
    rastEx=QgsRasterLayer(rast).extent()
    newExt=rastEx.toString().replace(" : ", ",")
    gdalExt=newExt.split(",")
    gdalExt=[gdalExt[i] for i in [0,3,1,2]]
    if format ==0:
        return newExt
    else:
        return " ".join(gdalExt)



def getExt(lay, ref, utmPix):
# 1.2 calculate extension of study sites
# input: lay - QGSRASTERLAYER - layer to use as reference
# ref - INTEGER- Destination UTM reference code ("EPSG") 
# utmPix - pixel resolution in meters

# output: QGSRASTERLAYER - reprojected and clipped layer
#         gdalExt - STRING (implicit) - new extension as global variable in gdal format
#         newExt - STRING (implicit) - new extension in Qgis format 
    if int(lay.crs().authid().split(":")[1]) != ref:
        print "the input layer is not in the correct CRS"
        raise SystemExit
    # size of layer
    xSize=lay.width()
    ySize=lay.height()
    # get extension of layer and clip input layer Xmin, Xmax, Ymin, Ymax

    xMin=round(Decimal(lay.extent().xMinimum()))
    yMin=round(Decimal(lay.extent().yMinimum()))
    # creating new extension based on minimum coordinates, pixel size and size of pixel
    newExt=[xMin,
            xMin+(xSize*utmPix), 
            yMin, 
            yMin+(ySize*utmPix)]

    global newExt
    print " new extension is %s\n" %(newExt)
    #resize input layer
    # TODO: remove as lsc might not be needed
#    newRast= os.path.join(finWod, "reprLay.tif")
#    #gdal command
#    cmd="gdalwarp -ot Byte -of GTiff -te %s %s %s %s -dstnodata -0 %s %s" %(newExt[0], newExt[2], newExt[1], newExt[3], lay.source(), newRast )
#    print cmd
#    os.system(cmd)

#    newLay=QgsRasterLayer(newRast)

#    if not newLay.isValid():
#       print "layer transformation not valid"
#        raise SystemExit
#    else:
#        return newLay


    # test crs of lay
def mkOvrOld(groupImg):
# create mosaic raster for each group of images
# input: imgList - LIST - list of paths to images to be grouped in one virtual layer
# output: path to virtual layer
    print "starting  creation of virtual raster"
    res=QgsRasterLayer(groupImg[0]).rasterUnitsPerPixelX()
    destRast=os.path.join(tempWod, "ND"+date+".tif")
    stList=" ".join(groupImg)
    ext=getExt(newLsc.source(), 1)
    cmd="gdalwarp -tr 10 10 -ot Float64 -tap %s %s" %( stList, destRast)
    print cmd
    #TODO remove if below
    if not QgsRasterLayer(destRast).isValid():
        os.system(cmd)
    if not QgsRasterLayer(destRast).isValid():
        print "problem creating virtual raster\n\n"
        raise SystemExit
    else:
        print "mosaic raster created"
        return destRast


def mkOvr(groupImg, ext):
# create mosaic raster for each group of images
# input: imgList - LIST - list of paths to images to be grouped in one virtual layer
#        ext - LIST - list of coordinates in xMin-xMax-yMin-yMax format
# output: path to mosaicked layer
    print " starting creation of virtual raster"
    res=QgsRasterLayer(groupImg[0]).rasterUnitsPerPixelX()
    destRast=os.path.join(tempWod, "ND"+date+".tif")
    stList=" ".join(groupImg)
    # excluded for now
    stExt="%s %s %s %s" %(ext[0], ext[3], ext[1], ext[2])
    cmd="gdal_merge.py -a_nodata -9999 -ot Float64 -o %s -of GTiff -ul_lr %s %s" %(destRast, stExt, stList)
    text_file= open(tempWod+"/cmd.txt", "w")
    text_file.write(cmd)
    text_file.close()
    print cmd
    global cmd
    os.system(cmd)
    
    # alternative using qgis plugin
   # stExt="%s,%s,%s,%s" %(ext[0], ext[1], ext[2], ext[3])
   # stList=",".join(groupImg)
   # p.runalg("grass7:i.image.mosaic", stList, stExt, 10, destRast)
   
   #check of results
    if not QgsRasterLayer(destRast).isValid():
        print "problem creating virtual raster\n\n"
        raise SystemExit
    else:
        print "virtual raster created"
        return destRast


#TODO discarded function?
def sameExt(tempWod,rast1,rast2):
# function to make rast1 same extent of rast2
#input: tempWod (STRING) - path to temporary working folder
#       rast1 (STRING) - path to raster that will change extension
#       rast2 (STRING) - path to raster that will serve as reference

#output: STRING - path to new raster 1 with modified extension

# New extent as Gdal formatted list (xmin, ymax, xmax, ymin)
    print "starting function sameExt"
    rast1ex=QgsRasterLayer(rast1).extent()
    rast2ex=QgsRasterLayer(rast2).extent()
    newExt=rast1ex.intersect(rast2ex).toString().replace(" : ", ",")
    newExtList=newExt.split(",")
    newExtList=[Decimal(i) for i in newExtList]
    gdalExt=[newExtList[i] for i in [0,3,1,2]]
    #getting size of raster in pixels
    xSize=round((newExtList[1]-newExtList[0])/Decimal(10))
    ySize=round((newExtList[3]-newExtList[2])/Decimal(10))
# setting new extension to raster 1
    newRast1 = os.path.join(tempWod, "newRast1.tif")
    # check if newRast1 is not already present and resized
    if not QgsRasterLayer(newRast1).isValid():
        if QgsRasterLayer(newRast1).extent() != newExt:
            cmd="gdal_translate -outsize %s %s -projwin %s %s %s %s %s %s " %(xSize, ySize, gdalExt[0], gdalExt[1], gdalExt[2], gdalExt[3],rast1, newRast1)
            print "resizing Raster 1 with command\n"+cmd
            os.system(cmd)
    # setting new extension to raster 2
    newRast2 = os.path.join(tempWod, "newRast2.tif")
    cmd="gdal_translate -outsize %s %s -projwin %s %s %s %s %s %s " %(xSize, ySize, gdalExt[0], gdalExt[1], gdalExt[2], gdalExt[3],rast2, newRast2)
    print cmd
    os.system(cmd)
    # check
    if not QgsRasterLayer(newRast1).isValid():
        print "problem resizing raster 1"
        raise SystemExit
    else:
        if not QgsRasterLayer(newRast2).isValid():
            print "problem resizing raster2"
            raise SystemExit
        else:
            print "completed translating raster to same extension\n"
            return [newRast1, newRast2]

def sameExt2(tempWod,rast1,rast2):
# function to make rast1 rast2 the same extent
#input: tempWod (STRING) - path to temporary working folder
#       rast1 (STRING) - path to raster that will change extension
#       rast2 (STRING) - path to raster that will serve as reference
#output gDal arrays?
    #get intersecting extension
    rast1ex=QgsRasterLayer(rast1).extent()
    rast2ex=QgsRasterLayer(rast2).extent()
    newExt=rast1ex.intersect(rast2ex).toString().replace(" : ", ",")
    newExtList=newExt.split(",")
    newExtList=[Decimal(i) for i in newExtList]
    gdalExt=[newExtList[i] for i in [0,3,1,2]]
    #getting parameters to read raster1
    
    xSize=round((newExtList[1]-newExtList[0])/Decimal(10))
    ySize=round((newExtList[3]-newExtList[2])/Decimal(10))
    ofX=round(abs(Decimal(rast1ex.xMinimum())-gdalExt[0])/Decimal(10))
    ofY=round(abs(Decimal(rast1ex.yMaximum())-gdalExt[1])/Decimal(10))
    print "\nRASTER 1\nxSize is %s,\nySize is %s,\nofX is %s,\nofY is %s" %(xSize, xSize, ofX, ofY)
    
    #read only relevant values from raster 1 as array
    gRast1=gdal.Open(rast1)
    band=gRast1.GetRasterBand(1)
    dRast1=gRast1.ReadAsArray(int(ofX),int(ofY),int(xSize),int(ySize))
    
    #getting parameters to read raster2
    xSize=round((newExtList[1]-newExtList[0])/Decimal(10))
    ySize=round((newExtList[3]-newExtList[2])/Decimal(10))
    ofX=round(abs(Decimal(rast2ex.xMinimum())-gdalExt[0])/Decimal(10))
    ofY=round(abs(Decimal(rast2ex.yMaximum())-gdalExt[1])/Decimal(10))
    print "\nRASTER2\nxSize is %s,\nySize is %s,\nofX is %s,\nofY is %s" %(int(xSize), int(ySize), int(ofX), int(ofY))
    
    #read only relevant values from raster 2 as array
    gRast2=gdal.Open(rast2)
    band=gRast2.GetRasterBand(1)
    dRast2=gRast2.ReadAsArray(int(ofX),int(ofY),int(xSize),int(ySize))
 
    #check arrays
    cShape=np.shape(dRast2)==np.shape(dRast1)
    cVals=len(np.unique(dRast1))>1 and len(np.unique(dRast2))>1
    if not cShape==True and cVals==True:
        print "problem extracting arrays :\ncShape is %s and cVals is %s" %(cShape, cVals)
        raise SystemExit
    else:
        print "Array extraction succesfull"
        return [dRast1,dRast2]

def imgMd(ipat, stArea):
#2 IMAGE METADATA AND LAYER
    # input: ipat (STRING)   - path to unprocessed image
    #      : stArea (VECTOR) - path to study area vector file
    # output: Dictionary with the following data:
    #{ key : value  : type }

    #  baseName : basename without extension : STRING
    #  orImg : Qgs layer of Image : LAYER
    #  extImg : Image extension : STRING
    #  imgCrs : reference system : CRSOBJECT
    global baseName
    global orImg
    global extImg
    global imgCrs
    baseName = os.path.basename(ipat)[0:-4]
    orImg = QgsRasterLayer(ipat, baseName)
    if not orImg.isValid():
        print "\nLayer failed to load!\n\n -->quitting script on image %s" % baseName
        raise SystemExit
    # image extension as string ADDED new extension to include only parts overlapping study site
    oldExt=orImg.extent()
    stExt=QgsVectorLayer(stArea, "stArea", "ogr").extent()
    newExt=oldExt.intersect(stExt)
    extImg=newExt.toString().replace(" : ",",")
    #check of reference system and reprojection
    imgCrs= orImg.crs()
    #output dictionary
    mdImg={}
    mdImg["baseName"]=baseName
    mdImg["orImg"]=orImg
    mdImg["extImg"]=extImg
    mdImg["imgCrs"]=imgCrs
    
    print "output of function 2 imgMd:\n\n %s" %str(mdImg)
    return mdImg

def sAr(finWod, date, myArray):
# check that folder and files exists and save array as .np file
# input: finWod -STRING - path to final folder
#        date   -INTEGER - date of image being processed
#        myArray-NUMPY ARRAY - array to be saved 
# output: NONE

# make directory to store numpy array (if needed)
    arFold=os.path.join(finWod, "database")
    if not os.path.exists(arFold):
        os.makedirs(arFold)
    # save array as date if not already present
    arPath=os.path.join(finWod, str(date)+".np")
    if not os.path.exists(arPath):
        np.save(arPath, myArray)
    #check
    if not np.array_equal(myArray, np.load(arPath)):
        print "problem saving array %s " %(date)
        #raise SystemExit
    else:
        print "saved array %s" %(date)

def main(date):
# fucntion to group all images froma certain date and clip them to the study area
# input: date - STRING - date to work on in YYYYMMDD format
#        lots of other implicit inputs
# output GTIff image ?

    print "\n\nworking on images  for date %s\n\n" %(str(date))
    # 2 merge images of the same date and clip to landscape 
    # 2.1 create groups of images based on newList

    # define groups based on dates 
    groupImg=sameDate(orImg[1], date)
    
    print "groups of images by date created\n"
    #2.2 combine raster images from the same date
    mosRast=mkOvr(groupImg,newExt)
    global mosRast
    print "mosaicked raster is %s " %(mosRast)
    #copy image to output folder
    shutil.copy2(mosRast, finWod)
    #TODO add transformation to Array

#ALTERNATIVE CODE using numpy arrays
   # lscToAr=gdal.Open(newLsc.source())
   # lscAr=np.array(lscToAr.GetRasterBand(1).ReadAsArray())
# 2 translate image to array
    #rast=gdal.Open(mosRast)
    #myArray=np.array(rast.GetRasterBand(1).ReadAsArray())
    #print "array shape is:image %s %s" %(myArray.shape[0], myArray.shape[1])


# 2.1 remove negative values from array
    #myArray[myArray < 0] = -100
# 2.2 remove values outside study area 
  #  for i in range(len(myArray)):
  #      if lscAr[i]!=1:
  #          myArray[i]==-100
# 3 save array as .np file

#save raster as array
    #arPath=os.path.join(finWod, date+".npy")
    #np.save(arPath, myArray)
    #return arPath
#-----------------------------CODE-------------------------------
# 0 Check for images to be processes
# 0.1 create working folders
mkWod(parent, orWodName)

# 0.2 get images List and establish date list
#TODO add check to work only on new image
orImg=makeList(imgFold, ".tif")

#check dates of processed images
procImg=makeList(finWod, ".tif")

#compare original images with processed images
try:
    uniqueDates=checkLst(orImg[0], procImg[0])
except:
    uniqueDates=list(set(orImg[0]))
    print "no processed images found"

print "unique dates are %s" %(uniqueDates)

#landscape preparation
#1.2 check compatibility of landscape image with sattelite image (crs)
newLsc=getExt(QgsRasterLayer(lsc), 32630, 10)
###main function to be applied for each image
# 1.1 start cycle through each image date available
for date in uniqueDates:
    main(date)

# closing cycle
print "all images processed and added to %s" %(finWod)
