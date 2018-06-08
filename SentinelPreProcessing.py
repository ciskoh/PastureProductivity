#!/usr/bin/env python

# -*- coding: utf-8 -*-

#----------------------------DESCRIPTION
#Script to preprocess sentinel images and obtain ONLY NDVI /SAVI / MSAVI values
#TODO:add function 9 to crop raster values to study area
#TODO 5.0 Altitude average uaing dem (as separate function)
#TODO add definition of output folder name within ndvi/savi/msavifunctions
#TODO: correct error below: "wrong parameter value: empty"

#------------------------------SETTINGS-------------------------
# Prepare the environment
import sys
import os
import qgis
from qgis.core import *
#from PyQt4.QtGui import *

QgsApplication.setPrefixPath("/usr", False)
app = QgsApplication([], False)
app.initQgis()
sys.path.append('/usr/share/qgis/python/plugins')
sys.path.append('/home/jkm2/GIT/QGIS-scripts/pyModule')

from processing.core.Processing import Processing
Processing.initialize()
import processing as p

from pyproj import Proj, transform
from multiprocessing import Pool
import datetime
import numpy as np
from osgeo import gdal
import shutil
import ogr2ogr

#----------------------------PARAMETERS
#input directory
# linux remote dir
ind="/mnt/cephfs/data/BFH/Geodata/World/Sentinel-2/S2MSI1C/GeoTIFF"
# linux local dir ind="/home/matt/Dropbox/ongoing/BFH-Pastures/gis data/Sentinel_preprocess/test/input"
#code of image tiles to be analysed
tiles=['T30STA', 'T30STB', 'T30SUA', 'T30SUB']

#years to be considered
yea=[2018]

#
#path to Dem for topographic correction
#remote linux
dem="/home/jkm2/GIS/DEM/ASTER_big/asterDemUTM/asterDemUTM_comp.tif"
#local linux dem="/home/matt/Dropbox/ongoing/BFH-Pastures/gis data/DEM/ASTER_big/asterDemUTM/asterDemUTM_comp.tif"

# output directory
#linux local dir outd="/home/matt/Dropbox/ongoing/BFH-Pastures/gis data/Sentinel_preprocess/test/output"
#linux remote dir
outd="/home/jkm2/GIS/Sentinel_preprocess/test/ndvi_output22"
# working folder for temporary folders
temp1="/home/jkm2/GIS/Sentinel_preprocess/test/WOD"

#path to cloud mask folder
clPath="/mnt/cephfs/data/BFH/Geodata/World/Sentinel-2/S2MSI1C/SAFE"

# path to study area vector file
stArea="/home/jkm2/GIS/bound/bb32630.shp"
#-------------------------------FUNCTIONS---------------------------------
#100.1 Stop and go function based on user input
#input lNum (INTEGER) - Line number
#      baseName (STRING) - name of image being processed
#output: either quit signal or nothing

def stopGo(lNum, basename):
    #prompt user signal
    cont= raw_input("script at line %s, \n CONTINUE? Y/n?" %lNum)
    if cont == "n":
        print "stop signal received from user for image %s at line %s" %(baseName, lNum)
        print "\n\n -->quitting script on user input"
        quit()

#### 1 TRIGGER

#1.1 function to search for images recursively
# Inputs:input iDir-folder empty list for root of images, tiles-code of tiles for area of interest)
# Output: double list with name and path of images in iDir

def srcImg(iDir, tiles):
    multifolds=(iDir +"/"+ tile for tile in tiles)
    bnlist=[]
    rlist=[]
    #search for tiff files in folder
    for i in multifolds:
        for r,d,files in os.walk(i):
            for f in files:
                if f.endswith(".tif"):
                    bnlist.append(f)
                    rlist.append(os.path.join(r,f))
    finlist=[bnlist, rlist]
    print "output of function 1.1 srcImg:\n\n %s"
    for i in finlist[0:50]:
        print i
    return finlist

#1.2 function to compare two lists and decide if the script continues
#input two double (full path and basename) image lists, unprocessed and processed:
#Orlist- Non processed, Outlist - Processed images
#Output: either a quit signal or a list of image to be processed

def checkLst(Orlist, Outlist):
    #comparing image lists
    RlistNP=[y for y in Orlist[1] if os.path.basename(y) not in Outlist[0]]
    #if sublist is of length 0 exit right now
    if len(RlistNP) == 0:
        print "no image to correct.\n\n -->quitting script on trigger"
        quit()
    else:
        logstr2= "\nfound %d unprocessed images" %(len(RlistNP))
        Lfile.write(logstr2)
        print logstr2
        print "\nthe following images will be processed"
        print "output of function checkList:\n\n"
        for i in RlistNP:
            print os.path.basename(i)
    return RlistNP

#2 IMAGE METADATA AND LAYER
# input: ipat (STRING)   - path to unprocessed image
#      : stArea (VECTOR) - path to study area vector file
# output: Dictionary with the following data:
#{ key : value  : type }

#  baseName : basename without extension : STRING
#  orImg : Qgs layer of Image : LAYER
#  extImg : Image extension : STRING
#  imgCrs : reference system : CRSOBJECT
#
def imgMd(ipat, stArea):
    global  baseName
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

# 2.1 create working directory
# input: baseName (STRING) : name of image to be processed
#        temp1 (STRING) : path to main working directory

# output : path to this image working directory (STRING)

def makeWod(baseName, temp1):
    directory=temp1+"/"+baseName+"/"
    if not os.path.exists(directory):
        os.makedirs(directory)
        print "output of function 2.1 makeWod, %s" %directory
    return directory

#3 DEM RESOLUTION CHECK AND TRANSFORM
#input: dem (STRING) - path to DEM,
#      imgCrs (CRSOBJECT) - Reference system of Image to be corrected,
#      wod (STRING) - path to working directory

#output: dem - path to new reprojected layer
def demCheck(dem, imgCrs, wod ):
    demCrs=QgsRasterLayer(dem).crs()
    print "\nImage CRS :%s\n dem crs:%s" %(imgCrs.description(), demCrs.description())

    if demCrs != imgCrs:
        print "\nreprojecting DEM to %s" %(imgCrs.description())
        ldem=QgsRasterLayer(dem)
        extdem="%f,%f,%f,%f" %(ldem.extent().xMinimum(),\
                ldem.extent().xMaximum(),\
                ldem.extent().yMinimum(),\
                ldem.extent().yMaximum())
        crsStr=imgCrs.authid()
        newdem=p.runalg("gdalogr:warpreproject",ldem,"",crsStr,"",30,0,False,extdem,"",5,4,75,6,1,False,0,False,"",None)
        LnDem=QgsRasterLayer(newdem['OUTPUT'])
        if LnDem.isValid():
            print "dem reprojected correctly"
            dem=newdem['OUTPUT']
            print "output of function 3 demCheck:\n\n %s" %dem
    else:
        print "\nall files are in the same CRS!"
    return dem

#4 BAND SEPARATION
# inputs: ipat (STRING) - path to multiband image to be processed,
#         bNum (INTEGER) - number of band to be extracted
#         baseName (STRING) : name of image to be processed
#         wod (STRING) - working directory to save geoTiff files

# output: path to single-band geotiff

def bndSep(ipat, bNum, baseName, wod):
    sBand={}
    #band name and path
    bname="band"+str(bNum)
    directory=wod
    path=directory+bname+".tif"
    # call to gdal translate with option
    gdal.Translate(path, ipat, bandList=[bNum])
    if not QgsRasterLayer(path).isValid():
        print "problem saving single band %s!\n\n -->quitting script on image %s, function bndSep" %(bNum,baseName)
        quit()
    print " output of function 4 bndSep:\n\n %s" %path
    return path

#TODO 5.0 Altitude average uaing dem (as separate function)
#    demStats=p.runalg("grass7:r.univar",dem,None,"","",False,extImg,None)
#    statFile=open(demStats['output'], "r")
#    line=statFile.read().splitlines()[1]

#### 5 PARAMETERS FOR S6 ATHMOSPHERIC CORRECTION
# input: orImg (QGSLAYER) - layer of multiband unprocessed image
#        baseName (STRING) - name of image to be processed
#        bNum (INTEGER) - number of band to be used
#        wod (STRING) - working directory to save geoTiff files
# output: (STRING) path to parameter file
def makePar(orImg, baseName, bNum, wod):

    # preparation parameters for 6S
    #path to parameter file
    band="band"+str(bNum)
    Spath=str(wod+band+"_.atcorParam.txt")
    #Writing paramater file
    Sfile=open(Spath, "w")

    #FIRST LINE : satellite type
    firstLine="25\n"
    Sfile.write(firstLine)
    #SECOND LINE
    #read end of basename for date and time
    dateS=baseName.split("_")[2]

    month=dateS[4:6]
    day=dateS[6:8]

    hours=dateS[9:11]
    minutes=dateS[11:13]
    # hours and day in decimal time
    dmin=round(100*float(minutes)/60,0)
    dtime=str(hours)+"."+str(dmin/100)[2:]

    #long lat (centre point) in WGS84
    lon=orImg.extent().xMinimum()+((orImg.extent().xMaximum()-orImg.extent().xMinimum())/2)
    lat=orImg.extent().yMinimum()+((orImg.extent().yMaximum()-orImg.extent().yMinimum())/2)
    # transform long and lat in WGS 84
    crsDest=QgsCoordinateReferenceSystem(4326)
    xform = QgsCoordinateTransform(orImg.crs(), crsDest)
    ncoord=xform.transform(QgsPoint(lon,lat))

    #print secondline
    secondLine=str(month)+" "+str(day)+" "+str(dtime)+" "+str(ncoord.x())+" "+str(ncoord.y())+"\n"
    Sfile.write(secondLine)

    #THIRD LINE of parameter file athmospheric model (1) continental
    if month < 3 or month > 10:
        atmMod=3
    else:
        atmMod=2
    thirdLine="%s \n" %atmMod
    #print third line
    Sfile.write(thirdLine)

    #FOURTH LINE: 1 continental model
    fourthLine="1\n"
    Sfile.write(fourthLine)

    #FIFTH LINE: Aerosol concentration model (visibility)
    fifthLine="-1\n"
    Sfile.write(fifthLine)

    #SIXTH LINE: Aerosol concentration model (visibility)
    sixthLine="0\n"
    Sfile.write(sixthLine)

    #SEVENTH LINE (average altitude in negative km)
#    missing: actual calculation of average altitude avgAlt
    avgAlt=1400
    avgDem=-1*float(avgAlt)/1000
    seventhLine=str(avgDem)+"\n"
    Sfile.write(seventhLine)
    # EIGHTH LINE: sensor altitude
    eighthLine="-1000\n"
    Sfile.write(eighthLine)

    # F : sensor band code
    bcodes=range(166, 179)
    ninthLine=str(bcodes[bNum-1])+"\n"
    Sfile.write(ninthLine)

    #output
    Sfile.close()
    print "output of function 5 makePar is:\n\n"
    print "1st line of 6S parameter file (SATELLITE TYPE) is %s\n" %firstLine
    print "2nd line of 6S parameter file is (month, day, hh.ddd, long, lat):\n %s" %secondLine
    print "3rd line of 6S parameter file is  (Athmospheric model):\n %s" %thirdLine
    print "4th line of 6S parameter file is (Aereosol model):\n %s" %fourthLine
    print "5th and 6th line of 6S parameter file is (Visibility model):\n %s and %s" %(fifthLine, sixthLine)
    print "7th and 8th lines of 6S parameter file (target and sensor altitude are:\n %s %s " %(seventhLine, eighthLine)
    print "9th line of 6S parameter file (band code) is %s " %(ninthLine)
    return Spath

#5.1 RANGE OF PIXEL VALUES: minimum and maximum pixel value
# input: singImg (STRING) - path to single band image
#        wod (STRING) - path to working directory

# output: STRING range of pixel values per band
def pixRange(singImg, wod):
    #range of values
    ds = gdal.Open(singImg)
    myarray = np.array(ds.GetRasterBand(1).ReadAsArray())
    pMin=np.nanmin(myarray)
    pMax=np.nanmax(myarray)
    pRange=str(pMin)+","+str(pMax)
    print "output of function 5.1 pixRange is:"
    print 'pixel range is %s' %(pRange)
    ds=None
    return pRange

#6 Athmospheric correction
#input: singImg (STRING) - path to single band image
#       baseName (STRING) - name of image being processed
#       par (STRING) - path to parameter file
#       pRange (STRING) - range of pixel values
#       bNum (INTEGER) - band number
#       extImg (STRING) - Image extension
#       wod (STRING) - path to working directory

#output: (STRING) path to single band corrected image
def atCorr(singImg, baseName, par, pRange, bNum, extImg, wod):
    # double check of input files
    if not QgsRasterLayer(singImg).isValid():
        print "problem with input image before atmospheric correction!\n\n -->quitting script on image %s" % baseName
        #stopGo(338, baseName)
    print "launching athmospheric correction on image "+baseName[-10:]
    #path to corrected image
    bname="band"+str(bNum)
    corPath=wod+"corr"+bname+".tif"
    corrImg=p.runalg("grass7:i.atcorr",singImg,pRange,None, None,par,pRange,False,True,False,False,extImg,0,corPath)
    if not QgsRasterLayer(corPath).isValid():
        print "problem correcting the image!\n\n -->quitting script on image %s" % baseName
        #stopGo(346, baseName)
    if corrImg['output'] != corPath:
        print "problem with corrected image location"
    else:
        print "output of function 6 atCorr is: \n\n"
        print corPath
    return corPath



#7 NDVI calculation
#input: Red (STRING) - path to Red band
#       Nir (STRING) - path to NIR band
#       wod (STRING) - path to working directory
#output: path to NDVI image
def ndCalc(Red, Nir, wod):
    fStr="(B-A)/(B+A)"
    ndPath=wod+"ndvi.tif"
    ndPath2=wod+"ndviClip.tif"
    nd=p.runalg("gdalogr:rastercalculator",Red,"1",Nir,"1",None,"1",None,"1",None,"1",None,"1",fStr,"",5,"",ndPath)
    # clip raster to study site area
#    prepare extension in gdal_translate format (xmin,ymin,xmax,ymax
    xmin,ymin,xmax,ymax=extImg.split(",")
    ext="%s %s %s %s" %(xmin,ymax,xmax,ymin)
    cmd="gdal_translate -ot Float32 -projwin %s %s %s" %(ext, ndPath, ndPath2)
    os.system(cmd)
    if not QgsRasterLayer(ndPath2).isValid():
        print "problem cerating NDVI!\n\n -->quitting script on image %s" % baseName
        raise SystemExit
    else:
        print "output of function 7 ndCalc is: \n\n"
        print ndPath2
        return ndPath2

#7.1 Alternative to NDVI: MSAVI calculation
#input: Red (STRING) - path to Red band
#       Nir (STRING) - path to NIR band
#       wod (STRING) - path to working directory
#output: path to MSAVI image
def saviCalc(Red, Nir, wod):
    fStr="((B-A)*(1+0. 7))/(B+A+0.7)"
    ndPath=wod+"SAVI.tif"
    nd=p.runalg("gdalogr:rastercalculator",Red,"1",Nir,"1",None,"1",None,"1",None,"1",None,"1",fStr,"",5,"",ndPath)
    if not QgsRasterLayer(ndPath).isValid():
        print "problem creating SAVI!\n\n -->quitting script on image %s" % baseName
        #stopGo(367, baseName)
    else:
        print "output of function 7 ndCalc is: \n\n"
        print ndPath 
    return ndPath

#7.2 Alternative to NDVI: MSAVI2 calculation
#input: Red (STRING) - path to Red band
#       Nir (STRING) - path to NIR band
#       wod (STRING) - path to working directory
#output: path to MSAVI image

def msavi2Calc(Red, Nir, wod):
    fStr="(2*B+1-sqrt (((2*B+1)^2)-(8*(B-A)))/2"
    ndPath=wod+"MSAVI2.tif"
    nd=p.runalg("gdalogr:rastercalculator",Red,"1",Nir,"1",None,"1",None,"1",None,"1",None,"1",fStr,"",5,"",ndPath)
    if not QgsRasterLayer(ndPath).isValid():
        print "problem creating MSAVI2!\n\n -->quitting script on image %s" % baseName
        raise SystemExit
    else:
        print "output of function 7 ndCalc is: \n\n"
        print ndPath
    return ndPath

#8 vector mask to exclude clouds
#input: img (STRING) - path to geotiff (NDVI in this case)
#       clPath (STRING) - path to cloud masks main location
#       baseName (STRING) - name of image being processed
#       wod (STRING) - path to working directory

#output (STRING) -path to masked ndvi (clouds excluded)
def clMask(img, clPath, baseName, wod):
    # find path to cloudmask
    year=baseName.split("_")[2][0:4]
    tile=baseName.split("_")[5]
    mainDir=os.path.join(clPath, tile, year, baseName+".SAFE")
    if not os.path.exists(mainDir):
        print "problem finding cloud mask folder"
        raise SystemExit
# look for cloudmask in remote folder
    for r, d, files in os.walk(mainDir):
        for f in files :
            if f.endswith("B00.gml"):
                clMpath=os.path.join(r,f)
    cpv=wod+"cloudVec.gml"
    shutil.copy2(clMpath, cpv)
# mask raster using vector
    cloudRes=QgsRasterLayer(img).rasterUnitsPerPixelX()
    cloudEx=extImg
    rMaskPath=wod+"clMaskNdvi.tif"
    
    print "starting cloud masking with the following parameters:"
    print "parameter 1: cloudMask vector: {}\n" .format(cpv)
    print "parameter 2: Raster image to mask:{} \n" .format(img)
    print "parameter 6: extension: {}\n" .format(cloudEx)
    print "parameter 7: resolution: {}\n" .format(cloudRes)
    print "parameter 10: output location: {} \n" .format(rMaskPath)
    x=datetime.datetime.now()
    #TODO: correct error below: "wrong parameter value: empty"
    p.runalg("grass7:r.mask.vect", cpv, img,"","",True, cloudEx, cloudRes, -1, 0.00001,rMaskPath)
    y=datetime.datetime.now()
    c=y-x
    d=divmod(c.days*86400+c.seconds, 60)
    print "finished cloud masking in %f minutes and %f seconds" %(d[0], d[1])
    if not QgsRasterLayer(rMaskPath).isValid():
        print "problem masking image %s function 8 clMask" %baseName
        raise SystemExit
    else:
        print"output of function 8 clMask is:\n %s" %rMaskPath
        return rMaskPath



#8.1 vector mask to exclude clouds using different procedure
#input: img (STRING) - path to geotiff (NDVI in this case)
#       clPath (STRING) - path to cloud masks main location
#       baseName (STRING) - name of image being processed
#       wod (STRING) - path to working directory

#output (STRING) -path to masked ndvi (clouds excluded)
def clMask2(img, clPath, baseName, wod):
    year=baseName. split("_")[2][0:4]
    tile=baseName.split("_")[5]
    mainDir=os.path.join(clPath, tile, year, baseName+".SAFE")
    if not os.path.exists(mainDir):
        print "problem finding cloud mask folder"
        raise SystemExit

# look for cloudmask in remote folder
    for r, d, files in os.walk(mainDir):
        for f in  files:
            if f.endswith("B00.gml"):
                clMpath=os.path.join(r,f)
    global clMpath
    cpv=wod+"cloudVec.shp"
    if not QgsVectorLayer(clMpath, "cl", "ogr").isValid():
        print "problem with original cloud mask or no clouds in the image"
        return img
    #translate gml to shapefile
    ogr2ogr.main(["","-f", "ESRI Shapefile", "-s_srs", "EPSG:32630", cpv, clMpath])
    if not QgsVectorLayer(cpv, "cl", "ogr").isValid():
        print "cloud vector translation didn't go well"
        raise SystemExit
    else:
        print "cloud vector translation did work!!"
    # rasterize cloud mask
    # parameters to rasterize
    fName=[field.name() for field in QgsVectorLayer(cpv, "cl", "ogr").pendingFields()][0]
    print "first field name is"
    print fName
    cloudRes=QgsRasterLayer(img).rasterUnitsPerPixelX()
    cloudRast=wod+"clMask.tif"
    #rasterization
    x=datetime.datetime.now()
    cmd="gdal_rasterize -burn 0 -a_nodata 1000 -a_srs %s -te %s -tr %s %s %s %s" %(imgCrs.authid(), extImg.replace(","," "), cloudRes, cloudRes, cpv, cloudRast)
    print cmd
    os.system(cmd)
    #p.runalg("gdalogr:rasterize", cpv, fName, 1, 10,10, extImg, False, 5,0,4,75, 6,1,False,0,"-burn 1 -a_srs 'EPSG:32630'",cloudRast)
    if not QgsRasterLayer(cloudRast).isValid():
        #turn null to 0 values in cloud mask and use them
        print "Cloud raster {cloudRast} has a problem \n"
        raise SystemExit
    # Mask out cloudy pixels
    rMaskPath=wod+"clMaskNDVI.tif"
    cmd="gdal_calc.py -A %s -B %s --outfile=%s --calc='A*(B>0)' " %(img, cloudRast, rMaskPath)
    print(cmd)
    os.system(cmd)
    #p.runalg("gdalogr:rastercalculator", img, "1", cloudRast, "1", None, "1", None, "1", None, "1", None, "1", "A*((-1*B)+1)", "", 5, "", rMaskPath)
    #TODO add gdal calc expression
    y=datetime.datetime.now()
    c=y-x
    d=divmod(c.days*86400+c.seconds, 60)
    print "finished cloud masking in %f minutes and %f seconds" %(d[0], d[1])
    if not QgsRasterLayer(rMaskPath).isValid():
        print "problem masking image %s function 8 clMask" %baseName
        raise systemExit
    else:
        print"output of function 8 clMask is:\n %s" %rMaskPath
        return rMaskPath
#9 crop to study area
#input: stArea

#100.2 Save output in appropriate folder
#input: Img (STRING) - Image to be saved
#       baseName (STRING) - name of image being processed
#       outd (STRING) - main output directory

#output: STRING path to exported image
def makeOut(Img, baseName, outd):
    print "starting makeOut"
    tile=baseName.split("_")[5]
    outDir=os.path.join(outd, tile)
    outPath=os.path.join(outDir, baseName+".tif")
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    shutil.copy2(Img, outPath)
    if not QgsRasterLayer(outPath).isValid():
        print "problem saving final image %s\n\n -->quitting script on function 100.2 makeOut"
        raise SystemExit
        
    print "output of function 100.2 makeOut is:\n %s" %outPath
    return outPath

#100.3 clean temporary files
#input : wod (STRING) - path to current image working directory
#        baseName(STRING) - name of current image
#output : none

def cleanUp(wod, baseName):
    shutil.rmtree(wod)
    print "removed working files for image %s" %baseName

# 100 MAIN FUNCTION to perfom athmospheric correction and NDVI calculation using all functions above
#input : ipat (STRING) - path to image to be processed
#        bnum ( LIST of integers) - band numbers to be processed
#        outd (STRING) - path to output folder

#output : string to final image (corrected and masked ndvi)
def main(ipat):
    print "\n\n***START PROCESSING IMAGE: %s***" %os.path.basename(ipat)
##   call function to get image metadata
    mdDic=imgMd(ipat, stArea)

    #output should be dictionary
    

    ## function to create working directory for this image
    wod=makeWod(baseName, temp1)
    global wod
    print "wod is %s" %wod
    #output is working folder path

    ## function to check for dem crs
    #parameters
        
    global dem
    dem=demCheck(dem, imgCrs, wod)
    #raise SystemExit
    #output is path to (new) dem

        ##------- RED IMAGE
    ## function to seprate multilayer image in single band (RED)
    bNum=4
    Red=bndSep(ipat, bNum,baseName, wod)
    #output is path to RED band geotiff

    ## function to create parameters for athm. corr. in single band image (RED)
    orImg=mdDic ["orImg"]
    redPar=makePar(orImg, baseName, bNum, wod)

    ## function to calculate pixel range in band (RED)
    redRange=pixRange(Red, wod)

    ## function to perform athm. corr. on single band image (RED)
    redCor=atCorr(Red, baseName, redPar, redRange, bNum, extImg, wod )

    # output should be path to corrected image (RED)
    print "corrected image is %s" %redCor

#    stopGo(478, baseName)

    ###------- NIR IMAGE
    ## function to seprate multilayer image in single band (NIR)
    bNum=8
    Nir=bndSep(ipat, bNum, baseName, wod)
    #output is path to NIR band geotiff

    ## function to create parameters for athm. corr. in single band image (NIR)
    nirPar=makePar(orImg, baseName, bNum, wod)

    ## function to calculate pixel range in band (NIR)
    nirRange=pixRange(Nir, wod)

    ## function to perform athm. corr. on single band image (NIR)
    nirCor=atCorr(Nir, baseName, nirPar, nirRange, bNum, extImg, wod )

    print "corrected image is %s" %nirCor
    # output should be path to corrected image (NIR)
#    stopGo(498, baseName)

    ## Function to calculate NDVI
    nd=ndCalc(redCor, nirCor, wod)

    print "ndvi should be in %s" %nd
    # output should be path to corrected image (NIR)
#    stopGo(504, baseName)

    ## function to mask ndvi using cloudmask of image

    msNdvi=clMask2(nd, clPath, baseName, wod)
#    raise SystemExi
    #stopGo(558, baseName)
    
    
    ##function to export final image (masked NDVI)
    fin=makeOut(msNdvi, baseName, outd)
    
    #if working files are not needed anymore
    
    #cleanUp(wod, baseName)

    #stopGo(514, baseName)

    return fin

#-------------------------------code---------------------------------

#
#####LOG
timestr=(
    '1{date:%Y-%m-%d %H:%M:%S}'.format( date=datetime.datetime.now() )
    )

logpath=outd+"/"+"log.txt"
Lfile=open(logpath, "w")
logtext= "\nAthmospheric correction using SentinelPreprocessing.py started at %s" %(timestr)
Lfile.write(logtext)
print logtext

#### 0- TRIGGER
#TRIGGER to check if there are new images

###Search in input directory

Orlist=srcImg(ind, tiles)
##Search in output directory

Outlist=srcImg(outd, tiles)
## get sublist of missing/non processed images
RlistNP=checkLst(Orlist, Outlist)

#sort list from most recent image
#procedural run of main function
for ipat in RlistNP:
    res=main(ipat)

#Parallel run of main function

##pool=Pool(44)
#res=pool.map(main,RlistNP)


##Log conclude
timestr2=(
    '1{date:%Y-%m-%d %H:%M:%S}'.format( date=datetime.datetime.now() )
    )

logstr="Corrected images are available in %s" %(outd)
logstr2="Script finished correctly at %s" %(timestr2)

Lfile.write(logstr)
Lfile.write(logstr2)
Lfile.close()
print logstr
print logstr2

# When your script is complete, call exitQgis() to remove the provider and
# layer registries from memory
app.exitQgis()

#
