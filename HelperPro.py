# ----------------------------------------------------------------------------------------
# HelperPro.py
# Version: ArcPro / Python 3+
# Creation Date: 2020-07-06
# Last Edit: 2020-11-06
# Creator:  Kirsten R. Hazler

# Summary:
# A library of generally useful helper functions. Many (most?) have been adapted from a previous script for Arc version 10.3/Python 2.x  

# ----------------------------------------------------------------------------------------

# Import modules
import os
import sys
import traceback
import time
import numpy
import math
from datetime import datetime as datetime 

try:
   arcpy
   print("arcpy is already loaded")
except:
   print("Initiating arcpy, which takes longer than it should...")
   import arcpy

from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")

# Set overwrite option so that existing data may be overwritten
arcpy.env.overwriteOutput = True

# Python version (for version 2/3 differences)
print('Using python interpreter version: ' + str(sys.version))
pyvers = sys.version_info.major


def getScratchMsg(scratchGDB):
   """Prints message informing user of where scratch output will be written"""
   if scratchGDB != "in_memory":
      msg = "Scratch outputs will be stored here: %s" % scratchGDB
   else:
      msg = "Scratch products are being stored in memory and will not persist. If processing fails inexplicably, or if you want to be able to inspect scratch products, try running this with a specified scratchGDB on disk."

   return msg


def printMsg(msg):
   arcpy.AddMessage(msg)
   print(msg)


def printWrng(msg):
   arcpy.AddWarning(msg)
   print('Warning: ' + msg)


def printErr(msg):
   arcpy.AddError(msg)
   print('Error: ' + msg)


def garbagePickup(trashList):
   """Deletes Arc files in list, with error handling. Argument must be a list."""
   for t in trashList:
      try:
         arcpy.Delete_management(t)
      except:
         pass
   return


def CleanFeatures(inFeats, outFeats):
   """Repairs geometry, then explodes multipart polygons to prepare features for geoprocessing."""

   # Process: Repair Geometry
   arcpy.RepairGeometry_management(inFeats, "DELETE_NULL")

   # Have to add the while/try/except below b/c polygon explosion sometimes fails inexplicably.
   # This gives it 10 tries to overcome the problem with repeated geometry repairs, then gives up.
   counter = 1
   while counter <= 10:
      try:
         # Process: Multipart To Singlepart
         arcpy.MultipartToSinglepart_management(inFeats, outFeats)

         counter = 11

      except:
         arcpy.AddMessage("Polygon explosion failed.")
         # Process: Repair Geometry
         arcpy.AddMessage("Trying to repair geometry (try # %s)" % str(counter))
         arcpy.RepairGeometry_management(inFeats, "DELETE_NULL")

         counter += 1

         if counter == 11:
            arcpy.AddMessage("Polygon explosion problem could not be resolved.  Copying features.")
            arcpy.CopyFeatures_management(inFeats, outFeats)

   return outFeats


def CleanClip(inFeats, clipFeats, outFeats, scratchGDB="in_memory"):
   """Clips the Input Features with the Clip Features.  The resulting features are then subjected to geometry repair
   and exploded (eliminating multipart polygons) """
   # # Determine where temporary data are written
   # msg = getScratchMsg(scratchGDB)
   # arcpy.AddMessage(msg)

   # Process: Clip
   tmpClip = scratchGDB + os.sep + "tmpClip"
   arcpy.Clip_analysis(inFeats, clipFeats, tmpClip)

   # Process: Clean Features
   CleanFeatures(tmpClip, outFeats)

   # Cleanup
   if scratchGDB == "in_memory":
      garbagePickup([tmpClip])

   return outFeats


def CleanErase(inFeats, eraseFeats, outFeats, scratchGDB="in_memory"):
   """Uses Eraser Features to erase portions of the Input Features, then repairs geometry and explodes any multipart
   polygons. """
   # # Determine where temporary data are written
   # msg = getScratchMsg(scratchGDB)
   # arcpy.AddMessage(msg)

   # Process: Erase
   tmpErased = scratchGDB + os.sep + "tmpErased"
   arcpy.Erase_analysis(inFeats, eraseFeats, tmpErased, "")

   # Process: Clean Features
   CleanFeatures(tmpErased, outFeats)

   # Cleanup
   if scratchGDB == "in_memory":
      garbagePickup([tmpErased])

   return outFeats


def countFeatures(features):
   """Gets count of features"""
   count = int((arcpy.GetCount_management(features)).getOutput(0))
   return count


def countSelectedFeatures(featureLyr):
   """Gets count of selected features in a feature layer"""
   desc = arcpy.Describe(featureLyr)
   count = len(desc.FIDSet)
   return count


def unique_values(table, field):
   """This function was obtained from:
   https://arcpy.wordpress.com/2012/02/01/create-a-list-of-unique-field-values/"""
   with arcpy.da.SearchCursor(table, [field]) as cursor:
      return sorted({row[0] for row in cursor})


def TabToDict(inTab, fldKey, fldValue):
   """Converts two fields in a table to a dictionary"""
   codeDict = {}
   with arcpy.da.SearchCursor(inTab, [fldKey, fldValue]) as sc:
      for row in sc:
         key = sc[0]
         val = sc[1]
         codeDict[key] = val
   return codeDict


def multiMeasure(meas, multi):
   """Given a measurement string such as "100 METERS" and a multiplier, multiplies the number by the specified
   multiplier, and returns a new measurement string along with its individual components """
   parseMeas = meas.split(" ")  # parse number and units
   num = float(parseMeas[0])  # convert string to number
   units = parseMeas[1]
   num = num * multi
   newMeas = str(num) + " " + units
   measTuple = (num, units, newMeas)
   return measTuple


def createTmpWorkspace():
   """Creates a new temporary geodatabase with a timestamp tag, within the current scratchFolder"""
   # Get time stamp
   ts = time.strftime("%Y%m%d_%H%M%S")  # timestamp

   # Create new file geodatabase
   gdbPath = arcpy.env.scratchFolder
   gdbName = 'tmp_%s.gdb' % ts
   tmpWorkspace = gdbPath + os.sep + gdbName
   arcpy.CreateFileGDB_management(gdbPath, gdbName)

   return tmpWorkspace


def tback():
   """Standard error handling routing to add to bottom of scripts"""
   tb = sys.exc_info()[2]
   tbinfo = traceback.format_tb(tb)[0]
   pymsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n " + str(sys.exc_info()[1])
   msgs = "ARCPY ERRORS:\n" + arcpy.GetMessages(2) + "\n"
   msgList = [pymsg, msgs]

   printErr(msgs)
   printErr(pymsg)
   printMsg(arcpy.GetMessages(1))

   return msgList


def clearSelection(fc):
   typeFC = (arcpy.Describe(fc)).dataType
   if typeFC == 'FeatureLayer':
      arcpy.SelectLayerByAttribute_management(fc, "CLEAR_SELECTION")


def Coalesce(inFeats, dilDist, outFeats, scratchGDB="in_memory"):
   """If a positive number is entered for the dilation distance, features are expanded outward by the specified
   distance, then shrunk back in by the same distance. This causes nearby features to coalesce. If a negative number
   is entered for the dilation distance, features are first shrunk, then expanded. This eliminates narrow portions of
   existing features, thereby simplifying them. It can also break narrow "bridges" between features that were
   formerly coalesced. """

   # If it's a string, parse dilation distance and get the negative
   if type(dilDist) == str:
      origDist, units, meas = multiMeasure(dilDist, 1)
      negDist, units, negMeas = multiMeasure(dilDist, -1)
   else:
      origDist = dilDist
      meas = dilDist
      negDist = -1 * origDist
      negMeas = negDist

   # Parameter check
   if origDist == 0:
      arcpy.AddError("You need to enter a non-zero value for the dilation distance")
      raise arcpy.ExecuteError

      # Set parameters. Dissolve parameter depends on dilation distance.
   if origDist > 0:
      dissolve1 = "ALL"
      dissolve2 = "NONE"
   else:
      dissolve1 = "NONE"
      dissolve2 = "ALL"

   # Process: Buffer
   Buff1 = scratchGDB + os.sep + "Buff1"
   arcpy.Buffer_analysis(inFeats, Buff1, meas, "FULL", "ROUND", dissolve1, "", "GEODESIC")

   # Process: Clean Features
   Clean_Buff1 = scratchGDB + os.sep + "CleanBuff1"
   CleanFeatures(Buff1, Clean_Buff1)

   # Process:  Generalize Features
   # This should prevent random processing failures on features with many vertices, and also speed processing in general
   # arcpy.Generalize_edit(Clean_Buff1, "0.1 Meters")

   # Eliminate gaps
   # Added step due to weird behavior on some buffers
   Clean_Buff1_ng = scratchGDB + os.sep + "Clean_Buff1_ng"
   arcpy.EliminatePolygonPart_management(Clean_Buff1, Clean_Buff1_ng, "AREA", "900 SQUAREMETERS", "", "CONTAINED_ONLY")

   # Process: Buffer
   Buff2 = scratchGDB + os.sep + "NegativeBuffer"
   arcpy.Buffer_analysis(Clean_Buff1_ng, Buff2, negMeas, "FULL", "ROUND", dissolve2, "", "GEODESIC")

   # Process: Clean Features to get final dilated features
   CleanFeatures(Buff2, outFeats)

   # Cleanup
   if scratchGDB == "in_memory":
      garbagePickup([Buff1, Clean_Buff1, Buff2])

   return outFeats


def ShrinkWrap(inFeats, dilDist, outFeats, smthMulti=8, scratchGDB="in_memory"):
   # Parse dilation distance, and increase it to get smoothing distance
   smthMulti = float(smthMulti)
   origDist, units, meas = multiMeasure(dilDist, 1)
   smthDist, units, smthMeas = multiMeasure(dilDist, smthMulti)

   # Parameter check
   if origDist <= 0:
      arcpy.AddError("You need to enter a positive, non-zero value for the dilation distance")
      raise arcpy.ExecuteError

      # tmpWorkspace = arcpy.env.scratchGDB
   # arcpy.AddMessage("Additional critical temporary products will be stored here: %s" % tmpWorkspace)

   # Set up empty trashList for later garbage collection
   trashList = []

   # Declare path/name of output data and workspace
   drive, path = os.path.splitdrive(outFeats)
   path, filename = os.path.split(path)
   myWorkspace = drive + path
   Output_fname = filename

   # Process:  Create Feature Class (to store output)
   arcpy.CreateFeatureclass_management(myWorkspace, Output_fname, "POLYGON", "", "", "", inFeats)

   # Process:  Clean Features
   # cleanFeats = tmpWorkspace + os.sep + "cleanFeats"
   cleanFeats = scratchGDB + os.sep + "cleanFeats"
   CleanFeatures(inFeats, cleanFeats)
   trashList.append(cleanFeats)

   # Process:  Dissolve Features
   # dissFeats = tmpWorkspace + os.sep + "dissFeats"
   # Writing to disk in hopes of stopping geoprocessing failure
   # arcpy.AddMessage("This feature class is stored here: %s" % dissFeats)
   dissFeats = scratchGDB + os.sep + "dissFeats"
   arcpy.Dissolve_management(cleanFeats, dissFeats, "", "", "SINGLE_PART", "")
   trashList.append(dissFeats)

   # Process:  Generalize Features
   # This should prevent random processing failures on features with many vertices, and also speed processing in general
   arcpy.Generalize_edit(dissFeats, "0.1 Meters")

   # Process:  Buffer Features
   # arcpy.AddMessage("Buffering features...")
   # buffFeats = tmpWorkspace + os.sep + "buffFeats"
   buffFeats = scratchGDB + os.sep + "buffFeats"
   arcpy.Buffer_analysis(dissFeats, buffFeats, meas, "", "", "ALL")
   trashList.append(buffFeats)

   # Process:  Explode Multiparts
   # explFeats = tmpWorkspace + os.sep + "explFeats"
   # Writing to disk in hopes of stopping geoprocessing failure
   # arcpy.AddMessage("This feature class is stored here: %s" % explFeats)
   explFeats = scratchGDB + os.sep + "explFeats"
   arcpy.MultipartToSinglepart_management(buffFeats, explFeats)
   trashList.append(explFeats)

   # Process:  Get Count
   numWraps = (arcpy.GetCount_management(explFeats)).getOutput(0)
   arcpy.AddMessage('Shrinkwrapping: There are %s features after consolidation' % numWraps)

   # Loop through the exploded buffer features
   counter = 1
   with arcpy.da.SearchCursor(explFeats, ["SHAPE@"]) as myFeats:
      for Feat in myFeats:
         arcpy.AddMessage('Working on shrink feature %s' % str(counter))
         featSHP = Feat[0]
         tmpFeat = scratchGDB + os.sep + "tmpFeat"
         arcpy.CopyFeatures_management(featSHP, tmpFeat)
         trashList.append(tmpFeat)

         # Process:  Repair Geometry
         arcpy.RepairGeometry_management(tmpFeat, "DELETE_NULL")

         # Process:  Make Feature Layer
         arcpy.MakeFeatureLayer_management(dissFeats, "dissFeatsLyr", "", "", "")
         trashList.append("dissFeatsLyr")

         # Process: Select Layer by Location (Get dissolved features within each exploded buffer feature)
         arcpy.SelectLayerByLocation_management("dissFeatsLyr", "INTERSECT", tmpFeat, "", "NEW_SELECTION")

         # Process:  Coalesce features (expand)
         coalFeats = scratchGDB + os.sep + 'coalFeats'
         Coalesce("dissFeatsLyr", smthMeas, coalFeats, scratchGDB)
         # Increasing the dilation distance improves smoothing and reduces the "dumbbell" effect. However, it can also cause some wonkiness which needs to be corrected in the next steps.
         trashList.append(coalFeats)

         # Merge coalesced feature with original features, and coalesce again.
         mergeFeats = scratchGDB + os.sep + 'mergeFeats'
         arcpy.Merge_management([coalFeats, "dissFeatsLyr"], mergeFeats, "")
         Coalesce(mergeFeats, "5 METERS", coalFeats, scratchGDB)

         # Eliminate gaps
         noGapFeats = scratchGDB + os.sep + "noGapFeats"
         arcpy.EliminatePolygonPart_management(coalFeats, noGapFeats, "PERCENT", "", 99, "CONTAINED_ONLY")

         # Process:  Append the final geometry to the ShrinkWrap feature class
         arcpy.AddMessage("Appending feature...")
         arcpy.Append_management(noGapFeats, outFeats, "NO_TEST", "", "")

         counter += 1
         del Feat

   # Cleanup
   if scratchGDB == "in_memory":
      garbagePickup(trashList)

   return outFeats


def DateStamp():
   return time.strftime('%Y%m%d')

def Hms():
   return time.strftime('%H:%M:%S')

def GetElapsedHours(t1, t2):
   """Gets the hours elapsed between the start time (t1) and the finish time (t2)."""
   hrs = round((t2 - t1) / 3600, 3)
   deltaString = str(hrs) + ' hours.'
   return deltaString
   
def GetElapsedTime (t0, t1):
   """Gets the time in days, hours, minutes, and seconds, elapsed between the start time (t0) and the finish time (t1).
   NOTE: This had to be modified from the function originally written for Python 2.x"""
   delta = t1 - t0
   (d, m, s) = (delta.days, int(delta.seconds/60), delta.seconds%60)
   h = int(m/60)
   deltaString = '%s days, %s hours, %s minutes, %s seconds' % (str(d), str(h), str(m), str(s))
   return deltaString
   
def createFGDB(FGDB):
   '''Checks to see if specified file geodatabase exists, and creates it if not.
   Parameters:
   - FGDB: full path to file geodatabase (e.g. r'C:\myDir\myGDB.gdb')
   '''
   gdbPath = os.path.dirname(FGDB)
   gdbName = os.path.basename(FGDB)
   
   if arcpy.Exists(FGDB):
      print("%s already exists." %gdbName)
      pass
   else:
      print("Creating new file geodatabase...")
      arcpy.CreateFileGDB_management(gdbPath, gdbName)
      print("%s created." %gdbName)
   return FGDB

def CompareSpatialRef(in_Data, in_Template):
   sr_In = arcpy.Describe(in_Data).spatialReference
   sr_Out = arcpy.Describe(in_Template).spatialReference
   srFacCode_In = sr_In.factoryCode
   # print("Input factory code: %s"%srFacCode_In)
   srFacCode_Out = sr_Out.factoryCode
   # print("Template factory code: %s"%srFacCode_Out)
   gcsFacCode_In = sr_In.GCS.factoryCode
   gcsFacCode_Out = sr_Out.GCS.factoryCode
    
   if srFacCode_In == srFacCode_Out:
      reproject = 0
      transform = 0
      geoTrans = ""
   else:
      reproject = 1
      
   if reproject == 1:
      if gcsFacCode_In == gcsFacCode_Out:
         transform = 0
         geoTrans = ""
      else:
         transList = arcpy.ListTransformations(sr_In, sr_Out)
         if len(transList) == 0:
            transform = 0
            geoTrans = ""
         else:
            transform = 1
            geoTrans = transList[0]
         
   return (sr_In, sr_Out, reproject, transform, geoTrans)
   
def ProjectToMatch_vec(in_Data, in_Template, out_Data, copy = 1):
   '''Check if input features and template data have same spatial reference.
   If so, make a copy. If not, reproject features to match template.
   
   Parameters:
   in_Data: input features to be reprojected or copied
   in_Template: dataset used to determine desired spatial reference
   out_Data: output features resulting from copy or reprojection
   copy: indicates whether to make a copy (1) or not (0) for data that don't need to be reprojected
   '''
   
   # Compare the spatial references of input and template data
   (sr_In, sr_Out, reproject, transform, geoTrans) = CompareSpatialRef(in_Data, in_Template)
   
   if reproject == 0:
      print('Coordinate systems for features and template data are the same.')
      if copy == 1: 
         print('Copying...')
         arcpy.CopyFeatures_management (in_Data, out_Data)
      else:
         print('Returning original data unchanged.')
         out_Data = in_Data
   else:
      print('Reprojecting features to match template...')
      if transform == 0:
         print('No geographic transformation needed...')
      else:
         print('Applying an appropriate geographic transformation...')
      arcpy.Project_management (in_Data, out_Data, sr_Out, geoTrans)
   return out_Data
   
def ProjectToMatch_ras(in_Data, in_Template, out_Data, resampleType = "NEAREST", cellSize = ""):
   '''Check if input raster and template raster have same spatial reference.
   If not, reproject input to match template.
   Parameters:
   in_Data = input raster to be reprojected
   in_Template = dataset used to determine desired spatial reference and cell alignment
   out_Data = output raster resulting from resampling
   resampleType = type of resampling to use (NEAREST, MAJORITY, BILINEAR, or CUBIC)
   '''
   
   # Compare the spatial references of input and template data
   (sr_In, sr_Out, reproject, transform, geoTrans) = CompareSpatialRef(in_Data, in_Template)
   
   if reproject == 0:
      print('Coordinate systems for input and template data are the same. No need to reproject.')
      return in_Data
   else:
      print('Reprojecting input raster to match template...')
      arcpy.env.snapRaster = in_Template
      if transform == 0:
         print('No geographic transformation needed...')
      else:
         print('Applying an appropriate geographic transformation...')
      arcpy.ProjectRaster_management (in_Data, out_Data, sr_Out, resampleType, cellSize, geoTrans)
      return out_Data
 
def Downscale_ras(in_Raster, in_Snap, out_Raster, resType = "BILINEAR", in_clpShp = "NONE"):
   '''Converts a lower resolution raster to one of higher resolution to match the cell size and alignment of the specified snap raster.
   
   Parameters:
   - in_Raster: input raster to be resampled. Enter NONE if using input points instead.
   - in_Snap: snap raster used to set output cell size and alignment; also acts as mask
   - out_Raster: output resampled raster
   - resType: raster resampling type (NEAREST, BILINEAR, CUBIC, or MAJORITY)
   - in_clpShp: input feature class used to clip the input raster. Enter NONE if no clipping is needed.
   '''
   
   # Set environment variables        
   arcpy.env.overwriteOutput = True
   arcpy.env.snapRaster = in_Snap
   arcpy.env.extent = in_Snap
   arcpy.env.mask = in_Snap
   cellSize = arcpy.GetRasterProperties_management(in_Snap, "CELLSIZEX").getOutput(0)
   scratchGDB = arcpy.env.scratchGDB
   
   if in_clpShp == "NONE":
      clpRast = in_Raster
   else:
      clpRast = scratchGDB + os.sep + "clpRast"
      print("Getting extents of clip shape...")
      desc = arcpy.Describe(in_clpShp)
      xmin = desc.extent.XMin
      xmax = desc.extent.XMax
      ymin = desc.extent.YMin
      ymax = desc.extent.YMax
      rect = "%s %s %s %s" %(xmin, ymin, xmax, ymax)
      print("Clipping raster...")
      arcpy.management.Clip(in_Raster, rect, clpRast, in_clpShp, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
   
   resRast = scratchGDB + os.sep + "resRast"
   tmpRast = ProjectToMatch_ras(clpRast, in_Snap, resRast, resType, cellSize)
   
   if tmpRast == clpRast:
      # If no re-projection occurred...
      print("Resampling...")
      arcpy.management.Resample(clpRast, resRast, cellSize, resType)
   else:
      pass

   print("Finalizing output and saving...")
   finRast = Con(in_Snap, resRast)
   finRast.save(out_Raster)
   return finRast

   print("Mission complete.")

def interpPoints(in_Points, valFld, in_Snap, out_Raster, in_clpShp = "NONE", interpType = "IDW", numPts = 9, maxDist = "", cellSize = ""):
   '''Converts a point dataset to a raster via specified interpolation method
   
   NOTES/LESSONS LEARNED: 
   - I tried multiple methods on the PMP point data. In the end I settled on the topographic method.
   - The topographic method is memory-intensive. It failed (memory allocation error) when I tried to interpolate the PMP data to 10-m. Instead, I ended up doing the topographic interpolation to a 250-m cell size, then resampling that output to get a 10-m resolution raster. 
   
   Parameters:
   - in_Points: input points with values to be interpolated
   - valFld: the field in the input points used to determine output raster values. 
   - in_Snap: snap raster used to set output cell size and alignment; also acts as mask
   - out_Raster: output resampled raster
   - in_clpShp: input feature class used to clip the input raster or input points. Enter NONE if no clipping is needed.
   - interpType: interpolation type (SPLINE, TREND2, TREND3, TOPO, or IDW). 
   - numPts: number of points used for interpolation. Ignored if TREND2 interpolation.
   - maxDist: maximum search radius for interpolation. Ignored if SPLINE or TREND2 interpolation.
   - cellSize: cell size of the output raster. If not specified, same as in_Snap raster.
   '''
   
   # timestamp
   t0 = datetime.now()
   
   # Set environment variables        
   arcpy.env.overwriteOutput = True
   arcpy.env.snapRaster = in_Snap
   arcpy.env.extent = in_Snap
   arcpy.env.mask = in_Snap
   scratchGDB = arcpy.env.scratchGDB
   
   if cellSize == "":
      cellSize = arcpy.GetRasterProperties_management(in_Snap, "CELLSIZEX").getOutput(0)
   
   if in_clpShp == "NONE":  
      clpPts = in_Points
   else:
      clpPts = scratchGDB + os.sep + "clpPts"
      print("Getting extents of clip shape...")
      desc = arcpy.Describe(in_clpShp)
      xmin = desc.extent.XMin
      xmax = desc.extent.XMax
      ymin = desc.extent.YMin
      ymax = desc.extent.YMax
      rect = "%s %s %s %s" %(xmin, ymin, xmax, ymax)
      print("Clipping feature class...")
      arcpy.Clip_analysis(in_Points, in_clpShp, clpPts)
   
   prjPts = scratchGDB + os.sep + "prjPts"
   tmpPts = ProjectToMatch_vec(clpPts, in_Snap, prjPts, copy = 1)
   
   if interpType == "IDW":
      # This is not as smooth as I'd like, but output makes more sense than spline for PMP points
      print("Interpolating points using inverse weighted squared distance...")
      radius = RadiusVariable(numPts, maxDist)
      finRast = Idw(tmpPts, valFld, cellSize, 2, radius)
   elif interpType == "SPLINE":
      # This did NOT work well with PMP data points - values WAY out of range on the periphery
      print("Interpolating points using spline...")
      finRast = Spline(tmpPts, valFld, cellSize, "REGULARIZED", "", numPts)
   elif interpType == "TREND2":
      print("Interpolating points using 2nd-order polynomial trend...")
      finRast = Trend(tmpPts, valFld, cellSize, 2, "LINEAR")
   elif interpType == "TREND3":
      print("Interpolating points using 3nd-order polynomial trend...")
      finRast = Trend(tmpPts, valFld, cellSize, 3, "LINEAR")
   elif interpType == "TOPO":
      print("Interpolating points using topographic algorithm...")
      ptFeats = TopoPointElevation([[tmpPts, valFld]])
      in_topo_features = [ptFeats]
      finRast = TopoToRaster(in_topo_features, cellSize, rect, "", "", "", "NO_ENFORCE", "SPOT")
   else:
      print("Interpolation type specification is invalid. Aborting.")
      sys.exit()
      
   print("Finalizing output and saving...")
   finRast.save(out_Raster)
   
   # timestamp
   t1 = datetime.now()
   ds = GetElapsedTime (t0, t1)
   print("Completed interpolation function. Time elapsed: %s" % ds)

def PolyToRaster(in_Poly, in_Fld, in_Snap, out_Rast):
   '''Converts polygons to raster based on specified field.
   
   Parameters:
   - in_Poly: input polygon feature class to be converted to raster
   - in_Fld: field in feature class used to determine raster values
   - in_Snap: input raster used to specify output coordinate system, processing extent, cell size, and alignment
   - out_Rast: output raster
   '''
   # Set overwrite to be true         
   arcpy.env.overwriteOutput = True
   
   # Specify scratch location
   scratchGDB = arcpy.env.scratchGDB
   
   # Set some output and environment variables
   out_Poly = scratchGDB + os.sep + "polyPrj"
   arcpy.env.snapRaster = in_Snap
   arcpy.env.extent = in_Snap
   arcpy.env.mask = in_Snap   # This shouldn't be necessary in this function.
   
   # Re-project polygons, if necessary
   out_Poly = ProjectToMatch_vec(in_Poly, in_Snap, out_Poly, copy = 0)
      
   # Convert to raster
   print("Rasterizing polygons...")
   arcpy.PolygonToRaster_conversion (out_Poly, in_Fld, out_Rast, "MAXIMUM_COMBINED_AREA", 'None', in_Snap)
   
   print("Rasterization complete.")
   return out_Rast


def JoinFast(ToTab, ToFld, FromTab, FromFld, JoinFlds):
   '''An alternative to arcpy's JoinField_management for table joins.
   Uses python dictionary and Search/Update cursors.
   Tested about 50x faster than arcpy.JoinFields_management for a 500k-row, one-to-one join.
   Adapted from: https://gis.stackexchange.com/questions/207943/speeding-up-join-in-arcpy

   Note: will remove existing fields from ToTab with the names of JoinFlds. This differs from JoinField_management.

   ToTab = The table to which fields will be added
   ToFld = The key field in ToTab, used to match records in FromTab
   FromTab = The table from which fields will be copied
   FromFld = the key field in FromTab, used to match records in ToTab
   JoinFlds = the list of fields to be added
   '''

   if type(JoinFlds) != list:
      JoinFlds = [JoinFlds]
   flds_info = [a for a in arcpy.ListFields(FromTab) if a.name in JoinFlds]
   if len(flds_info) == 0:
      print('No fields found, no changes made.')
      return ToTab
   else:
      flds = [f.name for f in flds_info]
      print('Joining [' + ', '.join(flds) + ']...')
   r = list(range(1, len(JoinFlds) + 1))
   joindict = {}
   with arcpy.da.SearchCursor(FromTab, [FromFld] + flds) as rows:
      for row in rows:
         joinval = row[0]
         joindict[joinval] = [row[a] for a in r]
   del row, rows
   tFlds = [a.name for a in arcpy.ListFields(ToTab)]
   # Add fields
   for j in JoinFlds:
      ft = [a.type for a in flds_info if a.name == j][0]
      if j in tFlds:
         arcpy.DeleteField_management(ToTab, j)
      if ft == 'String':
         arcpy.AddField_management(ToTab, j, "TEXT", field_length=8000)
      elif ft == 'Integer':
         arcpy.AddField_management(ToTab, j, "LONG")
      else:
         arcpy.AddField_management(ToTab, j, "DOUBLE")
   # Do updates
   with arcpy.da.UpdateCursor(ToTab, [ToFld] + flds) as recs:
      for rec in recs:
         keyval = rec[0]
         if keyval in joindict:
            for a in r:
               rec[a] = joindict[keyval][a - 1]
            recs.updateRow(rec)
   del rec, recs
   return ToTab


def fldQuantile(inTable, valFld, quantiles=10, byArea=False, qFld=None):

   if not qFld:
      qFld = valFld + '_quant' + str(quantiles)

   q = list(range(0, quantiles))
   vals = [a[0] for a in arcpy.da.SearchCursor(inTable, valFld) if a[0] is not None]

   if byArea:
      print('Calculating area-equalized classes (n=' + str(quantiles) + ') for ' + valFld + ' in field `' + qFld + '`...')
      arcpy.Sort_management(inTable, 'tmp_sort', [[valFld, "ASCENDING"]])
      area = sum([a[0] for a in arcpy.da.SearchCursor('tmp_sort', 'Shape_Area')]) / quantiles
      q_abrk = [(a + 1) * area for a in q]
      q_brk = [min(vals)]
      ra = 0
      cls1 = 0
      with arcpy.da.SearchCursor('tmp_sort', ["Shape_Area", valFld]) as sc:
         for i in sc:
            ra += i[0]
            cls = q[q_abrk.index(min([a for a in q_abrk if ra - i[0]/2 <= a]))]
            if cls1 != cls:
               # This is a break value.
               q_brk.append(i[1])
            cls1 = cls
   else:
      print('Calculating quantile classes (n=' + str(quantiles) + ') for ' + valFld + ' in field `' + qFld + '`...')
      q_brk = [numpy.quantile(vals, a / quantiles) for a in q]

   # Set up list with classes and value breaks
   q2 = list(set(q_brk))
   q2.sort()
   qf = [a + 1 for a in q]
   cls = [[a, qf[q_brk.index(a)]] for a in q2]

   if len(cls) != len(q):
      print('Not enough variation to produce ' + str(quantiles) + ' quantiles, some classes combined.')

   if qFld in [a.name for a in arcpy.ListFields(inTable)]:
      arcpy.DeleteField_management(inTable, qFld)
   arcpy.AddField_management(inTable, qFld, field_type="LONG")
   with arcpy.da.UpdateCursor(inTable, [valFld, qFld]) as curs:
      for r in curs:
         v = r[0]
         if v is not None:
            r[1] = max([c[1] for c in cls if v >= c[0]])
         else:
            r[1] = 0
         curs.updateRow(r)

   return inTable


def CopyFeats(inFeats, outFeats, dropFlds):
   '''
   :param inFeats: Input features
   :param outFeats: Output features
   :param dropFlds: Fields to exclude from output features
   :return: outFeats

   Used to copy a feature class with a subset of the input fields. Created to avoid using
   arcpy.DeleteField_management(), which can error and corrupt feature classes randomly.
   '''
   fields = [a.name for a in arcpy.ListFields(inFeats)]
   fieldinfo = arcpy.FieldInfo()
   f2 = [fieldinfo.addField(f, f, "HIDDEN", "") for f in fields if f in dropFlds]

   lyr = arcpy.MakeFeatureLayer_management(inFeats, field_info=fieldinfo)
   arcpy.CopyFeatures_management(lyr, outFeats)
   del lyr

   return outFeats


def nullToZero(table, field):
   if not isinstance(field, list):
      field = [field]
   for f in field:
      cb = '''def fn(x):
         if x is None:
            return 0
         else: 
            return x'''
      arcpy.CalculateField_management(table, f, 'fn(!' + f + '!)', code_block=cb)
   return table

# end
