# ---------------------------------------------------------------------------
# WtrshdImpact_Functions.py
# Version: ArcPro / Python 3+
# Creation Date: 2020-07-06
# Last Edit: 2022-05-24
# Creator: Kirsten R. Hazler
#
# Summary: Functions to produce the ConservationVision Watershed Impact Model.

# For background and references see:
# - Coastal Services Center. 2014. “Technical Guide for OpenNSPECT, Version 1.2.” National Oceanic and Atmospheric Administration (NOAA). https://coast.noaa.gov/digitalcoast/tools/opennspect
# - Renard, K.G., G.R. Foster, G.A. Weesies, D.K. McCool, and D.C. Yoder. 1997. “Predicting Soil Erosion by Water: A Guide to Conservation Planning with the Revised Universal Soil Loss Equation (RUSLE).” Agriculture Handbook 703. U.S. Department of Agriculture. 
# - Cronshey, R., R.H. McCuen, N. Miller, W. Rawls, S. Robbins, and D. Woodward. 1986. “Urban Hydrology for Small Watersheds (2nd Ed.).” Technical Release 55. Natural Resources Conservation Service, U.S. Department of Agriculture.

# NOTE: To set up inputs/outputs and string a series of functions together, use a separate workflow script that imports the functions from this one.
# ---------------------------------------------------------------------------

# Import modules
import HelperPro
from HelperPro import *

### General Functions
def getTruncVals(in_Raster, in_Mask = "NONE", numSD = 3):
   '''Based on raster statistics, calculates lower and upper cutoff values to be used for rescaling purposes.
   
   Parameters:
   in_Raster: The input raster to be analyzed
   in_Mask: Mask raster used to define the portion of the input raster to be analyzed
   numSD: The number of standard deviations used to determine the cutoff values
   '''
   
   try:
      r = Con(Raster(in_Mask), Raster(in_Raster))
   except:
      r = Raster(in_Raster)

   rMin = r.minimum
   rMax = r.maximum
   rMean = r.mean
   rSD = r.standardDeviation
   TruncMin = max(rMin,(rMean - numSD*rSD))
   TruncMax = min(rMax,(rMean + numSD*rSD))
   
   return (TruncMin, TruncMax)


### Functions for preparing gSSURGO data
def HydroGrp_vec(in_GDB):
   '''To the MUPOLYGON feature class, adds a field called "HYDROLGRP_DCD", containing the Hydrologic Soil Groups extracted from gSSURGO. Values range from A to D (with some compound classes possible, e.g., A/D). Also adds a field called "HydroGrpNum", which contains a numeric, simplified version of the hydrologic groups in which there are no compound groups and no nulls.
  
   This function modifies the input geodatabase by adding new tables and fields. It does not modify existing fields.
   
   Parameters:
   - in_GDB: input gSSURGO geodatabase
   
   Per OpenNSPECT guidance, numeric values for the groups are assigned as follows:
   - A = 1
   - B = 2
   - C = 3
   - D = 4
   
   Compound values (e.g., A/D) are assigned the latter group. Null values are assumed to be group D. 
   
   IMPORTANT: Prior to running this function, new data must be created within the input geodatabase, using a tool in the Soil Data Development Toolbox. Tools in this toolbox can only be run from within ArcMap (not ArcPro) and I haven't figured out a way to automate this with a script, so you need to do it manually.
   
   The Soil Data Development Toolbox must be added to ArcToolbox in ArcMap, and is available from https://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/home/?cid=nrcs142p2_053628#tools.
   
   TO CREATE THE DATA NEEDED FOR THIS FUNCTION:
   - Within ArcMap, add the MUPOLYGON feature class as a layer. 
     - NOTE: If you will be working on data from multiple databases, it is recommended that you rename the layer in the map (e.g., MUPOLYGON_VA for the data from Virginia). Alternatively, remove the layer once you are done with it, before adding the next one.
   - From the Soil Data Development Toolbox, gSSURGO Mapping Toolset, open the "Create Soil Map" tool
   - In the tool, set the following parameters:
     - Map Unit Layer: MUPOLYGON [or renamed layer]
     - SDV Folder = "Soil Qualities and Features"
     - SDV Attribute = "Hydrologic Soil Group"
     - Aggregation Method = "Dominant Condition"
   - Run the tool. A new layer symbolized on the Hydrologic Soil Group will appear.
   - Repeat as needed for MUPOLYGON data from different databases. 
   - Close ArcMap prior to attempting to run this function.
   
   The run of the above tool modifies the geodatabase by creating new tables with the prefix "SDV_". It does not modify existing tables.
   
   Given the above parameters, it creates a table named SDV_HydrolGrp_DCD, in which the field named HYDROLGRP_DCD contains the Hydrologic Soil Group code. If this is not the case, the function will fail.
   '''

   # Set up some variables
   mupolygon = in_GDB + os.sep + "MUPOLYGON"
   hydroTab = in_GDB + os.sep + "SDV_HydrolGrp_DCD"
   hydroFld = "HYDROLGRP_DCD"
   bname = os.path.basename(in_GDB)
   print("Working on %s" %bname) 
   
   # Process: Join Hydrologic Group value to MUPOLYGON
   # First check if field exists (due to prior processing) and delete if so
   fldList = arcpy.ListFields(mupolygon) 
   fldNames = [f.name for f in fldList]
   if hydroFld in fldNames:
      print("Deleting existing hydrologic group field in MUPOLYGON...")
      arcpy.DeleteField_management (mupolygon, hydroFld)
   print("Joining Hydrologic Group field to MUPOLYGON...")
   arcpy.JoinField_management(mupolygon, "MUKEY", hydroTab, "MUKEY", hydroFld)
   
   # Create and calculate a field in MUPOLYGON to store numeric version of hydrologic group, and calculate
   print("Adding numeric field for hydrologic group...")
   arcpy.AddField_management(mupolygon, "HydroGrpNum", "SHORT")
   
   print("Calculating HydroGrpNum field...")
   codeblock = '''def grpNum(fld):
      d = dict()
      d["A"] = 1
      d["B"] = 2
      d["C"] = 3
      d["D"] = 4
      
      if fld == None:
         key = "D"
      elif len(fld) > 1:
         key = fld[-1]
      else:
         key = fld   
      
      val = d[key]
      
      return val
   '''
   expression = "grpNum(!%s!)" %hydroFld
   arcpy.CalculateField_management (mupolygon, "HydroGrpNum", expression, 'PYTHON', codeblock)
   
   print("Mission complete for %s." %bname)

   return

def Kfactor_vec(in_GDB):
   '''To the MUPOLYGON feature class, adds a string field called "KFACTWS_DCD", containing the K-factor values extracted from gSSURGO. Also creates a new numeric field called "kFactor", needed for rasterization.
   
   This function modifies the input geodatabase by adding new tables and fields. It does not modify existing fields.
   
   Parameters:
   - in_GDB: input gSSURGO geodatabase
   
   K-factor values range from 0.02 to 0.69*. 
   
   * per https://dec.vermont.gov/sites/dec/files/wsm/stormwater/docs/StormwaterConstructionDischargePermits/sw_9020_Erodibility_%20Guidance.pdf

   IMPORTANT: Prior to running this function, new data must be created within the input geodatabase, using a tool in the Soil Data Development Toolbox. Tools in this toolbox can only be run from within ArcMap (not ArcPro) and I haven't figured out a way to automate this with a script, so you need to do it manually.
   
   The Soil Data Development Toolbox must be added to ArcToolbox in ArcMap, and is available from https://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/home/?cid=nrcs142p2_053628#tools.
   
   TO CREATE THE DATA NEEDED FOR THIS FUNCTION:
   - Within ArcMap, add the MUPOLYGON feature class as a layer. 
     - NOTE: If you will be working on data from multiple databases, it is recommended that you rename the layer in the map (e.g., MUPOLYGON_VA for the data from Virginia). Alternatively, remove the layer once you are done with it, before adding the next one.
   - From the Soil Data Development Toolbox, gSSURGO Mapping Toolset, open the "Create Soil Map" tool
   - In the tool, set the following parameters:
     - Map Unit Layer: MUPOLYGON [or renamed layer]
     - SDV Folder = "Soil Erosion Factors"
     - SDV Attribute = "K Factor, Whole Soil"
     - Aggregation Method = "Dominant Condition"
     - Top Depth (cm) = "0"
     - Bottom Depth (cm) = "10"
   - Run the tool. A new layer symbolized on the K-factor will appear.
   - Repeat as needed for MUPOLYGON data from different databases. 
   - Close ArcMap prior to attempting to run this function.
   
   The run of the above tool modifies the geodatabase by creating new tables with the prefix "SDV_". It does not modify existing tables.
   
   Given the above parameters, it creates a table named SDV_KfactWS_DCD_0to10, in which the field named KFACTWS_DCD contains the K-factor. If this is not the case, the function will fail.
   '''

   # Set up some variables
   mupolygon = in_GDB + os.sep + "MUPOLYGON"
   kfactTab = in_GDB + os.sep + "SDV_KfactWS_DCD_0to10"
   kfactFld = "KFACTWS_DCD"
   bname = os.path.basename(in_GDB)
   print("Working on %s" %bname) 
   
   # For some reason, the K-factor field created by the SSURGO toolbox is a string. 
   # Convert to double since this is needed for correct rasterization later.
   print("Converting string to double...")
   arcpy.AddField_management(kfactTab, "kFactor", "DOUBLE")
   expression = "float(!%s!)" %kfactFld
   arcpy.CalculateField_management (kfactTab, "kFactor", expression, 'PYTHON')
   kfactFld = "kFactor"
   
   # Process: Join K-factor value to MUPOLYGON
   # First check if field exists (due to prior processing) and delete if so
   fldList = arcpy.ListFields(mupolygon) 
   fldNames = [f.name for f in fldList]
   if kfactFld in fldNames:
      print("Deleting existing K-factor field in MUPOLYGON...")
      arcpy.DeleteField_management (mupolygon, kfactFld)
   print("Joining K-factor field to MUPOLYGON...")
   arcpy.JoinField_management(mupolygon, "MUKEY", kfactTab, "MUKEY", kfactFld)
   
   # Replace nulls in the K-factor field with the value 0.30, per the OpenNSPECT Technical Guide.
   print("Replacing nulls in K-factor field...")
   codeblock = '''def replaceNulls(fld):
      if fld == None:
         val = 0.3
      else:
         val = fld
      return val
   '''
   expression = "replaceNulls(!%s!)" %kfactFld
   arcpy.CalculateField_management (mupolygon, kfactFld, expression, 'PYTHON', codeblock)
   
   print("Mission complete for %s." %bname)

   return

def SSURGOtoRaster(in_gdbList, in_Fld, in_Snap, out_Raster):
   '''From one or more gSSURGO geodatabases, creates a raster representing values from a specified field in the MUPOLYGON feature class. 
   
   Parameters:
   in_gdbList: List of gSSURGO geodatabases containing added attributes
   in_Fld: field in MUPOLYGON feature class used to determine output raster values
   in_Snap: Input raster that determines output coordinate system, processing extent, cell size, and alignment
   out_Raster: Output raster 
   '''
   
   # Set overwrite to be true         
   arcpy.env.overwriteOutput = True
   
   # Specify scratch location
   scratchGDB = arcpy.env.scratchGDB
   
   # Empty list to contain raster paths
   rasterList = []
   
   # Work through loop converting polygons to rasters
   for gdb in in_gdbList:
      try:
         inPoly = gdb + os.sep + "MUPOLYGON"
         bname = os.path.basename(gdb).replace(".gdb","")
         print("Working on %s" %bname)
         outRast = scratchGDB + os.sep + bname
         PolyToRaster(inPoly, in_Fld, in_Snap, outRast)
         rasterList.append(outRast)
      except:
         print("Failed to rasterize %s" %bname)
   
   print("Finalizing output and saving...")
   finRast = CellStatistics(rasterList, "MAXIMUM", "DATA")
   finRast.save(out_Raster)
   
   print("Mission complete.")


### Functions for creating Soil Loss Potential, Runoff Potential, and Soil Sensitivity Scores
def SlopeTrans(in_Raster, inputType, transType, out_Trans, out_Slope, zfactor = 1):
   '''From a raster representing slope, creates a new raster representing a transformed representation of slope, depending on the transformation type (transType) specified. 
   
   The transformation types that may be specified are:
   - TRUNCLIN: A truncated linear function. Flat and nearly level slopes less than or equal to 1 degree (~2%) are scored 0, extreme slopes greater than or equal to 30 degrees (~58%) are scored 100, and values are scaled linearly in between the threshold values. This is a modification of the transformation used to derived the Slope Score in the 2017 edition of the ConservationVision Watershed Model.
   - TRUNCSIN: A truncated sine function. The sine of the angle is multiplied by 200 to get the score, but values above 100 are truncated, which happens at 30 degrees.
   - RUSLE: A stepwise sine function used to derive the slope steepness factor (S) in the RUSLE equation. (See equations 4-4 and 4-5 on page 107 of the RUSLE handbook.)

   Parameters:
   - in_Raster: input raster representing slope or elevation
   - inputType: indicates whether the input raster is slope in degrees (DEG or DEGREES), slope as percent grade (PERC or PERCENT), or elevation (ELEV or ELEVATION)
   - transType: the transformation function used to produce the output raster
     permitted values: TRUNCLIN, TRUNCSIN, RUSLE
   - out_Trans: output raster representing transformed slope
   - out_Slope: output raster representing slope as percent grade (ignored if input is a slope raster)
   - zfactor: Number of ground x,y units in one surface z-unit (ignored if input is a slope raster)
   '''
   
   # Make sure user entered valid parameters, and report what they are.   
   if inputType in ("DEG", "DEGREES"):
      slopeType = "DEGREES"
      print("Input is slope in degrees.")
   elif inputType in ("PERC", "PERCENT"):
      slopeType = "PERCENT"
      print("Input is slope as percent grade.")
   elif inputType in ("ELEV", "ELEVATION"):
      slopeType = "PERCENT"
      print("Input is elevation.")
   else:
      print("Input type specification is invalid. Aborting.")
      sys.exit()
      
   if transType == "TRUNCLIN":
      print("Appying the truncated linear transformation.")
   elif transType == "TRUNCSIN":
      print("Applying the truncated sine transformation.")
   elif transType == "RUSLE":
      print("Applying the RUSLE transformation to get the S-factor.")
   else:
      print("Transformation specification is invalid. Aborting.")
      sys.exit()
      
   # Set overwrite to be true         
   arcpy.env.overwriteOutput = True
   
   # Set scratch output location
   scratchGDB = arcpy.env.scratchGDB
   
   # Identify the slope raster or create it if necessary
   if inputType in ("ELEV", "ELEVATION"):
      print("Calculating slope from elevation...")
      in_Slope = Slope(in_Raster, "PERCENT_RISE", zfactor) 
      in_Slope.save(out_Slope)
   else:   
      in_Slope = Raster(in_Slope)
   
   if transType == "TRUNCLIN":
   # Set flat and nearly level slopes (LTE 1 degree) to 0. Set extreme slopes (GTE 30 degrees) to 100. Use linear function to scale between those values.
      minSlope = 1.0
      maxSlope = 30.0
      if slopeType == "PERCENT":
         print("Calculating score...")
         minSlope = 100*math.tan(minSlope*math.pi/180)
         maxSlope = 100*math.tan(maxSlope*math.pi/180)
         outRaster = Con(in_Slope <= minSlope, 0, Con((in_Slope > maxSlope), 100, 100 * (in_Slope - minSlope) / (maxSlope - minSlope)))
      else: 
         print("Calculating score...")
         outRaster = Con(in_Slope <= minSlope, 0, Con((in_Slope > maxSlope), 100, 100 * (in_Slope - minSlope) / (maxSlope - minSlope)))
   
   elif transType == "TRUNCSIN":
   # Take the sine, multiply by 200, and integerize. Upper values are truncated at 100 (which happens at 30 degrees).
      if slopeType == "PERCENT":
         print("Converting percent grade to radians and calculating score...")
         outRaster = Min(100, Int(0.5 + 200*Sin(ATan(in_Slope/100))))
      else: 
         print("Converting degrees to radians and calculating score...")
         outRaster = Min(100, Int(0.5 + 200*Sin(in_Slope * math.pi/180.0)))
         
   else:
   # Use RUSLE transformation equations
      inflect = 9.0
      if slopeType == "PERCENT":
         print("Converting percent grade to radians and calculating S-factor...")
         outRaster = Con(in_Slope < inflect, (10.8*(Sin(ATan(in_Slope/100))) + 0.03), (16.8*(Sin(ATan(in_Slope/100))) - 0.50))
      else: 
         inflect = math.atan(inflect/100)*180/math.pi
         print("Converting degrees to radians and calculating S-factor...")
         outRaster = Con(in_Slope < inflect, (10.8*(Sin(in_Slope * math.pi/180.0)) + 0.03), (16.8*(Sin(in_Slope * math.pi/180.0)) - 0.50))
      
   print("Saving output...")
   outRaster.save(out_Trans)
   
   print("Mission complete")
   
   return

def soilLoss_RKSC(in_Rfactor, in_Kfactor, in_Sfactor, in_Cfactor, out_RKSC):
   '''Multiplies the rasters representing three of the factors in the Revised Universal Soil Loss Equation (RUSLE), to produce a relative measure of the propensity for soil loss. Does not include the slope length (L)or supporting practices (P) factors. Inputs must have been first generated by previous functions to produce the input rasters.  

   Parameters:
   - in_Rfactor: input raster representing the rainfall/runoff erosivity factor
   - in_Kfactor: input raster representing the soil erodibility factor
   - in_Sfactor: input raster representing the slope steepness factor (output from SlopeTrans function with RUSLE option)
   - in_Cfactor: input raster or constant (float) representing the cover management factor
   - out_RKSC: output raster representing the product of R,K, and S factors
   
   NOTE: For the C-factor, the Watershed Impact Model uses a constant value (0.7) to get a measure of soil loss under a "worst-case scenario" of bare land. If you want soil loss under existing land cover conditions, first create a C-factor raster based on land cover (see OpenNSPECT documentation), and use that as the input. 
   
   This functions assumes all inputs are in the same coordinate system and properly aligned with each other.
   '''
   
   # Set overwrite to be true         
   arcpy.env.overwriteOutput = True
   
   # Calculate propensity for soil loss by multiplying the factors
   print("Calculating propensity for soil loss...")
   R = Raster(in_Rfactor)
   K = Raster(in_Kfactor)
   S = Raster(in_Sfactor)
   try:
      C = Raster(in_Cfactor)
   except:
      C = in_Cfactor
   RKSC = R*K*S*C
   
   print("Saving output...")
   RKSC.save(out_RKSC)
   
   print("Mission complete.")

def curvNum(in_LC, in_HydroGrp, out_CN):
   '''Given input land cover and hydrologic group, produces output raster representing runoff curve numbers.
   
   Curve numbers are assigned to combinations of land cover and soil types as specified in Table 1, page 6 of the OpenNSPECT Technical Guide. 

   Parameters:
   - in_LC: Input classified land cover raster, using standard NLCD land cover codes (updated with CCAP for code 32 = unconsolidated shore), OR an integer representing a desired land cover class (e.g., worst-case scenario bare land = 31)
   - in_HydroGrp: Input raster representing hydrologic groups (integer values must range from 1 = A to 4 = D)
   - out_CN: Output raster representing runoff curve numbers
   
   Note: If a land cover raster is used, the function modifies the input land cover attribute table, by adding and calculating a new field to store the curve numbers
   '''
   
   # Set overwrite to be true         
   arcpy.env.overwriteOutput = True
   
   # Set scratch output location
   scratchGDB = arcpy.env.scratchGDB
   
   # Initialize empty data dictionaries
   dictA = dict()
   dictB = dict()
   dictC = dict()
   dictD = dict()
   m = dict()
   
   # Populate dictionary for hydro group A, then append to list
   dictA[11] = 0
   dictA[21] = 49
   dictA[22] = 61
   dictA[23] = 77
   dictA[24] = 89
   dictA[31] = 77
   dictA[32] = 0
   dictA[41] = 30
   dictA[42] = 30
   dictA[43] = 30
   dictA[52] = 30
   dictA[71] = 30
   dictA[81] = 39
   dictA[82] = 67
   dictA[90] = 0
   dictA[95] = 0
   m["A"] = dictA
   
   # Populate dictionary for hydro group B
   dictB[11] = 0
   dictB[21] = 69
   dictB[22] = 75
   dictB[23] = 85
   dictB[24] = 92
   dictB[31] = 86
   dictB[32] = 0
   dictB[41] = 55
   dictB[42] = 55
   dictB[43] = 55
   dictB[52] = 48
   dictB[71] = 58
   dictB[81] = 61
   dictB[82] = 78
   dictB[90] = 0
   dictB[95] = 0
   m["B"] = dictB
   
   # Populate dictionary for hydro group C
   dictC[11] = 0
   dictC[21] = 79
   dictC[22] = 83
   dictC[23] = 90
   dictC[24] = 94
   dictC[31] = 91
   dictC[32] = 0
   dictC[41] = 70
   dictC[42] = 70
   dictC[43] = 70
   dictC[52] = 65
   dictC[71] = 71
   dictC[81] = 74
   dictC[82] = 85
   dictC[90] = 0
   dictC[95] = 0
   m["C"] = dictC
   
   # Populate dictionary for hydro group D
   dictD[11] = 0
   dictD[21] = 84
   dictD[22] = 87
   dictD[23] = 92
   dictD[24] = 95
   dictD[31] = 94
   dictD[32] = 0
   dictD[41] = 77
   dictD[42] = 77
   dictD[43] = 77
   dictD[52] = 73
   dictD[71] = 78
   dictD[81] = 80
   dictD[82] = 89
   dictD[90] = 0
   dictD[95] = 0
   m["D"] = dictD
      
   hydroGrps = ["A", "B", "C", "D"]
   in_HydroGrp = Raster(in_HydroGrp)
   
   if type(in_LC) == str:
      # Create and calculate curve number fields in the land cover attribute table
      for grp in hydroGrps:  
         fldName = "cn_%s" %grp
         d = m[grp]
         
         fldList = arcpy.ListFields(in_LC) 
         fldNames = [f.name for f in fldList]
         if fldName in fldNames:
            print("Deleting existing field %s..." %fldName)
            arcpy.DeleteField_management (in_LC, fldName)
         
         print("Adding field %s..." %fldName)
         arcpy.AddField_management(in_LC, fldName, "SHORT")
      
         print("Calculating field...")
         codeblock = '''def curvnum(code, dic):
            try:
               cn = dic[code]
            except:
               cn = 0
            return cn
            '''
         expression = "curvnum(!VALUE!, %s)" %d
         arcpy.CalculateField_management (in_LC, fldName, expression, 'PYTHON', codeblock)
      
      # Create a new raster from the curve number fields, based on soil type
      print("Creating curve number raster...")
      outRaster = Con(in_HydroGrp == 1, Lookup(in_LC, "cn_A"), Con(in_HydroGrp == 2, Lookup(in_LC, "cn_B"), Con(in_HydroGrp == 3, Lookup(in_LC, "cn_C"),Con(in_HydroGrp == 4, Lookup(in_LC, "cn_D")))))
   
   else:
      # Use the specified land cover constant with soil type to get the curve number
      outRaster = Con(in_HydroGrp == 1, dictA[in_LC], Con(in_HydroGrp == 2, dictB[in_LC], Con(in_HydroGrp == 3, dictC[in_LC],Con(in_HydroGrp == 4, dictD[in_LC]))))
   
   print("Saving output...")
   outRaster.save(out_CN)
   
   print("Mission complete.")

def eventRunoff(in_Raster, in_Rain, out_GDB, nameTag, cellArea = 1000000, inputType = "CN", convFact = 1, vol = 0):
   '''Produces an output raster representing event-based runoff depth in inches. Optionally produces runoff volume as well.
   
   Parameters:
   - in_Raster: input raster representing curve numbers OR maximum retention
   - in_Rain: input constant or raster representing rainfall
   - out_GDB: geodatabase for storing outputs
   - nameTag: tag to add to basenames to indicate land cover source
   - cellArea: area of cells in Curve Number raster, in square centimeters
   - inputType: indicates whether in_Raster is curve numbers (CN) or retention (RET)
   - convFact: conversion factor to convert input rainfall depth units to inches
   - vol: indicates whether a runoff volume raster should (1) or should not (0) be produced.
   '''
   # Set overwrite to be true         
   arcpy.env.overwriteOutput = True
   
   # Set scratch output location
   scratchGDB = arcpy.env.scratchGDB
   
   # Set up some variables
   in_Raster = Raster(in_Raster)
   try:
      in_Rain = Raster(in_Rain)
   except:
      pass
   out_Retention = out_GDB + os.sep + "Retention_%s" %nameTag
   out_runoffDepth = out_GDB + os.sep + "runoffDepth_%s" %nameTag
   out_runoffVolume = out_GDB + os.sep + "runoffVol_%s" %nameTag

   # Perform calculations
   # Result could be raster or a constant depending on input
   if convFact != 1:
      rain = convFact*in_Rain 
   else:
      rain = in_Rain
   
   if inputType == "CN":
      curvNum = in_Raster
      print("Calculating maximum retention...")
      # Have to deal with division by zero here.
      retention = Con(curvNum == 0, 1000, ((float(1000)/curvNum) - 10))
      print("Saving...")
      retention.save(out_Retention)
   else:
      retention = in_Raster
   
   print("Calculating runoff depth (inches)...")
   # Set runoff depth to zero if rainfall is less than initial abstraction
   runoffDepth = Con(curvNum == 0, 0, Con((rain - 0.2*retention) > 0,(rain - 0.2*retention)**2/(rain + 0.8*retention),0))
   print("Saving...")
   runoffDepth.save(out_runoffDepth)
   
   if vol == 1:
      print("Calculating runoff volume (liters)...")
      # 2.54 converts inches to cm
      # 0.001 converts cubic cm to liters
      volumeConversion = 0.00254*cellArea
      runoffVolume = volumeConversion*runoffDepth
      print("Saving...")
      runoffVolume.save(out_runoffVolume)
   
   print("Mission accomplished.")

def calcSoilSensScore(in_SoilLoss, in_Runoff, out_GDB, in_Mask = "NONE"):
   '''Creates a "Soil Sensitivity Score" raster, representing relative potential for impacts due to soil loss and rain runoff, on a scale from 1 (low sensitivity) to 100 (high sensitivity). Also creates component rasters from which this is derived, namely "Soil Loss Potential Score" and "Runoff Potential Score" rasters.
   
   NOTE: The input rasters representing potential for soil loss and runoff should be based on a "worst case scenario" land cover type, e.g., bare soil. These should have been created using functions in the procSSURGO.py script.
   
   Parameters:
   - in_SoilLoss: input raster representing relative potential for soil loss under "worst case" land cover conditions (output from soilLoss_RKSC function)
   - in_Runoff: input raster representing relative potential for runoff under "worst case" land cover conditions (output from eventRunoff function)
   - out_GDB: geodatabase to store outputs
   - in_Mask: Mask raster used to define the processing area
   '''
   
   # Set processing environment
   if in_Mask != "NONE":
      arcpy.env.mask = in_Mask
   
   # Set up outputs
   soilLossScore = out_GDB + os.sep + "soilLoss_Score"
   runoffScore = out_GDB + os.sep + "runoff_Score"
   sensScore = out_GDB + os.sep + "soilSens_Score"

   # Get truncation values
   print("Calculating raster cutoff values...")
   (slTruncMin, slTruncMax) = getTruncVals(in_SoilLoss, in_Mask)
   (roTruncMin, roTruncMax) = getTruncVals(in_Runoff, in_Mask)

   # Rescale the SoilLoss raster
   print("Rescaling soil loss potential...")
   Fx = TfLinear ("", "", slTruncMin, 1, slTruncMax, 100) 
   slScore = RescaleByFunction(in_SoilLoss, Fx, 1, 100)
   # print("Saving...")
   slScore.save(soilLossScore)
   
   # Rescale the Runoff raster
   print("Rescaling runoff potential...")
   Fx = TfLinear ("", "", roTruncMin, 1, roTruncMax, 100) 
   roScore = RescaleByFunction(in_Runoff, Fx, 1, 100)
   print("Saving...")
   roScore.save(runoffScore)

   # Take the average of the rescaled values
   print("Calculating soil sensitivity score...")
   sens = (slScore + roScore)/2
   print("Saving...")
   
   sens.save(sensScore)
   
   print("Mission accomplished.")


### Functions for creating Overland Flow, Karst Prevalence, and Landscape Position Scores
def makeHdwtrsIndicator(in_FlowLines, in_Catchments, in_BoundPoly, in_Mask, out_Hdwtrs):
   '''Creates a "Headwaters Indicator" raster, representing presence in a headwater (1) or non-headwater (0) catchment. This is a component of the Landscape Position Score.
   
   NOTE: Flowlines and catchments are assumed to come from NHDPlus-HR, which includes fields to identify headwater streams and link them to their corresponding catchments. It assumes that the field indicating headwater status is already attached to the NHDPlus feature class. If this is not the case this needs to be done first.
   
   Parameters:
   - in_FlowLines: Input NHDPlus Flowlines (line features)
   - in_Catchments: Input NHDPlus Catchments (polygon features)
   - in_BoundPoly: Input polygon feature class delimiting the area of interest
   - in_Mask: Input raster used to define the processing area, cell size, and alignment
   - out_Hdwtrs: Output raster representing headwater status
   '''
   
   # Set environment variables
   arcpy.env.snapRaster = in_Mask
   arcpy.env.cellSize = in_Mask
   arcpy.env.mask = in_Mask
   
   scratchGDB = arcpy.env.scratchGDB

   # Select the catchments intersecting in_BoundPoly, and save them to a temp feature class
   print
   tmpCatch = scratchGDB + os.sep + "tmpCatch"
   print("Selecting catchments within area of interest...")
   arcpy.MakeFeatureLayer_management(in_Catchments, "catch_lyr")
   arcpy.SelectLayerByLocation_management("catch_lyr", "intersect", in_BoundPoly)
   print("Copying subset...")
   arcpy.CopyFeatures_management("catch_lyr", tmpCatch)
   
   # Attach the headwaters indicator field to the catchment subset, then rasterize
   fldID = "NHDPlusID"
   fldHead = "StartFlag"
   print("Joining headwaters indicator field to catchments...")
   arcpy.JoinField_management(tmpCatch, fldID, in_FlowLines, fldID, fldHead)
   print("Rasterizing...")
   PolyToRaster(tmpCatch, fldHead, in_Mask, out_Hdwtrs) 
   
   print("Mission complete.")
   
   return out_Hdwtrs

def calcFlowScore(in_FlowLength, out_FlowScore, in_Hdwtrs = "NONE", minDist = 50, maxDist = 500, discount = 0.9):
   '''Creates a "Flow Distance Score" raster, in which cells within a specified flow distance to water are scored 100, and cells farther than a specified flow distance are scored 1. A headwaters indicator raster may optionally be used to "discount" cells not within headwater catchments. This is a component of the Landscape Position Score.
   
   NOTE: The FlowLength raster must have been derived from an overland flow direction raster (e.g., provided by NHDPlus.) 
   
   Parameters:
   - in_FlowLength: Input raster representing the overland flow distance to water
   - out_FlowScore: Output raster in which flow lengths have been converted to scores
   - in_Hdwtrs: Input raster indicating whether cells are within a headwater catchment (1) or not (0), used to "discount" non-headwater values. Can be set to "NONE" if no discount is desired.
   - minDist: The flow distance threshold below which the (non-discounted) score is set to 100.
   - maxDist: The flow distance threshold above which the (non-discounted) score is set to 1.
   - discount: A value multiplied by the initial score to get the final score (ignored if no headwaters raster is specified). If the initial score is 100, but the cell is not in a headwater catchment and the discount value is 0.9, the final score will be 90.
   '''
   
   # Set environment variables
   if in_Hdwtrs != "NONE":
      arcpy.env.mask = in_Hdwtrs
   
   # Rescale the flow length raster to scores 
   print("Rescaling flow lengths to scores...")
   Fx = TfLinear ("", "", minDist, 100, maxDist, 1) 
   flowScore = RescaleByFunction(in_FlowLength, Fx, 100, 1)
   
   # Discount scores
   if in_Hdwtrs == "NONE":
      finScore = flowScore
   else:
      print("Discounting non-headwater scores...")
      finScore = Con(Raster(in_Hdwtrs)==0, discount*flowScore, flowScore)
   print("Saving...")
   finScore.save(out_FlowScore)
   
   return finScore

def calcSinkScore(in_SinkPolys, fld_Area, procMask, clipMask, out_GDB, searchRadius = 5000):
   '''From input sinkhole polygons, generates three outputs:
   - A point feature class containing sinkhole centroids
   - A raster representing sinkhole density
   - A raster representing sinkhole scores from 1 to 100, derived from sinkhole density
   
     Parameters:
   - in_SinkPolys: Input polygons representing sinkhole features
   - fld_Area: Field representing the sinkhole area; used for "population" in kernel density calculation
   - procMask: Mask raster used to define the processing area, cell size, and alignment
   - clipMask: Mask raster used to define final output area
   - out_GDB: Geodatabase to store output products
   - searchRadius: Search radius used to calculate kernel density
   
   Note: Prior to running this function, make sure the sinkhole data are "clean", i.e., no overlaps/duplicates. There must also be a field representing the sinkhole area in desired units.
   '''
   
   # Set environment variables
   arcpy.env.mask = procMask
   arcpy.env.extent = procMask
   arcpy.env.snapRaster = procMask
   arcpy.env.cellSize = procMask
   
   # Set up outputs
   sinkPoints = out_GDB + os.sep + "sinkPoints"
   sinkPoints_prj = out_GDB + os.sep + "sinkPoints_prj"
   sinkDens = out_GDB + os.sep + "sinkDens"
   sinkScore = out_GDB + os.sep + "sinkScore"
   
   # Generate sinkhole centroids
   print("Generating sinkhole centroids...")
   arcpy.FeatureToPoint_management(in_SinkPolys, sinkPoints)
   
   # Run kernel density
   print("Calculating kernel density...")
   pts = ProjectToMatch_vec(sinkPoints, procMask, sinkPoints_prj, copy = 0)
   kdens = KernelDensity(pts, fld_Area, procMask, searchRadius, "HECTARES", "DENSITIES", "PLANAR")
   print("Saving...")
   kdens.save(sinkDens)
   
   # Convert kernel density to score
   print("Calculating truncation values...")
   msk = Con(kdens > 0, 1)
   (TruncMin, TruncMax) = getTruncVals(kdens, msk)
   TruncMax = int(TruncMax)
   Fx = TfLinear ("", "", 0, 1, TruncMax, 100)
   print("Truncation values set to 0, %s." %TruncMax)
   print("Converting kernel density to scores...")
   arcpy.env.mask = clipMask
   Score = RescaleByFunction(sinkDens, Fx, 1, 100)
   # print("Saving...")
   Score.save(sinkScore)
    
   return Score
   
   print("Mission accomplished.")

def calcKarstScore(in_KarstPolys, procMask, clipMask, out_GDB, in_SinkScore, minDist = 100, maxDist = 5000):
   '''From karst polygons and an optional sinkhole score raster, generates three or four outputs:
   - A raster representing karst polygons
   - A raster representing distance to karst
   - A raster representing distance scores from 1 to 100 (omitted if no density score raster is used)
   - A raster representing the final karst score from 1 to 100 
   
   Parameters:
   - in_KarstPolys: Polygons representing karst geology
   - procMask: Mask raster used to define the processing area, cell size, and alignment
   - clipMask: Mask raster used to define final output area
   - out_GDB: Geodatabase to store output products
   - in_SinkScore: Input raster representing a score from 0 to 100 based on density of sinkhole features. May be omitted (set to "NONE") for a simpler karst score based only on distance to karst geology.
   - minDist: Minimum distance to karst, below which the score is 100
   - maxDist: Maximum distance to karst, above which the score is 1
   
   Justification for minDist and maxDist default values:
   I calculated distance between sinkhole centroids and karst polygons. 97.3% of points were within 100 m of polygons. Only one point was greater than 5000 m from a polygon; the distance for that point was ~5200 m.
   '''

   # Set environment variables
   arcpy.env.mask = procMask
   arcpy.env.snapRaster = procMask
   arcpy.env.cellSize = procMask
   
   # Set up outputs
   karst_Raster = out_GDB + os.sep + "karst_Raster"
   karst_eDist = out_GDB + os.sep + "karst_eDist" 
   karst_distScore = out_GDB + os.sep + "karst_distScore" 
   KarstScore = out_GDB + os.sep + "KarstScore" 

   # Convert karst polygons to raster
   print("Converting karst polygons to raster...")
   PolyToRaster(in_KarstPolys, "OBJECTID", procMask, karst_Raster)

   # Get Euclidean Distance and Distance Score
   print("Getting Euclidean distance to karst...")
   edist = EucDistance(karst_Raster, "", arcpy.env.cellSize)
   print("Saving...")
   edist.save(karst_eDist)
   print("Converting distances to scores...")
   arcpy.env.mask = clipMask
   Fx = TfLinear ("", "", minDist, 100, maxDist, 1) 
   edistScore = RescaleByFunction(karst_eDist, Fx, 100, 1)
 
   if in_SinkScore != "NONE":
      print("Saving...")
      edistScore.save(karst_distScore)

      # Get final Karst Score
      print("Calculating final Karst Score from distance and density scores...")
      combinedScore = CellStatistics([in_SinkScore, edistScore], "MEAN", "DATA")
      print("Saving...")
      combinedScore.save(KarstScore)
      return combinedScore
   
   else:
      print("Karst score is based only on Euclidean distance. Saving...")
      edistScore.save(KarstScore)
      return edistScore
   print("Mission accomplished.")

def calcPositionScore(in_FlowScore, in_KarstScore, out_PositionScore):
   '''Creates a "Landscape Position Score" raster, representing relative importance to stream health based on position in the landscape. It is the maximum of the Karst Score and the Flow Distance Score.

   Parameters:
   - in_FlowScore: Input raster representing the Flow Distance Score
   - in_KarstScore: Input raster representing the Karst Score
   - out_PositionScore: Output raster representing the Landscape Position Score
   '''
   
   # Calculate score
   print("Calculating Landscape Position Score...")
   score = CellStatistics([in_FlowScore, in_KarstScore], "MAXIMUM", "DATA")
   print("Saving...")
   score.save(out_PositionScore)
   
   print("Mission accomplished.")
   return score

def calcImpactScore(in_PositionScore, in_SoilSensScore, out_ImpactScore):
   '''Creates a raster representing the potential impact, based on landscape position and soil sensitivity. 
   
   Parameters:
   - in_PositionScore: Input raster representing relative importance based on landscape position
   - in_SoilSensScore: Input raster representing relative importance based on soil sensitivity
   - out_ImpactScore: The output raster representing either the Impact Score (if there is no input Importance Score) or the General Priority Score (if there is an input Importance Score used to adjust the Impact Score).
   '''

   print("Calculating Impact Score...")
   score = CellStatistics([in_PositionScore, in_SoilSensScore], "MEAN", "DATA")
      
   print("Saving...")
   score.save(out_ImpactScore)
   
   print("Mission accomplished.")
   return out_ImpactScore
   
def calcImportanceScore(in_raList, in_Snap, out_Raster):
   '''Calculates an Importance Score, based on polygon features identifying areas impacting resources of interest (e.g., catchments for Healthy Waters sites, assessment zones for drinking water intakes, or a Stream Conservation Site delineation).
   
   TO DO (?): Add a weighting field for individual features (e.g., count of EOs in SCS, number served by drinking water sources.) Will then also need to rescale each input, or else require rescale prior to input.
   
   Parameters:
   - in_raList: Input list of tuples. Each tuple consists of a polygon feature class representing the impact areas (which may overlap) and the weight assigned to those features.
   - in_Snap: An input raster used to set output cell size and alignment
   - out_Raster: Output raster representing the importance for protection/restoration
   
   NOTE: Areas of the output raster not covered by any polygon will be null. This could be appropriate if you are ONLY interested in the drainage areas of specific resources (e.g., for Healthy Waters prioritization). If you want a non-null, non-zero value throughout the "background" cells of the study area (e.g., for a seamless statewide Watershed Model), you can add a bounding polygon to cover the entire area; typically it would make sense to give a weight of 1 to this layer.
   '''
   
   # Set environment
   scratchGDB = arcpy.env.scratchGDB
   
   procList = list()
   for t in in_raList:
      (polys, weight) = t
      nameTag = os.path.basename(polys)
      print("Working on %s..."%nameTag)
      outCount = scratchGDB + os.sep + "polyCounts_%s"%nameTag
      countRast = scratchGDB + os.sep + "countRast_%s"%nameTag
      weightRast = scratchGDB + os.sep + "weightRast_%s"%nameTag
      
      # Count the overlapping features
      print("Counting overlaps...")
      arcpy.CountOverlappingFeatures_analysis(polys, outCount)
      
      # Rasterize
      print("Rasterizing...")
      PolyToRaster(outCount, "COUNT_", in_Snap, countRast)
      
      # Multiply by weight
      if weight != 1:
         print("Weighting...")
         wtCt = weight*Raster(countRast)
         wtCt.save(weightRast)
      else:
         wtCt = countRast
      
      # Append to processing list
      procList.append(str(wtCt))
   
   # Get weighted sum of resources
   if len(procList) > 1:
      print("Calculating weighted sum of resources...")
      wtSum = CellStatistics(procList, "SUM", "DATA")
   else:
      wtSum = Raster(procList[0])
      
   # Rescale weighted sums to scores
   rMax = wtSum.maximum
   
   print("Rescaling to scores...")
   score = 100.0*wtSum/rMax
   
   print("Saving...")
   score.save(out_Raster)
   
   print("Mission complete.")
   return out_Raster

def ScenarioScore(in_Case, in_WorstCase, in_BestCase, priorType, out_Score, in_Mask = "NONE"):
   '''Creates a raster representing a "Scenario Score", depending on how the values in the input raster compare to best- and worst-case scenarios for the same variable. 
   
   *** Not sure I will ever use this function.
   
   Parameters:
   - in_Case: Input raster representing the values of a particular variable of interest
   - in_WorstCase: Input raster representing the "worst case" scenario of the variable of interest
   - in_BestCase: Input raster representing the "best case" scenario of the variable of interest
   - priorType: Indicates whether the score is developed from a conservation (CONS) or restoration/stormwater management (REST) perspective
   - out_Score: Output raster representing a score for the variable of interest
   - in_Mask: Mask raster used to define the processing area
   '''
   
   # Set processing environment
   if in_Mask != "NONE":
      arcpy.env.mask = in_Mask
   
   # Make raster objects
   Case = Raster(in_Case)
   WorstCase = Raster(in_WorstCase)
   BestCase = Raster(in_BestCase)
   
   # Calculate
   if priorType == "CONS":
      print("Calculating scenario score for conservation...")
      score = 100.0*(WorstCase - Case)/(WorstCase - BestCase)
   else:
      print("Calculating scenario score for restoration or management...")
      score = 100.0*(Case - BestCase)/(WorstCase - BestCase)
   
   adjScore = Con(score > 100, 100, Con(score < 0, 0, score))
   adjScore.save(out_Score)

def calcPriorityScores(in_ImpactScore, in_ImportanceScore, in_ConsMask, in_RestMask, in_MgmtMask, out_GDB, rescale = "SLICE", slice = 10, nameTag = "NONE"):
   '''Calculates priorities for conservation, restoration, and stormwater management, depending on landcover type.
   
   Parameters:
   - in_ImpactScore: Input raster representing potential impact
   - in_ImportanceScore: Input raster representing relative importance based on the number of resources of interest that could be impacted by changes at each location. Alternatively, enter 1 to apply constant importance throughout.
   - in_ConsMask: Input raster representing lands that should get conservation priorities
   - in_RestMask: Input raster representing lands that should get restoration priorities
   - in_MgmtMask: Input raster representing lands that should get stormwater management priorities
   - out_GDB: Geodatabase to hold final outputs
   - rescale: Indicates whether to first rescale the Impact Score raster. Options: "SLICE" (to slice into specified number of quantiles), "STANDARD" (to do a standard linear rescale), or "NONE" (to use the raw impact scores).
   - slice: If the SLICE option is used for rescaling, the number of slices (i.e, quantiles)
   - nameTag: String to add as a suffix to standard output names. Set to "NONE" if no suffix is desired.
   '''
   
   # Set up outputs
   gPrior = out_GDB + os.sep + "genPriority"
   cPrior = out_GDB + os.sep + "consPriority"
   rPrior = out_GDB + os.sep + "restPriority"
   mPrior = out_GDB + os.sep + "mgmtPriority"
   
   if nameTag != "NONE":
      for p in [gPrior, cPrior, rPrior, mPrior]:
         p = p + "_%s" %nameTag
   #bName = os.path.basename(in_ImpactScore)
   
   # Calculate the general priority score
   print("Calculating General Priority Score...")
   if in_ImportanceScore == 1:
      p_Score = Raster(in_ImpactScore)
   else:
      p_Score = Raster(in_ImportanceScore)/100.0 * Raster(in_ImpactScore)
      p_Score.save(gPrior)
   
   # Rescale the general priority scores, if specified
   if rescale == "SLICE":
      print("Slicing priority scores into quantiles...")
      score = Slice(p_Score, slice, "EQUAL_AREA", 1)
      outPath = out_GDB + os.sep + "genPriority_slice"
      score.save(outPath)
   elif rescale == "STANDARD":
      print("Calculating minimum and maximum impact scores...")
      r = p_Score
      rMin = r.minimum
      rMax = r.maximum
      print("Rescaling...")
      Fx = TfLinear ("", "", rmin, 1, rmax, 100) 
      score = RescaleByFunction(r, Fx, 1, 100)
      outPath = out_GDB + os.sep + "genPriority_rscl"
      score.save(outPath)
   else:
      score = p_Score
      
   # Create priority rasters
   print("Creating conservation priority raster...")
   cons = Con(in_ConsMask, score)
   cons.save(cPrior)
   
   print("Creating restoration priority raster...")
   rest = Con(in_RestMask, score)
   rest.save(rPrior)
   
   print("Creating stormwater management priority raster...")
   mgmt = Con(in_MgmtMask, score)
   mgmt.save(mPrior)
   
   print("Mission accomplished.")
   return (cons, rest, mgmt)  
