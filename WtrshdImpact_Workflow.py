# ---------------------------------------------------------------------------
# WtrshdImpact_Workflow.py
# Version: ArcPro / Python 3+
# Creation Date: 2020-07-06
# Last Edit: 2022-05-25
# Creator: Kirsten R. Hazler
#
# Summary: Workflow to produce the ConservationVision Watershed Impact Model. This is intended as a guide for how to reproduce the model. Do NOT attempt to run this all at once. You would need to download the latest data, perform some manual operations, and replace hard-coded data paths with your own. It's recommended to comment out sections as needed, and just run small sections at a time. 

# For background and references see:
# - Coastal Services Center. 2014. “Technical Guide for OpenNSPECT, Version 1.2.” National Oceanic and Atmospheric Administration (NOAA). https://coast.noaa.gov/digitalcoast/tools/opennspect
# - Renard, K.G., G.R. Foster, G.A. Weesies, D.K. McCool, and D.C. Yoder. 1997. “Predicting Soil Erosion by Water: A Guide to Conservation Planning with the Revised Universal Soil Loss Equation (RUSLE).” Agriculture Handbook 703. U.S. Department of Agriculture. 
# - Cronshey, R., R.H. McCuen, N. Miller, W. Rawls, S. Robbins, and D. Woodward. 1986. “Urban Hydrology for Small Watersheds (2nd Ed.).” Technical Release 55. Natural Resources Conservation Service, U.S. Department of Agriculture.

# Also refer to methods in the associated technical report which can be found at https://www.dcr.virginia.gov/natural-heritage/vaconviswater

# ---------------------------------------------------------------------------
# Import modules
import HelperPro
from HelperPro import *

import WtrshdImpact_Functions
from WtrshdImpact_Functions import *

def main():
   
   ### Input Data ###
   
   # Processing geodatabase and masks
   procGDB = "N:\ProProjects\WatershedModels\CVWIM_Data\ProcessedData.gdb" # For most processed data
   clpShp = procGDB + os.sep + "HW_template_buff13k_noHoles"
   MaskNoWater = procGDB + os.sep + "mskNoWater_2016"
   procMask = r"Y:\SpatialData\SnapMasks\procMask50_conus.tif"
   jurisbnd = r"N:\SpatialData\VDOT\FromVDOT.gdb\JURISDICTIONS_albers"
   
   # Soils geodatabases, downloaded by state from gSSURGO site
   # Also need the Soil Data Development Toolbox
   # *** Manual operations required prior to function runs ***  
   # --> In each soils GDB, create the table "SDV_KfactWS_DCD_0to10" assigning K-factor values to soil map units. See the documentation in the commented section of the Kfactor_vec function for how to do this. 
   # --> In each soils GDB, create the table "SDV_HydrolGrp_DCD" assigning hydrologic soil groups to soil map units. See the documentation in the commented section of the HydroGrp_vec function for how to do this.
   dc_gdb = r"Y:\SpatialData\NRCS\SSURGO\gSSURGO_DC\gSSURGO_DC.gdb"
   de_gdb = r"Y:\SpatialData\NRCS\SSURGO\gSSURGO_DE\gSSURGO_DE.gdb"
   ky_gdb = r"Y:\SpatialData\NRCS\SSURGO\gSSURGO_KY\gSSURGO_KY.gdb"
   md_gdb = r"Y:\SpatialData\NRCS\SSURGO\gSSURGO_MD\gSSURGO_MD.gdb"
   nc_gdb = r"Y:\SpatialData\NRCS\SSURGO\gSSURGO_NC\gSSURGO_NC.gdb"
   pa_gdb = r"Y:\SpatialData\NRCS\SSURGO\gSSURGO_PA\gSSURGO_PA.gdb"
   tn_gdb = r"Y:\SpatialData\NRCS\SSURGO\gSSURGO_TN\gSSURGO_TN.gdb"
   va_gdb = r"Y:\SpatialData\NRCS\SSURGO\gSSURGO_VA\gSSURGO_VA.gdb"
   wv_gdb = r"Y:\SpatialData\NRCS\SSURGO\gSSURGO_WV\gSSURGO_WV.gdb"
   gdbList = [dc_gdb, de_gdb, ky_gdb, md_gdb, nc_gdb, pa_gdb, tn_gdb, va_gdb, wv_gdb] # for use in loops
   
   # Rainfall points from Dam Safety's PMP tool 
   # Toolbox and associated data from https://www.dcr.virginia.gov/dam-safety-and-floodplains/pmp-tool
   # *** Manual operation required prior to function run ***
   # --> For this model I used the PMP tool from within ArcGIS Pro to generate the rainfall points. I set the extent to run statewide, specified a 24-hour storm duration, and used the "General" output points
   pmpPath = r"Y:\SpatialData\DCR_DamSafety\PMP\pmpEvalTool_v2\Output\General\PMP_64457.gdb"
   pmpPts = pmpPath + os.sep + "General_PMP_Points_64457" # output from the PMP tool
   pmpFld = "PMP_24" # field in pmpPts
   
   # NHDPlus data and derivatives
   NHDPlus_path = r"N:\SpatialData\NHD_Plus\HydroNet\VA_HydroNetHR"
   FlowLines = NHDPlus_path + os.sep + r"VA_HydroNetHR.gdb\HydroNet\NHDFlowline"
   Catchments = NHDPlus_path + os.sep + r"VA_HydroNetHR.gdb\NHDPlusCatchment"
   
   # Other inputs
   in_Elev = r"N:\SpatialData\3DEP\Elev_cm.gdb\elev_cm_VA" # derived from 3DEP data
   in_Rfactor = r"N:\ProProjects\WatershedModels\CVWIM_Data\R_Factor\R-Factor_CONUS.tif" # downloaded from https://coast.noaa.gov/data/digitalcoast/zip/R-Factor-CONUS.zip
   KarstPolys = r"Y:\SpatialData\USGS\USKarstMap\USKarstMap.gdb\Contiguous48\Carbonates48" # downloaded from https://pubs.usgs.gov/of/2014/1156/
   SinkPolys = r"Y:\SpatialData\DMME\Sinkholes_VaDMME_2020SeptCurrent\Sinkholes_VaDMME.shp" # obtained by request from DMME

   ### End Input Data ###


   ### Soil Loss Potential procedures ###
   print("Starting procedures for soil loss potential...")
   
   ## Prepare the R-factor raster
   # Requires coarse-scale R-factor raster (in_Rfactor) from NOAA
   print("Downscaling R-factor raster...")
   Rfactor = procGDB + os.sep + "rusleR" # downscaled raster
   Downscale_ras(in_Rfactor, in_Elev, Rfactor, "BILINEAR", clpShp)
   arcpy.management.BuildPyramids(Rfactor)
   print("Downscaled R-factor raster complete.")
   
   ## Prepare the K-factor raster
   # Requires pre-processed SSURGO geodatabase(s) listed in gdbList
   for gdb in gdbList:
      print("Processing K-factor for %s..."%gdb)
      Kfactor_vec(gdb)
   print("K-factors complete.")
   Kfactor = procGDB + os.sep + "rusleK"
   print("Rasterizing K-factors...")
   SSURGOtoRaster(gdbList, "kFactor", in_Elev, Kfactor)
   arcpy.management.BuildPyramids(Kfactor)
   print("K-factor raster complete.")
   
   ## Prepare the S-factor raster
   # Requires an elevation raster (in_Elev)
   # In this case we used a DEM with elevation in cm, hence the zfactor below.
   Sfactor = procGDB + os.sep + "rusleS"
   slope_perc = procGDB + os.sep + "slope_perc"
   print("Transforming slope to S-factor...")
   SlopeTrans(in_Elev, "ELEV", "RUSLE", Sfactor, slope_perc, zfactor = 0.01)
   arcpy.management.BuildPyramids(Sfactor)
   print("S-factor raster complete.")
   
   ## No need for C-factor raster; we are assuming worst-case scenario, bare soil, so use a constant value (from OpenNSPECT)
   Cfactor = 0.7
   
   ## Prepare the soil loss raster based on RUSLE factors created above, assuming bare soil
   soilLoss_bare = procGDB + os.sep + "rusleRKSC_bare"
   print("Generating barren land soil loss raster...")
   soilLoss_RKSC(Rfactor, Kfactor, Sfactor, Cfactor, soilLoss_bare)
   arcpy.management.BuildPyramids(soilLoss_bare)
   print("Barren land soil loss raster complete.")
   
   
   ### Runoff Potential procedures ###
   print("Starting procedures for runoff potential...")
   
   ## Prepare the rainfall raster: Probable Maximum Precipitation (PMP)
   # Requires rainfall points generated from Dam Safety's PMP tool (pmpPts, with attributes in pmpFld)
   # Note: Interpolated to coarse scale first b/c memory was failing otherwise.
   maxPrecip250 = procGDB + os.sep + "maxPrecip_gen24_topo250"
   maxPrecip10 = procGDB + os.sep + "maxPrecip_gen24_topo10"
   print("Interpolating rainfall points...")
   interpPoints(pmpPts, pmpFld, in_Elev, maxPrecip250, clpShp, "TOPO", "", "", 250) # interpolate 
   print("Downscaling interpolated rainfall raster...")
   Downscale_ras(maxPrecip250, in_Elev, maxPrecip10, "BILINEAR", clpShp) # downscale
   arcpy.management.BuildPyramids(maxPrecip10)
   print("Rainfall raster complete.")
   
   
   ## Prepare the hydro group raster
   # Requires pre-processed SSURGO geodatabase(s) listed in gdbList
   for gdb in gdbList:
      print("Processing hydro group for %s..."%gdb")
      HydroGrp_vec(gdb)
   print("Hydro groups complete.")
   hydroGrp = procGDB + os.sep + "HydroGroup"
   print("Rasterizing hydro groups...")
   SSURGOtoRaster(gdbList, "HydroGrpNum", in_Elev, hydroGrp)
   arcpy.management.BuildPyramids(hydroGrp)
   print("Hydro group raster complete.")
   
   ## Prepare the curve number raster
   # Curve number based on soil hydro group only, assuming bare soil (NLCD code 31)
   curvNum_bare = procGDB + os.sep + "curvNum_bare"
   print("Creating barren land curve number raster...")
   curvNum(31, hydroGrp, curvNum_bare) 
   arcpy.management.BuildPyramids(curvNum_bare)
   print("Barren land curve number raster complete.")
   
   ## Prepare the runoff raster based on SCS curve-number method, assuming bare soil
   print("Creating barren land runoff raster...")
   eventRunoff(curvNum_bare, maxPrecip10, procGDB, "bare")
   # --> Relevant output is runoffDepth_bare
   runoffDepth_bare = procGDB + os.sep + "runoffDepth_bare"
   arcpy.management.BuildPyramids(runoffDepth_bare)
   print("Barren land runoff raster complete.")

   
   ### Calculate Soil Loss Potential, Runoff Potential, and Soil Sensitivity Scores ###
   print("Creating score rasters for soil loss, runoff, and soil sensitivity...")
   calcSoilSensScore(soilLoss_bare, runoffDepth_bare, procGDB, MaskNoWater)
   # --> Function outputs in procGDB are soilLoss_Score, runoff_Score, and soilSens_Score
   soilLossScore = procGDB + os.sep + "soilLoss_Score_bare"
   runoffScore = procGDB + os.sep + "runoff_Score_bare"
   SoilSensScore = procGDB + os.sep + "soilSens_Score_bare"
   for ras in [soilLossScore, runoffScore, SoilSensScore]:
      arcpy.management.BuildPyramids(ras)
   print("Soil score rasters complete.")
   
   
   ### Overland Flow procedures ###
   print("Starting procedures for overland flow...")
   
   ## Prepare the flow distance raster
   FlowLength = procGDB + os.sep + "overlandFlowLength"
   #[David - insert function call(s) to produc FlowLength here, and function definition(s) in Wtrshd_Functions script]
   
   ## Prepare the headwaters raster
   print("Creating headwaters raster...")
   Headwaters = procGDB + os.sep + "Hdwtrs"
   makeHdwtrsIndicator(FlowLines, Catchments, clpShp, MaskNoWater, Headwaters)
   arcpy.management.BuildPyramids(Headwaters)
   print("Headwaters raster complete.")
   
   ## Calculate Overland Flow score ##
   FlowScore = procGDB + os.sep + "FlowScore"
   calcFlowScore(FlowLength, FlowScore, Headwaters)
   arcpy.management.BuildPyramids(FlowScore)
   
   
   ### Karst procedures ###
   
   ## Prepare sinkholes
   # --> Manual operation required: Prior to running density scoring function, make sure the sinkhole data are "clean", i.e., no overlaps/duplicates. There must also be a field representing the sinkhole area in desired units (i.e., square meters).
   
   ## Calculate the sinkhole density score
   calcSinkScore(SinkPolys, "SqMeters", procMask, MaskNoWater, procGDB, searchRadius = 5000)
   # --> Function outputs in procGDB are sinkPoints, sinkPoints_prj, sinkDens, and sinkScore
   sinkPoints = procGDB + os.sep + "sinkPoints"
   sinkPoints_prj = procGDB + os.sep + "sinkPoints_prj"
   sinkDens = procGDB + os.sep + "sinkDens"
   sinkScore = procGDB + os.sep + "sinkScore"
   
   ## Calculate Karst Prevalence Score
   calcKarstScore(KarstPolys, procMask, MaskNoWater, procGDB, SinkScore, minDist = 100, maxDist = 5000)
   # --> Function outputs in procGDB are karst_Raster, karst_eDist, karst_distScore, karst_Score
   karst_Raster = procGDB + os.sep + "karst_Raster"
   karst_eDist = procGDB + os.sep + "karst_eDist" 
   karst_distScore = procGDB + os.sep + "karst_distScore" 
   KarstScore = procGDB + os.sep + "KarstScore" 
   arcpy.management.BuildPyramids(KarstScore)
   
   
   ### Landscape Position Score ###
   PositionScore = procGDB + os.sep + "PositionScore"
   calcPositionScore(FlowScore, KarstScore, PositionScore)
   arcpy.management.BuildPyramids(PositionScore)
   
   
   ### Potential Impact Score ###
   ImpactScore = procGDB + os.sep + "ImpactScore"
   calcImpactScore(PositionScore, SoilSensScore, ImpactScore)
   arcpy.management.BuildPyramids(ImpactScore)
   
   
   ### Finalize TIF Products ###
   mask = jurisbnd
   outPath = r"N:\ProProjects\WatershedModels\CVWIM_Data\cvwimTIF"
   procList = [runoffScore, soilLossScore, soilSens_Score, FlowScore, KarstScore, PositionScore, ImpactScore]
   
   for ras in procList:
      finalize_gdbRas2Tif(ras, mask, outPath)
   
if __name__ == "__main__":
   main()
   
