# ---------------------------------------------------------------------------
# WtrshdImpact_Workflow.py
# Version: ArcPro / Python 3+
# Creation Date: 2020-07-06
# Last Edit: 2022-05-24
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
   
   # Inputs - Soils geodatabases, downloaded by state from gSSURGO site
   # Also need the Soil Data Development Toolbox
   dc_gdb = r"E:\SpatialData\SSURGO\gSSURGO_DC\gSSURGO_DC.gdb"
   de_gdb = r"E:\SpatialData\SSURGO\gSSURGO_DE\gSSURGO_DE.gdb"
   ky_gdb = r"E:\SpatialData\SSURGO\gSSURGO_KY\gSSURGO_KY.gdb"
   md_gdb = r"E:\SpatialData\SSURGO\gSSURGO_MD\gSSURGO_MD.gdb"
   nc_gdb = r"E:\SpatialData\SSURGO\gSSURGO_NC\gSSURGO_NC.gdb"
   pa_gdb = r"E:\SpatialData\SSURGO\gSSURGO_PA\gSSURGO_PA.gdb"
   tn_gdb = r"E:\SpatialData\SSURGO\gSSURGO_TN\gSSURGO_TN.gdb"
   va_gdb = r"E:\SpatialData\SSURGO\gSSURGO_VA\gSSURGO_VA.gdb"
   wv_gdb = r"E:\SpatialData\SSURGO\gSSURGO_WV\gSSURGO_WV.gdb"
   gdbList = [dc_gdb, de_gdb, ky_gdb, md_gdb, nc_gdb, pa_gdb, tn_gdb, va_gdb, wv_gdb]
   
   # NHDPlus data and derivatives
   NHDPlus_path = r"N:\SpatialData\NHD_Plus\HydroNet\VA_HydroNetHR"
   FlowLines = NHDPlus_path + os.sep + r"VA_HydroNetHR.gdb\HydroNet\NHDFlowline"
   Catchments = NHDPlus_path + os.sep + r"VA_HydroNetHR.gdb\NHDPlusCatchment"
   
   # Other inputs
   in_Elev = r"N:\SpatialData\3DEP\Elev_cm.gdb\elev_cm_VA" # derived from 3DEP data
   KarstPolys = r"Y:\SpatialData\USGS\USKarstMap\USKarstMap.gdb\Contiguous48\Carbonates48" # downloaded from https://pubs.usgs.gov/of/2014/1156/
   SinkPolys = r"Y:\SpatialData\DMME\Sinkholes_VaDMME_2020SeptCurrent\Sinkholes_VaDMME.shp" # obtained by request from DMME

   ### End Input Data ###


   ### Soil Loss Potential procedures ###
   
   ## Prepare the R-factor raster
   # - Download/extract R-factor raster (in_Rfactor) provided by NOAA: https://coast.noaa.gov/data/digitalcoast/zip/R-Factor-CONUS.zip
   # - Downscale the R-factor raster to 10-m resolution
   in_Rfactor = r"N:\ProProjects\WatershedModels\CVWIM_Data\R_Factor\R-Factor_CONUS.tif" # input downloaded raster
   Rfactor = procGDB + os.sep + "rusleR" # downscaled raster
   Downscale_ras(in_Rfactor, in_Elev, Rfactor, "BILINEAR", clpShp)
   
   ## Prepare the K-factor raster
   # Manual operation --> In each soils GDB, create the table "SDV_KfactWS_DCD_0to10" assigning K-factor values to soil map units. See the documentation in the commented section of the Kfactor_vec function for how to do this. Then run the functions below.
   for gdb in gdbList:
      Kfactor_vec(gdb)
   Kfactor = procGDB + os.sep + "rusleK"
   SSURGOtoRaster(gdbList, "kFactor", in_Elev, Kfactor)
   
   ## Prepare the S-factor raster
   # Get a digital elevation model. In this case we used one with elevation in cm, hence the zfactor below.
   Sfactor = procGDB + os.sep + "rusleS"
   slope_perc = procGDB + os.sep + "slope_perc"
   SlopeTrans(in_Elev, "ELEV", "RUSLE", Sfactor, slope_perc, zfactor = 0.01)
   
   ## No need for C-factor raster; we are assuming worst-case scenario, bare soil, so use a constant value (from OpenNSPECT)
   Cfactor = 0.7
   
   ## Prepare the soil loss raster based on RUSLE factors, assuming bare soil
   soilLoss_bare = procGDB + os.sep + "rusleRKSC_bare"
   soilLoss_RKSC(Rfactor, Kfactor, Sfactor, Cfactor, soilLoss_bare)
   
   
   ### Runoff Potential procedures ###
   
   ## Prepare the rainfall raster: Probable Maximum Precipitation (PMP)
   # Manual operation --> For this model I used the PMP tool(https://www.dcr.virginia.gov/dam-safety-and-floodplains/pmp-tool) from within ArcGIS Pro to generate the points used for interpolation. I set the extent to run statewide, specified a 24-hour storm duration, and used the "General" output as input to the interpolation function below.
   pmpPts = procGDB + os.sep + "General_PMP_Points_64457" # output from the PMP tool
   pmpFld = "PMP_24" # field in pmpPts
   maxPrecip250 = procGDB + os.sep + "maxPrecip_gen24_topo250"
   maxPrecip10 = procGDB + os.sep + "maxPrecip_gen24_topo10"
   interpPoints(pmpPts, pmpFld, in_Elev, maxPrecip250, clpShp, "TOPO", "", "", 250) # interpolate points
   Downscale_ras(maxPrecip250, in_Elev, maxPrecip10, "BILINEAR", clpShp) # downscale
   # Note: Interpolated to coarse scale first b/c memory was failing otherwise.
   
   ## Prepare the hydro group raster
   # Manual operation --> In each soils GDB, create the table "SDV_HydrolGrp_DCD" assigning hydrologic soil groups to soil map units. See the documentation in the commented section of the HydroGrp_vec function for how to do this. Then run the functions below.
   for gdb in gdbList:
      HydroGrp_vec(gdb)
   hydroGrp = procGDB + os.sep + "HydroGrp"
   SSURGOtoRaster(gdbList, "HydroGrpNum", in_Elev, hydroGrp)
   
   ## Prepare the curve number raster
   curvNum_bare = procGDB + os.sep + "curvNum_bare"
   curvNum(31, hydroGrp, curvNum_bare) # Curve number under assumption of worst-case land cover: bare soil = 31
   
   ## Prepare the runoff raster based on SCS curve-number method, assuming bare soil
   eventRunoff(curvNum_bare, maxPrecip10, procGDB, "bare")
   # --> Relevant output is runoffDepth_bare
   runoffDepth_bare = procGDB + os.sep + "runoffDepth_bare"
   
   
   ### Calculate Soil Loss Potential, Runoff Potential, and Soil Sensitivity Scores ###
   calcSoilSensScore(soilLoss_bare, runoffDepth_bare, procGDB, MaskNoWater)
   # --> Function outputs in procGDB are soilLoss_Score, runoff_Score, and soilSens_Score
   soilLossScore = procGDB + os.sep + "soilLoss_Score_bare"
   runoffScore = procGDB + os.sep + "runoff_Score_bare"
   SoilSensScore = procGDB + os.sep + "soilSens_Score_bare"
   
   
   ### Overland Flow procedures ###
   
   ## Prepare the flow distance raster
   FlowLength = procGDB + os.sep + "overlandFlowLength"
   #[David - insert function call(s) to produc FlowLength here, and function definition(s) in Wtrshd_Functions script]
   
   ## Prepare the headwaters raster
   Headwaters = procGDB + os.sep + "Hdwtrs"
   makeHdwtrsIndicator(FlowLines, Catchments, clpShp, MaskNoWater, Headwaters)
   
   ## Calculate Overland Flow score ##
   FlowScore = procGDB + os.sep + "FlowScore"
   calcFlowScore(FlowLength, FlowScore, Headwaters)
   
   
   ### Karst procedures ###
   
   ## Prepare sinkholes
   # --> Manual operation: Prior to running density scoring function, make sure the sinkhole data are "clean", i.e., no overlaps/duplicates. There must also be a field representing the sinkhole area in desired units (i.e., square meters).
   
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
   
   
   ### Landscape Position Score ###
   PositionScore = procGDB + os.sep + "PositionScore"
   calcPositionScore(FlowScore, KarstScore, PositionScore)
   
   
   ### Potential Impact Score ###
   ImpactScore = procGDB + os.sep + "ImpactScore"
   calcImpactScore(PositionScore, SoilSensScore, ImpactScore)
   
   
   ### Finalize TIF Products ###
   arcpy.env.mask = jurisbnd
   fldr = r"N:\ProProjects\WatershedModels\CVWIM_Data\cvwimTIF"
   procList = [runoffScore, soilLossScore, soilSens_Score, FlowScore, KarstScore, PositionScore, ImpactScore]
   
   
if __name__ == "__main__":
   main()
   
