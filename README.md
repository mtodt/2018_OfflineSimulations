# 2018_OfflineSimulations
Code to derive correction factors, implement correction factors in CLM4.5, post-process simulations, and create figures shown in the following publication:
Todt, M., Rutter, N., Fletcher, C. G., and Wake, L. M.: Simulated single-layer forest canopies delay Northern Hemisphere snowmelt, The Cryosphere Discuss., https://doi.org/10.5194/tc-2018-270, in review, 2019.

MATLAB code to derive correction factors and create Figure 3 in Todt et al. (2019) is given in CalculationOfCorrection_Todt2019_Cryosphere.m, which uses Multiple Linear Regression performed by MLR_SkySWin.m. Figure 2 in Todt et al. (2019) is created by RegressionGraph_Todt2019_Cryosphere.m. Creation of input data files MLRdata_XXXXX.mat is described in https://github.com/mtodt/2018_ToyModel.

Modified source code of CLM4.5, output configurations, and post-processing are included in two .zip-files, CORR.zip and CTRL.zip, which refer to the scenarios described in Todt et al. (2019). Note that post-processing files use paths and file names for scenarios, CLM45_NewControl_CRUNCEP_MR (CTRL) and CLM45_FinalLightbulb_CRUNCEP_MR (CORR), that need to be adjusted when recreating model runs. Model output was created a) hourly for years 2004 to 2007 in order to compare simulations and measurements for the forest stand at Alptal,  Switzerland, which was also used to create Figure 6; and b) daily for snow seasons 1981/82 to 2015/16 in order to assess impact on longwave enhancement and snow cover across the Northern Hemisphere.

MATLAB code to create Figures 4 and 5 shown in Todt et al. (2019) is given in Analysis_Alptal_Todt2019_Cryosphere.m

MATLAB code to create Figures 1, 6, 7, and 8 shown in Todt et al. (2019) is given in Analysis_Global_Todt2019_Cryosphere.m. This necessitates surfdata_0.9x1.25_simyr2000_c130418.nc, the surface data set used by CLM4.5, which can be downloaded at https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/surfdata_map/. Maps featuring polar projections are created using mapping toolbox M_Map available at https://www.eoas.ubc.ca/~rich/map.html.
