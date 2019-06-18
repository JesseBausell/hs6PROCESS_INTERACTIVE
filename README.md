# hs6PROCESS_INTERACTIVE
Consistent with up to date protocols , hs6PROCESS_INTERACTIVE processes raw backscattering data as sampled in natural water bodies using HOBI Labs Hydroscat-6 (hs6). hs6PROCESS_SEABASS should is also designed to be compatible with Hydroscat-4 and 2.

This Matlab script processes raw backscattering coefficients, as sampled in natural water bodies using HOBI Labs Hydroscat-6 Backscattering Meter and Fluorometer (hs6). Consistent with up to date protocols, hs6PROCESS_INTERACTIVE processes raw backscattering data (engineering units) sampled in natural water bodies using hs6. hs6PROCESS_INTERACTIVE is designed to be compatible with Hydroscat-4 (hs4) and Hydroscat-2 (hs2), however it remains untested with these two instrument variants. We also point out that, the expectation for two fluorescence (fl) channels is hardcoded into the hs6PROCESS_SEABASS script. Thus (uncorrected) input hs6 data MUST contain two fl channels. Otherwise it may not work properly. The wavelengths of these hardcoded  fl channels are as follows: 470 excitation/510 emission (channel 1) and 442 excitation/700 emission (channel 2). Nevertheless, if these wavelength values are incorrect, user should not worry. hs6PROCESS_SEABASS does not perform any processing steps on fl values that require knowing the correct fl wavelengths. This means that user can simply "correct" fl wavelengths in the output files after processing.

hs6PROCESS_INTERACTIVE uses up to date processing protocols (as of June 2019) to process backscattering coefficients, which it outputs as both individual as well as depth-binned spectra. All output data products are formatted as to be consistent with Hydrolight 5.0. User should run hs6PROCESS_INTERACTIVE AFTER running acs_PROCESS_INTERACTIVE, if he/she wishes to perform Doxaran et al. (2016) sigma-correction.

Inputs:
metadata_HeaderFile_hs6.txt - ascii file containing metadata required to process hs6 data 

Outputs:
Station_#_bbp.mat - MAT/HDF5 file containing individual particulate backscattering (bbp) + fl spectra. No sigma correction applied.
Station_#_bb_bin#.txt* - Seabass-formatted ascii file(s) containing sigma-corrected and depth-binned bbp spectra, as well as depthb-binned fl spectra

Required Matlab Scripts and Functions:
Doxarian_SIGMA.m
hs6_binFUNCTION.m
HS6_fileREADER_MS.m
lambda_INTERPOLATE.m
metaData_Reader_hs6.m
morel_Correction_pw.m

Required data files:
Seabass_header_ACS5.mat

Program Description:
hs6PROCESS_SEABASS processes raw field-collected backscattering coefficients following a series of steps. It is outfitted to process raw data contained in Hydrosoft-output ascii (.dat) files. Hydrosoft is a free software package provided by HOBI Labs for preliminary processing of its optical oceanographic instrumentation. 
  1. Reads ascii data. Accepts uncorrected backscattering coefficients and uncorrected fluorescence.
  2. Calculates particulate backscattering coefficients (bbp)
    a. Calculates wavelength-dependent backscattering coefficients of pure-water using methods devised by Morrel (1974)
    b. Subtracts pure-water backscattering coefficients from total (uncorrected) backscattering coefficients
  3. Time stamps bbp and fl spectra. It can also synchronize bbp spectra with ac-s data processed by acsPROCESS_SEABASS 
  4. QA/QC hs6 data. Individual bbp spectra + fluorescence measurements are flagged for removal if bbp are less than zero or greater than      0.5 /m
  5. Produces SeaBASS-formatted ascii (.txt) file containing time-stamped bbp spectra/fluorescence values with depths at which they were        sampled
  6. Depth-bin bbp and fl spectra.
    a. Sigma-correct bbp spectra according to Doxaran et al. (2016) using depth-binned absorption spectra measured with ac-s. A binned 
    absorption spectrum is chosen for each bbp spectrum using nearest neighbor approach. 
    b. bbp and fl are depth-binned at using the same bin size as the absorption spectra used to sigma-correct them
  7. Produces SeaBASS-formatted ascii (.txt) file(s) containing depth-binned bbp and fl average spectra and standard deviations. 

User Instructions:
  1. Fill out metadata_HeaderFile_hs6.txt (as specified below)
  2. Run hs6PROCESS_SEABASS using Matlab command window.
  3. Select appropriate metadata_HeaderFile_hs6.txt file when prompted. 
  6. Time-stamp hs6 data for potential syncrhonization with ac-s (if ac-s was deployed simultaneously)
    a. Matlab will produce a "time series" plot indicating hs6 position (depth) over "time" (spectum index).
    b. User selects a reference point on this plot by entering the index of the desired point into the command window. Data cursor is a 
    useful resource.
    c. User is asked to confirm his/her selection on command window with y/n keys. If user rejects his/her selection, he/she will be           prompted to try again.
    d. User is now asked to enter a "reference time". This will be used to determine the exact time of day that the reference point (step
    6b) was measured by hs6. This time of day should be GMT and formatted as military time (HH:MM:SS). This reference time is used to 
    back and forward-calculate time of day from the reference point.***
    
 ***If user wishes to syncrhonize hs6 time stamps with those of an ac-s cast, he/she should select a reference point in the same position as that of ac-s cast adn enter the time stamp appearing at the top of the figure (see Station_#_ACS.fig) when prompted to do so. As with acsPROCESS_SEABASS, this time should be GMT and formatted at military time.
 
Filling out metadata_HeaderFile_hs6.txt:
hs6PROCESS_SEABASS relies on a metadata header to process hs6 data. All information pertaining to the specific hs6 cast should be included in this header. A header template (metadata_HeaderFile_hs6.txt), indicating important fields, is provided in the hs6PROCESS_SEABASS repository on GitHub. When filling out this header file, the first three headers (indicating user instructions) should be left alone. Required information fields contain = signs. USER SHOULD ONLY ALTER TEXT APPEARING ON THE RIGHT HAND SIDE OF =. User should indicate unavailability of desired information with "NA". DO NOT DELETE ROWS! Below are fields contained in metadata_HeaderFile_hs6.txt and instructions on how to fill them out. Spaces should never be used in header fields; use underscore instead (_).

data_file_name=indicate name of ascii (.dat) file containing unprocessed hs6 data. This file is generated using HydroSoft software created by HOBI Labs to support their instrument platforms (e.g. HydroScat, a-Beta, c-Beta & Abyss-2) by converting raw signals to engineering units. Please note: When using HydroSoft to output .dat files, user should be sure NOT to apply pure-water correction; hs6PROCESS_SEABASS applies this step already. It is ok for user to output both sigma-corrected and uncorrected backscttering coefficients as hs6PROCESS_SEABASS can differentiate between them.

data_file_name=pathway for aforementioned HydroSoft-generated ac-s .dat file (data_file_name). This pathway should include the folder in which sits, and should be ended using "/" or "\" for mac and pc respectively. 

affiliations=name of company or research institution with which investigators are affiliated. 

investigators=lists of investigators. Multiple names should be separated by commas and "_" should be used in place of spaces.

contact=email of principle investigator

experiment=name of experiment or field campaign 

station=field station number 

latitude=latitude of field station. This should be indicated in decimal form. DO NOT format in minutes or seconds. Do not include Roman letters. South should be indicated with a negative sign.

longitude=longitude of field station. This should be indicated in decimal form. DO NOT format in minutes or seconds. Do not include Roman letters. West should be indicated with a negative sign.

documents=additional documents user wishes to submit to SeaBASS. DO NOT INDICATE kudelalab_HS6_readme.pdf. This is printed automatically in output files.

water_depth=bottom depth of the field station in meters. Numerals only. Do not include units.

calibration_file=name of original factory-supplied HOBI Labs calibration file. This file contains instrument-specific coefficients used to convert raw signals (measured by hs6) into engineering units. Although this file is not used by hs6PROCESS_SEABASS to process hs6 data (it was used in Hydrosoft), SEABaSS encourages submitters to disclose all ancillary calibration files used in data processing. hs6PROCESS_SEABASS will correctly place this file into headers of output files (e.g. Station_#_bb.txt and Station_#_bb_bin#.txt).

date(yyyymmdd)=indicate date on which ha6 was deployed.

kexp_vals=wavelength-dependent constants (kexp) used to calculate sigma coefficients of attenutation, Kbb, (see HydroSoft user manual or kudelalab_HS6_readme.pdf for details). These coefficients are instrument-specific and can be found inside the factory-supplied calibration file (header field=calibration_file). They should be listed in the same order as wavelengths of backscattering coefficients appear in the header of the HydroSoft-generated .dat file containing the unprocessed hs6 data (header field=data_file_name). In the event that instrument-specific Kexp_vals are unavailable for a particular Hydroscat, the constant 0.14 can be used for all wavelengths. In the event that it's used, this value should still be listed one time for each backscattering channel (e.g. 6 times for hs6).

apg_bin_files=binned ac-s measured absorption files with which to sigma-correct bbp spectra. These files must be formatted for SEABaSS submission and would ideally contain data measured at the same time and location as hs6. DO NOT specify depth bin sizes in the metadata header file. hs6PROCESS_SEABASS bins hs6 data according to the bin sizes used in absorpiton files listed in hte header field. To prevent hs6PROCESS_SEABASS from depth-binning bbp and fl data, place "NA" after the = sign (e.g. apg_bin_files=NA).

apg_bin_path=pathway for apg_bin_files. This pathway must be the same for all binned absorption files should more than one be listed.

Metadata Header File Example:
hs-6 metadata template
Template contains information necessary for the processing of hs-6 data files (.dat) output using the HOBI Labs software program HydroSoft. Use commas to separate names of investigators and files, but DO NOT leave ANY spaces between words. If a space is unavoidable, use an underscore between words (like_this). Unknown or unavailable information should be indicated with NA. Latitude and longitude should be in decimal degrees and water depth should be in meters. Do not include units of measurement. These will be added later by the program. 
#### DO NOT ALTER HEADER FIELDS####
data_file_name=COAST18.dat
data_file_path=/Users/JBausell/Documents/acs_data/
affiliations=UC_Santa_Cruz
investigators=Jesse_T_Bausell,_Easter_B_Bunny,Kris_B_Kringle
contact=kudela@ucsc.edu
experiment=COAST
station=18
latitude=36.8972
longitude=-121.8859
documents=NA
water_depth=24
calibration_files=HS990216_v3.ca 
date(yyyymmdd)=20181025
kexp_vals=[0.1450003,0.14499998,0.14200014,0.1439971,0.145993,0.147154]
apg_bin_files=COAST_18_ac-s_bin_0.5.txt,COAST_18_ac-s_bin_1.txt,COAST_18_ac-s_bin_2.txt
apg_bin_path=/Users/JBausell/Documents/acs_data/

Bibliography:

Doxaran, D., E. Leymarie, B. Nechad, A. Dogliotti, K. Ruddick, P. Gernex, and E. Knaeps, Improved correction methods for field measurements of particulate light backscattering in turbid waters. Optics express, 2016. 24: p. 3615-3637.

Morel, A, Optical Aspects of Oceanography. in Optical Aspects of Oceanography, M. G.
Jerlov, and E. S. Nielsen, eds. (Academic Press Inc, 1977).
