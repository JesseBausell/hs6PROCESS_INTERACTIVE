% hs6PROCESS_INTERACTIVE
% Jesse Bausell
% September 13, 2017
%
% Updated: June 11, 2019
%
% This matlab script takes raw hs6 data files that are produced by
% HydroSoft, a commercial software package produced by HOBI Labs, parent
% company of hs6. hs6PROCESS_SEABASS re-organizes data into a format that
% is compatible with NASA's SEABASS data portal. Data, should NOT be
% sigma-corrected in HydroSoft, as this script performs sigma-correction if
% indicated by the metadata file (see readme). The script outputs files
% into Hydrolight (version 5) format and allows user to choose a subset of
% the multicast, as well as to interactively QA/QC spectral bbp.
%
% Required matlab scripts and functions:
% dataSNATCH5.m
% depthSELECTOR.m
% Doxarian_SIGMA.m
% hs6_binFUNCTION_HE53.m
% HS6_dataFLAGGER5.m
% HS6_fileREADER_MS.m
% lambda_INTERPOLATE.m
% metaData_Reader_hs6.m
% morel_Correction_pw.m
%
% Required data files:
% Seabass_header_ACS5.mat
%
% Input:
% metadata_headerfile.txt (see readme)

clear all; close all; clc;

%% 1. Read in HS6 data
% Read in the selected HS6 file and creates the variables and create variables
% to be used in the rest of the program.
HS6_fileREADER_MS; % Script to uipload data and create variables associated 
% with particulate backscattering.

% Before any subsequent steps, user will select the subsets of the hs6
% multicast that he/she wants to process (below). 
[deptH_IND, deptH] = depthSELECTOR(deptH); % subset cast depths & indices

% Index matrices and arrays associated with depth
HS6_data = HS6_data(deptH_IND,:); % depth-subset of backscattering
HS6_data_fl = HS6_data_fl(deptH_IND,:); % depth-subset of fluorescence
timE = timE(deptH_IND); % depth-subset of time (decimal date)

%% 2. QA/QC HS6 data.

[deptH, deptH_IND] = sort(deptH); % sort depth in ascending order
HS6_data = HS6_data(deptH_IND,:); % use depth order indices to sort bbp
HS6_data_fl = HS6_data_fl(deptH_IND,:); % use depth order indices to sort fl
timE = timE(deptH_IND); % use depth order indices to sort time

[lambdA, lamda_IND] = sort(lambdA); % sort wavelength in ascending order
HS6_data = HS6_data(:,lamda_IND); % use wavelength order indices to sort bbp
%HS6_data_fl = HS6_data_fl(:,lamda_IND); % use wavelength order indices to sort fl
kexp = kexp(lamda_IND);

flaggedROWS = HS6_dataFLAGGER5(HS6_data,lambdA,NaN);
% Allows user to look at HS6 depth profiles and flag data that he/she
% finds to be questionable. Data will be reviewed in step 4.

max_2BIN = ceil(max(deptH)); % Create a max integer for 2 m depth bin
if rem(max_2BIN,2) > 0
    % If max integer is odd, add 1 to make it even
    max_2BIN = max_2BIN + 1;
end

flaggedROWS = dataSNATCH5(HS6_data,deptH,flaggedROWS,lambdA,{0:2:max_2BIN},'bbp',1); 

if isempty(flaggedROWS)
    HS6_data_fl(flaggedROWS,:) = []; % fl
    HS6_data(flaggedROWS,:) = []; % bbp
    deptH(flaggedROWS) = []; % depth
    timE(flaggedROWS) = []; % time
end
% Rename HS6_data and lambdA variables to avoid variable confusion when re-opening .mat/hdf5 file
HS6_data_bbp = HS6_data;
lambdA_bbp = lambdA;

save([in_DIR experiment '_' station  '_bbp'],'HS6_data_bbp','HS6_data_fl','deptH','timE','lambdA_bbp','lambdA_fl','-v7.3');
%BIN_acDATA_HE53; % Bin ac-s data and print it in Hydrolight-compatible .txt files


%% 2. Sigma-correct bbp and produce binned files
% If user specifies binned absorption files in the metadata header, this
% section will bin and sigma-correct bbp data according to Doxaran et al.
% (2016).

if ~isempty(binned_ABS(1).FILES)
    % If user selects to run Doxaran sigma-correction by listing binned
    % ac-s absorption files in metadata header.
    diR_acs = binned_ABS(1).PTHWY;
    for ii = 1:length(binned_ABS)
        % This for-loop cycles through binned absorption files listed in
        % the metadata header. It performs a Doxaran sigma-correction for
        % each of them.
        filE_acs = [diR_acs binned_ABS(ii).FILES]; %Give ac-s file and pathway it's own variable
        filE_hs6 = [in_DIR experiment '_' station  '_bbp.mat']; % Give hs6 file and pathway it's own variable
        hs6_binFUNCTION_HE53(filE_acs,filE_hs6,in_FILE,kexp); % Create binned hs6 file with sigma-corrected bbp
    end
end
        
