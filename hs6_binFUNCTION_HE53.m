function hs6_binFUNCTION_HE53(filE_acs,filE_hs6,in_FILE,kexp)
% hs6_binFUNCTION_HE53
% Jesse Bausell
% June 9, 2019
%
%
% This function is designed to bin HS6 data and perform sigma-correction on
% particulate backscattering data using methods described in Doxaran et al.
% (2016). The function uses binned absorption data for the
% sigma-correction. This matlab script is written to as a nested function
% for hs6PROCESS_INTERACTIVE, however, it can also opperate by itself.
%
% Doxaran, D., E. Leymarie, B. Nechad, et al. (2016) Improved correction
% methods for field measurements of bbp in turbid waters. Optics Express
% 24(4): 3615 - 3637.
%
% Inputs:
% HS6_mat - unbinned HS6 data incudes depth, bbp, and flourescence
% filE_acs - file pathway to ac-s data for use in Doxaran prescribed
% sigma-correction. These files must be binned and formatted for hydrolight
% filE_hs6 - file pathway to recently processed hs6 data, used to collect
% file header.
% kexp - wavelength-dependent kexp coefficients used in sigma-correction.
% These change with individual instruments.
%
% Outputs:
% Binned Seabass-compatible hs6 data file.
%
% The function assumes that input file will be formatted up as follows:
%
%
% hs6 (individual readings)
% /begin_header
% ...
% /fields=time,depth,bbp###, ...,stimf_ex###_em###,stimf_ex###_em###
% ...
% /end_header
% -DATA-
%
% ac-s (depth-binned)
% /begin_header
% ...
% /fields=depth,agp###, ...,agp###_SD....
% ...
% /end_header
% -DATA-
%% 1. Read in HS6 data file into matlab and organize data variables

load(filE_hs6); % Load hs6 data from .mat file

% Part a. Order data matrices by depth (vertical) 
[deptH, d_order] = sort(deptH); % Order depths by ascending values with indices
HS6_data_bbp = HS6_data_bbp(d_order,:); % Use ascending depth indices (d_order) to re-order bbp matrix
HS6_data_fl = HS6_data_fl(d_order,:); % Use ascending depth indices (d_order) to re-order fl matrix

% Part b. Order data matrices by wavelength (horizontal) 
[lambdA_bbp, l_ordeR] = sort(lambdA_bbp); % Order bbp wavelengths by ascending values with indices
HS6_data_bbp = HS6_data_bbp(:,l_ordeR); % Use ascending wavelength indices (l_ordeR) to re-order fl matrix
kexp = kexp(l_ordeR); % sort kexp values so that they remain in the same order as bbp wavelengths

%% 2. Read in absorption data for binning.
% This section of code takes the input binned absorption data and reads it
% into matlab. This data file should be binned absorption formatted to be
% Seabass-comaptible.

fid_acs = fopen(filE_acs); 
%opens the file and provides a file identifyer (fid)
txtscn_fodder = '%f';

while 1
    % This while-loop opens the binned absorption file and creates all
    % necessary variables that are used in the subsequent parts of this
    % script. It like the while-loop in section 1, it catalogues the file 
    % header one line at a time, and uses textscan to read in absorption
    % data.
    linE = fgetl(fid_acs); % Examines one line of binned absorption file header
    if ~isempty(regexpi(linE,'Depth'))
        % If the file header line is the line directly above absorption
        % data, read in absorption data and close the binned absorption
        % file. Find the 10th line that says "Depth."
        hdr_linE = linE; % Change variable name of header line
        LTR_ind = regexpi(hdr_linE,'[a-z]'); % Index all non-numers in the line
        hdr_linE(LTR_ind) = []; % Eliminate all non-numbers, leaving a string of numbers
        lambdA_abs = str2num(hdr_linE); % Convert string into array with wavelengths
        l_wv = length(lambdA_abs)/2; % Divide array by 2, wavelengths are repeated
        lambdA_abs = lambdA_abs(1:l_wv); % Make sure to take just one set of wavelengths (get rid of repitition)

        for hh = 1:2
            % This for-loop creates the format specifications for textscan.
            % It goes through two iterations, the first one lays down float
            % specifications to read in .txt data. The second uses '*f'
            % specifications, which instruct textscan NOT to read in the
            % data. *f are for attenuation fields.
            for ii = 1:l_wv
                if isequal(hh,1)
                    txtscn_fodder = [txtscn_fodder '%f'];
                else
                    txtscn_fodder = [txtscn_fodder '%*f'];
                end
            end
        end
        abs_MATRIX = textscan(fid_acs,[txtscn_fodder '\n'],'Delimiter','\t'); % Read data into cell array
        fclose(fid_acs); % Close binned absorption data file
        deptH_abs = abs_MATRIX{1}(2:end); % Array of binned depths
        abs_MATRIX = cell2mat(abs_MATRIX); % Convert absorption cell array to matrix
        abs_MATRIX = abs_MATRIX(2:end,2:l_wv+1); % Cut off depth, leaving absorption values
        break % break the while loop
    end
end


%% 3. Perform sigma-correction on bbp using binned absorptions
% Using the binned absorption data read into matlab (see above section),
% hs6_binFUNCTION now interpolates absorption values according to bbp
% wavelengths. It then procedes to sigma-correct bbp data.

A_particulate = nan(length(deptH_abs),length(lambdA_bbp));  
% Creates a nan matrix to hold interpolated absorption values corresponding
% with bbp wavelengths

for ii = 1:length(deptH_abs)
    % Interpolates absorption spectrum for HS6 wavelengths.
    A_particulate(ii,:) = lambda_INTERPOLATE(lambdA_abs,abs_MATRIX(ii,:),lambdA_bbp);
end
    
for ii = 1:length(deptH)
    % This loop performes the Doxaran sigma-correction one HS6 line at
    % a time. It matches the HS6 line with an absorption line of
    % similar (or closest) depth, and then performs all steps of
    % sigma-correction. 
    HS6_data_bbp(ii,:) = Doxarian_SIGMA(HS6_data_bbp(ii,:),deptH(ii),kexp,A_particulate,deptH_abs);
end

%% 4. Bin data according to depth bin of absorption data

depth_BINSIZE = min(abs(deptH_abs(1:end-1)-deptH_abs(2:end))); 
% Find the depth bin size for binned absorption data

max_DEPTH = ceil(max(deptH)); % maximum depth of hs6. Round up to the nearest whole number
if ~isequal(rem(max_DEPTH,2),0)
    % If the max_DEPTH is an odd number
    max_DEPTH = max_DEPTH+1; 
    %add one to the already rounded up maximum depth to make it an even
    %number
end

% Find the dimensions of the bbp data matrix
[v,h] = size(HS6_data_bbp); 
% create empty nan matrix to hold binned bbp and fl means and standard
% deviations
%BINNED_HS6 = nan(round(max_DEPTH/depth_BINSIZE),1+2*h+2*length(lambdA_fl));  
BINNED_HS6 = nan(round(max_DEPTH/depth_BINSIZE),1+h);  

for jj = 0:depth_BINSIZE:max_DEPTH-depth_BINSIZE
    % This for-loop cycles through the hs6 data (bbp & fl) and calculates
    % means and standard deviations for each depth bin using nearest
    % neighbor approach. The loop creates binned values one depth-bin at a
    % time.
    
    % Calculate depth bin and place it in the binned hs6 matrix
    vert_IND = length(0:depth_BINSIZE:jj); % binned depth index
    BINNED_HS6(vert_IND,1) = (jj+jj+depth_BINSIZE)/2; % bined depth value
    
    % Calculate average of bbp and fl spectra calculated across depth bin
    BIN_IND = find(deptH >=jj & deptH <jj+depth_BINSIZE); % index for min and max depth-bin values
    BINNED_HS6(vert_IND,2:h+1) = nanmean(HS6_data_bbp(BIN_IND,:),1); % mean bbp spectrum for depth-bin
    %BINNED_HS6(vert_IND,h+2:h+3) = nanmean(HS6_data_fl(BIN_IND,:),1); % mean fl spectrum for depth-bin
    
    % Calculate standard deviation of bbp and fl spectra calculated across depth bin
    %BINNED_HS6(vert_IND,1+h+length(lambdA_fl)+1:1+2*h+length(lambdA_fl)) = nanstd(HS6_data_bbp(BIN_IND,:),0);
    %BINNED_HS6(vert_IND,1+2*h+length(lambdA_fl)+1:end) = nanstd(HS6_data_fl(BIN_IND,:),0);
end

nNAN_IND = find(~isnan(BINNED_HS6(:,2))); % Find index of all non NAN rows
BINNED_HS6 = BINNED_HS6(nNAN_IND,:); % Eliminate Nan rows from hs6 data matrix

%% 5. Create Hydrolight-compatable binned hs6 .txt file.

fid_HS6 = fopen([filE_hs6(1:end-4) '_bin' num2str(depth_BINSIZE) '_HE53.txt'],'w');
% Create a new file in which to place newly-created header and binned,
% sigma-corrected data.
    % Header lines 1-9. Hydrolight requires 10 header rows.
    fprintf(fid_HS6,'%s\n','Particulate Back Scattering Coefficients: Binned Sigma-Corrected');
    fprintf(fid_HS6,'%s\n',['Corrected: ' datestr(clock)]);
    fprintf(fid_HS6,'%s\n',['Instrument: Hobi Labs Hydroscat-' num2str(length(lambdA_bbp))]);
    fprintf(fid_HS6,'%s\n',['File name: ' in_FILE]);
    fprintf(fid_HS6,'%s\n','Doxaran Sigma-correction: yes');
    IND_acs = sort([regexpi(filE_acs,'/') regexpi(filE_acs,'\')]);
    fprintf(fid_HS6,'%s\n',['ac-s file name: ' filE_acs(IND_acs(end)+1:end)]);
    fprintf(fid_HS6,'%s\n','C-interpolated: yes');
    fprintf(fid_HS6,'%s\n',['bin=' num2str(depth_BINSIZE)]);
    fprintf(fid_HS6,'%s\n','Data have been processed using code written and made avaiable by Jesse Bausell (email: jbausell@ucsc.edu, GitHub: JesseBausell).');

for hh = 1:2  
    % This for-loop is responsible for printing lines 10 and 11 onto the
    % .txt file.
    if isequal(hh,1)
        % For header line 10
        fprintf(fid_HS6,'%s\t','Depth'); % start header line with "Depth"
        frnt = 'bb'; % All subsequent headers have bb in front of them
    else
        % For header line 11
        fprintf(fid_HS6,'%s\t',num2str(length(lambdA_bbp))); % start header line with number of bb channels
        frnt = ''; % Do not put anything infront of subsequent headers
    end
    FMT_SPEC = '%4.2f\'; % Beginning of string of format specifier string
    for ii = lambdA_bbp
        % This for-loop determines number of float ("f") format specifier
        % based on the number of numerical fields. It determines this as double
        % the number of wavelengths
        FMT_SPEC = [FMT_SPEC 't%10.8f\']; % Append new format specifiers 
            if isequal(ii,lambdA_bbp(end))
                % At the end of the lines 10 & 11
                fprintf(fid_HS6,'%s\n',[frnt num2str(ii)]); % Print wavelength with end of line character
            else
                fprintf(fid_HS6,'%s\t',[frnt num2str(ii)]); % Print wavelength with tab delimiter
            end               
    end
end
    
for ii = 1:length(BINNED_HS6)    
    % Print data into binned hs6 .txt file one line at a time
    fprintf(fid_HS6,[FMT_SPEC 'n'],BINNED_HS6(ii,:));
end
fclose(fid_HS6); % Close the file to cease editing.