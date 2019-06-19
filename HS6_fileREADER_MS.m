% HS6_fileREADER_MS
% Jesse Bausell
% September 15, 2017
%
% Updated: June 6, 2019
% 
% This script is designed to run with GSR_hs6PROCESS. It reads a raw hs6
% file into matlab and delineates important variables that will be used in
% hs6 processing. It also performs pure-water scattering correction on
% raw data to convert coefficients from backscattering to particulate 
% backscattering.

%% 1. Extract data from user-supplied header file

metaData_Reader_hs6; % Processes the metadata header file required to run 
% hs6PROCESS_SEABASS.

%% 2. Extract data from hs-6 raw data file (.dat) output by HydroSoft

fid = fopen([in_DIR in_FILE]); % Opens raw hs6 file and assigns it a file identifier
bb0 = []; %Empty variable that will hold the pure water instrument calibration value.
while 1
    % This while loop cycles through all of the header files in the raw hs6
    % file that is output by Hydrosoft (HOBI Lab software package). The
    % loop will find necessary information contained inside the file
    % header.  
    linE = fgetl(fid); % Find each header line, one per loop iteration
    % Search for pure-water instrument calibration value and variable
    % header (below)
    if ~isempty(regexpi(linE,'bb0='))
        % If pure-water instrument calibration value is located.
        bb0_IND = regexpi(linE,'='); % Index the numerical pure water value
        bb0 = str2double(linE(bb0_IND+1:end)); % Extract numerical value and convert to float
    elseif ~isempty(regexpi(linE,'Time'))
        % If variable header line is located       
        fgetl(fid); % Skip one line (there is a space between the data).
        linE = [',' linE ',']; %book end the .dat field headers line with commas
        txtscn_fodder = '%n %f'; % beginning of format specifiers string
        titlE = ',Time,Depth,'; % beginning of the new header string
        % Save the previous line and put a comma ahead and behind it.
        commA = regexpi(linE,','); % Find all commas in the header line

        for ii = 1:length(commA)-1
            % This for-loop will creates the format specification string
            % for textscan (reading raw hs6 data into matlab). It also
            % creates a new header with hs6 fields that were read into
            % matlab (the program doesn't read in everything). This header
            % will later be used to order data columns (form desired
            % fields) into the hs6 backscatter and fluorescence data matrices.
            headER = linE(commA(ii)+1:commA(ii+1)-1); % Extract one field header at a time
            if ~isempty(regexpi(headER,'bb\d{3}uncorr')) && ... % if header contains "bb" and "uncorr"
                    isempty(regexpi(headER,'beta')) % and the header doesn't contain "beta"
                txtscn_fodder = [txtscn_fodder '%f']; % add float format specifier
                titlE = [titlE headER ',']; % Also add the header onto the field header onto new header string 
            elseif ~isempty(regexpi(headER,'[fl]\d{3}')) && ...  % if field header contains "fl"
                    isempty(regexpi(headER,'uncorr')) && ...  % and it doesn't contain "uncorr"
                    isempty(regexpi(headER,'beta')) % and it doesn't contain "beta"
                txtscn_fodder = [txtscn_fodder '%f']; % add float format specifier
                titlE = [titlE headER ',']; % add field header onto new header string
            elseif isempty(regexpi(headER,'TIME')) && isempty(regexpi(headER,'DEPTH')) % any other scenario
                txtscn_fodder = [txtscn_fodder '%*f']; % Add a format specifier for unwanted float
            end
        end
            % Before while-loop ends, read hs6 data into matlab
            Data_Grid = textscan(fid,[txtscn_fodder ' %*[^\n]'],'Delimiter',',');
            fclose(fid); % close the raw data file to prevent its corruption
            break % Break the while loop.    
    end    
end

%% 3. Formulate backscattering and fluorescence data matrices

% Using the newely constructed header, titlE ....
commA_title = regexpi(titlE,','); % Find all comma indices
uncorr_IND = regexpi(titlE,'bb...uncorr'); %  find indices of backscatter headers
uncorr_IND_fl = regexpi(titlE,'fl'); % find indices of fluorescence

% Create empty nan matrices in which to place hs6 raw data and wavelengths
HS6_data = nan(length(Data_Grid{1}),length(uncorr_IND)); % backscatter data
HS6_data_fl = nan(length(Data_Grid{1}),length(uncorr_IND_fl)); % fluorescence data
lambdA = nan(1,length(uncorr_IND)); % backscattering wavelengths
lambdA_fl = nan(1,length(uncorr_IND_fl)); %fluorescence wavelengths


for ii = 1:length(lambdA) 
    % This for-loop will take empty backscattering variables (HS6_data and
    % lambdA) and fill them with appropriate data. It will do this one
    % wavelength (column) at a time.
    
    % find column index for bb wavelength
    poinT = abs(commA_title-uncorr_IND(ii));     
    placE = find(poinT == min(poinT)); 
    % use the column index to extract data from Data_Grid (the hs6 cell array
    % created by textscan) and place it into the newly created nan matrix.
    % Do the same for the wavelength (below).
    HS6_data(:,ii) = Data_Grid{placE}; % backscatter
    waveSTRING = titlE(commA_title(placE):commA_title(placE+1)); % header string    
    wave_IND = regexpi(waveSTRING,'\d'); % wavelength index within header string
    lambdA(ii) = str2double(waveSTRING(wave_IND)); % wavelength as a number        
end 


for ii = 1:length(lambdA_fl)
    % This for-loop will take empty fluorescence variables (HS6_data_fl and
    % lambdA_fl) and fill them with appropriate data. It will do this one
    % wavelength (column) at a time.

    % find column index for fl wavelength
    poinT = abs(commA_title-uncorr_IND_fl(ii));     
    placE = find(poinT == min(poinT));
    % use the column index to extract data from Data_Grid (the hs6 cell array
    % created by textscan) and place it into the newly created nan matrix.
    % Do the same for the wavelength (below).
    HS6_data_fl(:,ii) = Data_Grid{placE}; % fluorescence    
    waveSTRING = titlE(commA_title(placE):commA_title(placE+1)); % header string   
    wave_IND = regexpi(waveSTRING,'\d'); % wavelength index within header string
    lambdA_fl(ii) = str2double(waveSTRING(wave_IND)); % wavelength as a number
end

% Create reference variables: Depth and time
deptH = Data_Grid{2}; % depth
timE = Data_Grid{1}; % time

%% 4. Calculate and subtract pure-water backscatter to create particulate backscatter

% Convert backscattering coefficients to particulate backscattering using 
% Morel, A. (1974)
HS6_data = morel_Correction_pw(HS6_data,lambdA,bb0); % particulate backscattering