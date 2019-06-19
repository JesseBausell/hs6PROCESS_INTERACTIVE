% metaData_Reader_hs6
% Jesse Bausell
% June 6, 2019
%
% This matlab script will upload data from a metadata header hs6 file and
% assign variables to user-input headers as necessary. These metadata are
% used to facilitate the processing of raw hs6 data into Seabass
% submittable or Hydrolight compatible files.

[in_FILE, in_DIR] = uigetfile('.txt'); 
%prompt user to select a file
fid = fopen([in_DIR in_FILE]); 
%opens the file and provides a file identifyer (fid)
for i = 1:3
    % Removes the first three header lines
    fgetl(fid);
end
%% 1. Assign variables to metadata by indexing the = sign

% 1. in_FILE
linE = fgetl(fid);
eQ = regexpi(linE,'=');
in_FILE = linE(eQ+1:end); %replaces in_FILE with hs6 data (previously metadata file name)
% 2. in_DIR
linE = fgetl(fid);
eQ = regexpi(linE,'=');
in_DIR = linE(eQ+1:end); %replaces in_DIR with hs6 data pathway (previously metadata pathway)
% 3. affiliations
linE = fgetl(fid);
eQ = regexpi(linE,'=');
affiliations = linE(eQ+1:end);
% 4. investigators
linE = fgetl(fid);
eQ = regexpi(linE,'=');
investigators = linE(eQ+1:end);
% 5. contact
linE = fgetl(fid);
eQ = regexpi(linE,'=');
contact = linE(eQ+1:end);
% 6. experiment
linE = fgetl(fid);
eQ = regexpi(linE,'=');
experiment = linE(eQ+1:end);
% 7. station
linE = fgetl(fid);
eQ = regexpi(linE,'=');
station = linE(eQ+1:end);
% 8. latitude
linE = fgetl(fid);
eQ = regexpi(linE,'=');
lat = linE(eQ+1:end);
% 9. longitude
linE = fgetl(fid);
eQ = regexpi(linE,'=');
lon = linE(eQ+1:end);
% 10. document
linE = fgetl(fid);
eQ = regexpi(linE,'=');
doC = linE(eQ+1:end);
% 11. water depth
linE = fgetl(fid);
eQ = regexpi(linE,'=');
D = linE(eQ+1:end);
% 12. calibration files
linE = fgetl(fid);
eQ = regexpi(linE,'=');
cal_FILE_ac = linE(eQ+1:end);
% 13. date
linE = fgetl(fid);
eQ = regexpi(linE,'=');
dat = linE(eQ+1:end);

%% 2. Assign variables to lines 14-16

% 14. Kexp
linE = fgetl(fid);
eQ = regexpi(linE,'=');
kexp = str2num(linE(eQ+1:end));
% 15-16. Define binned absorption file names and their pathway
linE = fgetl(fid); % Extract line 14 from the header file
eQ = regexpi(linE,'='); % find = sign index (as before)
linE = linE(eQ+1:end);
if strcmpi(linE,'NA') % If no absorption files are listed,
    binned_ABS.FILES = ''; % Create structure for abs files, but leave it blank
    binned_ABS.PTHWY = '';
    fclose(fid); % close the header file.
else % If there are absorption files listed,
    linE = [',' linE ',']; % add extra comma to the end of the header line
    commA = regexpi(linE,','); % Find indices of all commas
    for ii = 1:length(commA)-1 
        % For-loop finds indices of binned absorption files. It then counts
        % the binned absorption files and places them into a structure. One
        % element per listed file.
        binned_ABS(ii).FILES = linE(commA(ii)+1:commA(ii+1)-1);
    end
    linE = fgetl(fid); % Extract the last line of the metadata file
    eQ = regexpi(linE,'=');
    % Places the pathway of the binned absorption files into binned_ABS
    % structure, first element only (see below)
    binned_ABS(1).PTHWY = linE(eQ+1:end); 
    fclose(fid); % Close the metadata file
end