function bb_MAT_CORR = morel_Correction_pw(bb_MAT,bb_lambdA,bb0)
% morel_Correction_pw
% Jesse Bausell
% June 7, 2019
% 
% This function goes through backscattering data one wavelength (column) at a
% time and calculates the wavelength-dependent pure-water
% backscattering coefficient. It uses methods from Morel (1974)
% 
% Morel, A. (1974) "Optical Properties of Pure Water and Pure Sea Water".
% In: Optical Aspects of Oceanography. M. G. Jerlov and E. S. Nielsen
% (eds.), Academic Press, New York, 1974, chap. 1, pp 1-24.
% 
% Inputs: 
% bb_MAT - uncorrected backscattering coefficients. Oriented be depth
% (rows) and wavelengths (columns)
% bb_lambdA - wavelengths. Array must be horizontal and oriented the same
% as bb_MAT columns.
% bb0 - instrument-specific pure-water constant (extracted from
% HydroSoft-generated .dat file header)
%
% Outputs:
% bb_MAT_CORR - particulate backscattering coefficients (pure-water
% backscattering subtracted).

for ii = bb_lambdA    
    % This for-loop goes through backscattering data one wavelength at a
    % time and calculates the wavelength-dependent pure-water
    % backscattering coefficient using methods from Morel (1974). See
    % readme for more information. 
    
    bw = bb0*(525/ii)^4.32; % wavelength dependednt pure-water backscattering coefficient 
    hor_IND = (bb_lambdA==ii); % find wavelength column index for HS6_data
    
    % subtract pure-water backscattering coefficient from backscattering
    % matrix (one wavelength/column at a time)
    bb_MAT(:,hor_IND) = bb_MAT(:,hor_IND) - bw; 
end

bb_MAT_CORR = bb_MAT; % Rename (corrected) particulate backscattering variable
