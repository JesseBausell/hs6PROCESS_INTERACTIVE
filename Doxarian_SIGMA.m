function [sigmaDOX_bb] = Doxarian_SIGMA(HS6_uncorrected,HS6_deptH,kexp,A_particulate,A_particulate_DeptH)
% Doxarian_SIGMA
% Jesse Bausell
% November 15, 2016
%
%
% This matlab function performs a sigma-correction on hs6-measured
% particulate backscattering coefficients using methods described in
% Doxaran et al. (2016), namely using field-collected particulate/Gelbstoff
% absorption coefficients.
%
% Inputs:
% HS6_uncorrected - non-corrected bbp spectrum. 
% HS6_deptH - depth value assigned to bbp spectrum
% kexp - instrument-specific attenuation coefficient (used to calculate
% sigma)
% A_particulate - matrix of binned absorption coefficients interpolated to bbp wavelenths
% A_particulate_DeptH - absorption coefficient depth bins
% 
% Outputs:
% sigmaDOX_bb - sigma-corrected bbp matrix

%% 1. Perform several checks to make sure that the input variables are correct.

[aHS6, bHS6] = size(HS6_uncorrected); % Vertical and horizontal dimentions of bbp matrix
[aACS, bACS] = size(A_particulate); % Vertical and horizontal dimentions of binned agp matrix

    if ~isequal(bACS,bHS6)
        % bbp and binne absorption have different numbers of wavelengths
        % (columns)
        error('Make sure that ac-s and bbp have the same wavelengths'); 
    elseif ~isequal(length(HS6_deptH),1) || ~isequal(1,aHS6)
        % bbp data is a 2D matrix and not a 1D linear array
        error('Only one row (reading) of HS6 data allowed per iteration');
    end

%% 2. Select appropriate particulate absorption value depth-bin for bbp spectrum
% Doxaran_SIGMA uses nearest neighbor approach to find the appropriate bin
depth_DIFF = abs(A_particulate_DeptH-HS6_deptH);
% Find the differences between the depth of the HS6 reading (a single
% scalar). This will tell us to which binned absorption depth the depth 
% of the single bbp spectrum (HS6_deptH) is closest.

depth_IND = find(depth_DIFF == min(depth_DIFF));
% Determins the index of the binned absorption specra that will be used 
aPART = A_particulate(depth_IND(1),:);
% Extracts the correct particulate absorption spectrum with which to
% calculate coefficient of attenuation (Kbb).

%% 3. Calculate Sigma for each bbp wavelength and multiply it by bbp spectrum

% Calculate Kbb spectrum via Doxaran et al. (2016) 4.34 is a coefficient
% derived by authors.
Kbb = aPART + 4.34*HS6_uncorrected; 
% Calculate Sigma using Kbb spectrum
sigmA = exp(kexp.*Kbb);
% Sigma-correct bbp by mutliplying Sigma and bbp spectrum
sigmaDOX_bb = HS6_uncorrected.*sigmA;
end