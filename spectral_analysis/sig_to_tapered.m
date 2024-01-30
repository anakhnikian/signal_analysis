function [sig_tapered,bw_factor] = sig_to_tapered(sig,taper,fs)

%FUNCTION: Subroutine called by functions that require tapered time series.
%Handles all 3 cases (single, multi, untapered).

% INPUT:
%
%       sig is a times by trials by channels array
%
%       taper is a vector of tapering coefficients matching the number of
%       rows in sig
%
% OUTPUT:
%
%       sig_tapered is an array of tapered times series. An extra dimension
%       is added for each Slepian sequence in the multitaper case
%
%       bw_factor is the equivalent noise bandwidth factor for converting
%       spectra (unit^2/Hz^2) to spectral densities (units^2/Hz).
%
% Notes: Untapered spectra are useful for demonstrating bias/variance
% tradeoffs of various estimators. I've included them in this function as
% a placeholder so that the enbw scaling is always returned by the same
% function and I don't have to code special cases for periodograms


if size(taper,2) == 1 %single taper
    sig_tapered = bsxfun(@times,sig,taper); 
    bw_factor   = fs.*sum(taper.^2); %SS*fs for standard taper

elseif size(taper,2)>2 %multitaper, adds tapers to 4th dim for vector 
    %time series and a 3rd for scalar (1 channel)
    n_tapers        = size(taper,2);
    sig_rep         = repmat(sig,1,1,1,n_tapers); %adds a dummy singleton in dim 2 for scalar data
    taper_prods     = (bsxfun(@times,permute(sig_rep,[1,4,2,3]),taper)); %match dims and multiply tapers
    tps_repermuted  = permute(taper_prods, [1,3,4,2]); %undo dim matching
    sig_tapered     = squeeze(tps_repermuted); %remove dummy singletons if necessary
    bw_factor       = fs; %dpss is unit SS normalized, bwf = 1*fs

elseif isempty(taper) %untapered
    sig_tapered = sig;
    bw_factor   = fs*size(sig,1); %SS = duration for all ones window

end