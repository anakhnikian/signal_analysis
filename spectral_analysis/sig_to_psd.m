function [psd,frq_ax] = sig_to_psd(sig,fs,taper)

%FUNCTION: Returns the power spectral density of a single time series witout
%computing full spectral matrix. Quick and easy program for assessing
%individual channels on the fly
%
%INPUT:
%
%       sig: A times by trials matrix
%
%       fs: sampling rate in Hz
%
%       taper: a taper returned by get_taper or other function
%
%OUTPUT:
%
%       psd: a one-side power spectrum
%
%       frq_ax: The frequency axis in Hz
%
%A. Nakhnikian 2024

[s, bw_factor, frq_ax] = sig_to_spectrum(sig,fs,taper); 
energy = s.*conj(s);
psd = mean(energy,[2,3])./bw_factor;
psd(2:end) = 2*psd(2:end); %double everything except last frequency