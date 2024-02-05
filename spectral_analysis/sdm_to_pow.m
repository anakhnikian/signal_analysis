function [psd] = sdm_to_pow(s)

%FUNCTION: Extracts the power spectral densities along the diagonals of the
%spectral density matrix s
%
%INPUT: s is a channels by channels by frequencies array of spectral
%matrices
%
%OUTPUT:psd is a channels by frequencies matrix of one-sided power spectral
%densities.
%
%A. Nakhnikian, 2024

n_chans = size(s,1);
n_freqs = size(s,3);

psd = zeros(n_chans,n_freqs);
for freq_ind = 1:n_freqs
    psd(:,freq_ind) = diag(s(:,:,freq_ind));
end
