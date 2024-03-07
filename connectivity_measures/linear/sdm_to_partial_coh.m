function [pcoh] = sdm_to_partial_coh(s,rm_diag)

%FUNCTION: Converts a cross-spectral density matrix to partial coherence
%spectra. Each spectral matrix is inverted and the entries are entered into
%the standard coherence function (abs(s_xy)^-2/(s_xx^-1*s_yy^-1)). This is
%equivalent to computing partial correlation from the precision
%matrix in the time domain.
%
%INPUT:
%
%   s: Complex-valued coherency matrices in a channels by channels by
%   frequencies array
%
%   rm_diag: Binary flag for whether to set diagonal to zero. This is
%   necessary if using coherence matrices as network weights where
%   self-connections are exluced. Default is zero (standard matrix where
%   self-coherence is all ones).
%
%
%OUTPUT:
%
%   pcoh: Partial coherence matrices in channels by channels by
%   frequencies arrays. Each entry is the pair-wise coherence with the
%   coherence among all other channels regressed out. 
%
%A. Nakhnikian 2024


arguments
    s
    rm_diag = 0;
end

n_chans = size(s,1);
n_freqs = size(s,3);

pcoh = zeros(n_chans,n_chans,n_freqs);
for freq_ind = 1:n_freqs
        s_xy_inv = inv(squeeze(s(:,:,freq_ind)));
        s_xx_inv = diag(s_xy_inv)*diag(s_xy_inv)';
        pcoh(:,:,freq_ind) = abs(s_xy_inv).^2./(s_xx_inv);       
end

if rm_diag
    pcoh = pcoh-eye(n_chans);
end
