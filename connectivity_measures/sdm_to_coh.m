function [coh] = sdm_to_coh(s,rm_diag)

%FUNCTION: Converts complex coherency matrices to coherence spectra
%
%INPUT:
%
%   s: Complex-valued, coherency matrices in channels by channels by
%   frequencies array
%
%   rm_diag: Binary flag for whether to set diagonal to zero. This is
%   necessary if using coherence matrices as network weights where
%   self-connections are exluced.
%
%
%OUTPUT:
%
%   coh: Coherence matrices in channels by channels by
%   frequencies arrays  
%
%A. Nakhnikian 2024

arguments
    s
    rm_diag = 0;
end
n_chans = size(s,1);
n_freqs = size(s,3);

coh = zeros(n_chans,n_chans,n_freqs);
for freq_ind = 1:n_freqs
    cross_spectra = squeeze(s(:,:,freq_ind));
    auto_spectra = diag(cross_spectra)*diag(cross_spectra)';
    coh(:,:,freq_ind) = abs(cross_spectra).^2./(auto_spectra);
end

if rm_diag
    coh = coh-eye(n_chans);
end
