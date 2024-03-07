function [lag_coh, lag_lin_dep] = sdm_to_lag_coh_uvar(s)

%FUNCTION: Computes pair-wise lagged coherence among all channels in the
%spectral matrix s. Lagged coherence eliminates zero-lag effects but does
%not account for multivariate interactions, ie shared input from a third
%node to a coherent pair is not regressed out.
%
%INPUT:
%
%       s is a spectral density matrix returned by
%       spectrum_to_sdm, channels by channels by frequencies
%
%
%OUTPUT:
%
%       lag_coh: univariate lagged coherence
%
%       lag_lin_dep: log scaled lagged linear dependence
%
%
%REF: https://arxiv.org/ftp/arxiv/papers/0711/0711.1455.pdf
%
%A. Nakhnikian 2024
n_chans = size(s,1);
n_freqs = size(s,3);
lag_coh = zeros(n_chans, n_chans, n_freqs);
for chan_x = 1:n_chans
    %generate full matrix, shortening the loop isn't worth the time suck to
    %add back the lower half Hermitian for trial sizes around 2000 or more
    for chan_y = 1:n_chans
        cross_spectra = squeeze(s(chan_x,chan_y,:));
        auto_spectra_x = squeeze(s(chan_x,chan_x,:));
        auto_spectra_y = squeeze(s(chan_y,chan_y,:));
        lag_coh(chan_x,chan_y,:) = (imag(cross_spectra).^2)./(auto_spectra_x.*auto_spectra_y-real(cross_spectra).^2);
    end
end
lag_lin_dep = -log(1-lag_coh);