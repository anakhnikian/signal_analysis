function [psd] = sdm_to_pow(s)

n_chans = size(s,1);
n_freqs = size(s,3);

psd = zeros(n_chans,n_freqs);
for freq_ind = 1:n_freqs
    psd(:,freq_ind) = diag(s(:,:,freq_ind));
end