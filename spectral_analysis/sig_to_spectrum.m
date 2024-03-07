function [s_trials, bw_factor, frq_ax] = sig_to_spectrum(sig,fs,taper)

%FUNCTION: Returns complex-valued single-trial spectra in channels by
%frequencies by trials arrays. These are passed to other functions to
%convert the trial data to spectral densities.
%
%INPUT:
%
%   sig: A scalar or vector-valued signal. Dims are times by trials by channels.
%   Scalar data are in times by trials matrices.
%
%   fs: Sampling rate in Hz
%
%   taper: A single windowing function or collection of orthogonal
%   functions (multitaper case) returned by get_taper
%
%OUTPUT:
%
%   s_trials: Channels by n_frequencies by n_trials by n_tapers arrays where each
%   matrix in the rows and columns is a Hermitian matrix with raw complex coherency
%   values
%
%   bw_factor: The equivalent noise-bandwidth coefficient returns by
%   sig_to_tapered. Output here to be passed to spectral density estimator
%   functions
%
%   frq_ax: The frequency axis in Hz
%
%A. Nakhnikian 2024

[sig_tapered,bw_factor] = sig_to_tapered(sig,taper,fs); 
n_times                 = size(sig,1);
df                      = fs/n_times;

%Get index of last positive, non-Nyquist, frequency
if rem(n_times,2) %odd length, no Nyquist freq
    last_freq_index = (n_times-1)/2;
else %Nyquist defined
    last_freq_index = n_times/2;
end
frq_ax = 0:df:floor(fs/2)-df;

%Fourier transform, remove negative frequencies
sig_ft   = fft(sig_tapered);
s_trials = sig_ft(1:last_freq_index,:,:,:);


