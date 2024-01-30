function [sdm] = spectrum_to_sdm(spectra, bw_factor)

%FUNCTION: Converts single trial coherency into spectral density matrices
%
%INPUT:
%
%   spectra: Channels by frequenies by trials by tapers
%            array.
%
%   bw_factor: Equivalent-noise bandwidth scaling factor for converting
%              spectra into spectral densities.
%
%OUTPUT:
%
%   sdm: A channels by channels by frequecies array of spectral density
%   matrices to be passed to connectivity estimation functions. Each matrix
%   in the rows and cols for the i_th page contains spectral densities for
%   the i_th frequency.
%
%
%
%A. Nakhnikian 2024

%Dev notes: Vectorized implementation slows down more quickly with size
%than nested loops. Increased cost for small data is around 1 second. That
%will add up if an analysis requires repeating the computation over many
%small data sets. In that case, adding a vectorize option for
%small data would be worth it

n_chans = size(spectra,1);
n_freqs   = size(spectra,2);
n_trials  = size(spectra,3);
n_tapers  = size(spectra,4);

%nested for loops are slightly slower than vectorizing for small data but
%slow down less for larger sets (eg 500 versus 2000 time points). 
trial_array = zeros(n_chans,n_chans,n_freqs,n_trials,n_tapers);
for chan1 = 1:n_chans
    %generate full matrix, shortening the loop isn't worth the time suck to
    %add back the lower half Hermitian for trial sizes around 2000 or more
    for chan2 = 1:n_chans 
        x = squeeze(spectra(chan1,:,:,:));
        y = squeeze(spectra(chan2,:,:,:));
        trial_array(chan1,chan2,:,:,:) = x.*conj(y);
    end
end

%average over relevant dimensions, squeeze out book keeping singletons
%left over for scalar data and/or single tapers
sdm  = squeeze(mean(trial_array./bw_factor,[4,5]));
sdm(:,:,2:end) = 2*sdm(:,:,2:end); %double everything except DC for two-sided spectra



