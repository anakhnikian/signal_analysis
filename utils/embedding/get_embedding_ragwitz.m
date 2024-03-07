function [opt_d, opt_tau, mses, mse_surf] = get_embedding_ragwitz(data, k, nrm, d_range, tau_range)

%FUNCTION: searches for the optimal Takens embedding parameters for a
%collection of scalar valued time series by jointing minimizing the MS
%prediction error over all channels for the given test values of embedding
%dimension and delay

z_data = zscore(data);
n_delays = length(d_range);
n_lags = length(tau_range);
n_chans = size(data,3);
mses = zeros(n_delays, n_lags, n_chans);

for d_ind = 1:n_delays
    for tau_ind = 1:n_lags
        [lag_data] = data_to_embedding(z_data,d_range(d_ind),tau_range(tau_ind));
        parfor chan_ind = 1:n_chans
            chan_data = squeeze(lag_data(:,:,:,chan_ind));
            mses(d_ind,tau_ind,chan_ind) = get_prediction_error(chan_data, k, nrm);
        end
    end
end

%get maximum dim and lag over channels
mse_surf = mean(mses,3);
min_ind = find(mse_surf(:) == min(mse_surf(:)));
[r,c] = ind2sub(size(mse_surf),min_ind);
opt_d = d_range(r);
opt_tau = tau_range(c);
