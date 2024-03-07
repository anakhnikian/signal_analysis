function [mse] = get_prediction_error(data, k, nrm)

%FUNCTION: Computes the MSE between the future state of a system and the
%mean of the future of its k nearest neighbors. Uses the locally constant
%predictor of Ragwitz and Kantz, where the predicted and obtained are the
%current state of the underlying scalar time series (the first entry in
%each delay vector).
%
%INPUT:
%
%   data: a dims by times by trials vector-valued time series
%
%   k: number of nearest neighbors to search for
%
%   nrm: norm to be used for nn search
%
%OUPUT: 
%
%   mse: the mean squared error between obtained and predicted future for
%   all points in the data set
%
%REF: 
%
%   Ragwitz and Kantz (2002) "Markov models from data  by simple nonlinear
%   time series predictors in delay embedding spaces" Phys Rev E 65, 056201
%
% A. Nakhnikian 2024


n_dims = size(data,1);
n_times = size(data,2);
n_trials = size(data,3);

%collapse trials for nn search, remove t(end) from search
full_data = reshape(data(:,1:n_times-1,:),n_dims,(n_times-1)*n_trials); 
wrapped_times = repmat(1:n_times-1,n_trials,1); %for time indexing full data

err = zeros(n_times,1);
for time_ind = 1:n_times-1
    trial_data = squeeze(data(:,time_ind,:));

    for tr_ind = 1:n_trials
        point = trial_data(:,tr_ind);
        [~, nn_inds,~] = knn_find(point', full_data', k, nrm);

        %get time and trial indices for nn in full phase space
        time_nn = wrapped_times(nn_inds);
        trial_nn = ceil(nn_inds/n_times);

        %get nn and computed error from LC prediction
        %index the full array to avoid "advancing" time over trial gaps
        future_nn = squeeze(data(:, [time_nn+1,trial_nn]));
        predicted = mean(future_nn(1,:)); %locally constant pred
        obtained = data(1,time_ind+1,tr_ind);
        err(time_ind) = predicted-obtained;
    end
    
end

mse = mean(err.^2);