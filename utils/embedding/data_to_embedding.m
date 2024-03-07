function [lag_vects] = data_to_embedding(data, d, tau)

%FUNCTION: Returns a vector valued-time series given a scalar input and
%embedding parameters
%
%INPUT: 
%
%   data: scalar data times by trials by channels
%
%   d: embedding dimension
%
%   tau: embedding lag
%
%OUTPUT:
%
%   lag_vects: a d by times by trials by channels array where each column
%   is the reconstructed vector-valued series for a given
%   time/trial/channel combination

n_trials = size(data,2);
n_channels = size(data,3);
lag_data = data(1:tau:end,:,:); %time indices of interest
n_vects = floor(size(lag_data,1)/d);
embedding_data = lag_data(1:(d*n_vects),:,:); %truncate data to specified embedding
lag_vects = reshape(embedding_data,d, n_vects,n_channels,n_trials);
lag_vects = permute(lag_vects,[1,2,4,3]); 