function [data] = multi_mode_ar2(fs,tr_len,tr_num,network_params)

% Generates an n_nodes*n_trials by frequencies array to mimic the
% connectivity data used our network study. Uses AR2
% to generate the data.
%
% fs: sampling rate in Hz
%
% tr_len: trial length in sec
%
% tr_num: total number of trials
%
% network_params has fields: couplings - cell with entry for each node pair set
%                        n_nodes
%                        freqs
%
% example:
%
% two network model with 10 nodes
%
% network_params.n_nodes = 10;
%
% Specify autonomous network dynamics, these are backgrounds activity over
% which coupled dynamics are superimposed
%
% Set Peak Frequencies
% network_params.AR.freq(1) = 10;
% network_params.AR.freq(2) = 10;
%
% Set Filter Bandwidth
% network_params.AR.modulus(1) = 0.9
% network_params.AR.modulus(2)  = 0.1
%
% Set Maximum Power Values (Default is unit power for all series)
% network_params.AR.max_power(1) = 1;
% network_params.AR.max_power(2) = 1.5;
%
% Set Couplings
%
% Set pairs matrices, row 1 is senders, row 2 is receivers, row 3 is
% coupling strength, 4th row is delay
% Set to empty for uncoupled network
% network_params.couplings{1}.pairs = [1, 2, 0.1, 0.2; 4, 10, 0.1, 0.2];
% network_params.couplings{2}.pairs = [5, 7, 0.2, 0.5; 1, 9, 0.1, 0.1];
%
% Set frequency parameters for senders
% network_params.drivers.freq(1) = 21;
% network_params.drivers.freq(2) = 27;
% network_params.drivers.modulus(1) = 0.5;
% network_params.drivers.modulus(2) = 0.9;
%
% Set transmission delays in seconds row 1 is forward row 2 is reverse
% network_params.couplings{1}.delays = [0.1,0.1];
% network_params.couplings{2}.delays = [0.2,0.2];
%
%
% Ref:
% Box, G. Jenkins, G. Reinsel, G. Ljung, G. (2016) Time Series
% Analysis: Forecasting and Control 5th Edition. Wiley and Sons, Hoboken NJ

%Preliminaries and input parsing
tr_len_samp = fs*tr_len; %trial length in samples
dt          = 1/fs; %sampling interval
n_nodes     = network_params.n_nodes;
n_nets      = length(network_params.couplings);
cent_freqs  = 2*pi*network_params.AR.freq./fs; %convert to radians per sample
if isfield(network_params,'drivers')
    drive_freqs = 2*pi*network_params.drivers.freq./fs;
end
moduli      = network_params.AR.modulus;
if n_nets ~= length(cent_freqs)
    error('Number of central frequencies must match number of subnets')
end
%Set flag for user specified maximum power levels
if ~isfield(network_params.AR,'max_power')
    max_default = 1;
else
    max_default = 0;
end

%Build coupling matrices
strength_mats = zeros(n_nodes,n_nodes,n_nets);
delay_mats = zeros(n_nodes,n_nodes,n_nets);
for net_ind = 1:n_nets
    coup = network_params.couplings{net_ind};
    if ~isempty(coup.pairs)
        %Get indices and convert to linear
        inds = coup.pairs(:,1:2);
        strengths = coup.pairs(:,3);
        delays = coup.pairs(:,4);
        n_pairs = size(inds,1);
        dim_3 = net_ind*ones(n_pairs,1); %which page (network) we're on
        lin_ind = sub2ind(size(strength_mats),inds(:,1),inds(:,2),dim_3);
        strength_mats(lin_ind) = strengths; %Assign ~0 couplings
        delay_mats(lin_ind) = round(delays/dt); %assign delays in time indices
    end
end
max_lag = max([2;delay_mats(:)]); %time horizon for all models

%Define the update matrices
ar_mats = zeros(n_nodes,n_nodes,max_lag,n_nets); %preallocate vector AR coefficient matrices
noise_level = zeros(n_nodes,n_nets);
for net_ind = 1:n_nets
    drivers = [];
    coup = network_params.couplings{net_ind};
    if ~isempty(coup.pairs)
        drivers = coup.pairs(:,1);
        receivers = coup.pairs(:,2);
    end

    for node_ind = 1:n_nodes
        
        %define AR params
        if any(drivers == node_ind) %define parameters for driver node
            mod = network_params.drivers.modulus(net_ind);
            freq = drive_freqs(net_ind);
        else %set parameters for background activity
            mod = moduli(net_ind);
            freq = cent_freqs(net_ind);
        end
        param2 = -mod^2; %Box and Jenkins 3.2.22
        param1 = cos(freq)*2*sqrt(-param2); %Box and Jenkins 3.2.23

        %Set noise level. The following is derived from the above equations
        %and the defintion of the PSD for an AR(2) process
        if max_default
            peak_pow = 1; %unit max power for all time series
        else
            peak_pow = network_params.AR.max_power(net_ind);
        end
        denom = ... %Denom of Box and Jenkins 3.2.29 with angular freq
            (1+param1^2+param2^2-2*param1*(1-param2)*cos(freq)-...
            2*param2*cos(2*freq));
        noise_level(node_ind,net_ind) = ...
            sqrt(0.5*fs*peak_pow*denom); %SD of the noise for one-sided spectrum
        %Place intrinsic dynamic params on diagonals
        ar_mats(node_ind,node_ind,1,net_ind) = param1;
        ar_mats(node_ind,node_ind,2,net_ind) = param2;

        %couplings are the rest of the matrices
        coup_mat = squeeze(strength_mats(:,:,net_ind)); %get couplings for subnet
        for pair_ind = 1:length(drivers) %loop over coupled pairs
            r = drivers(pair_ind);
            c = receivers(pair_ind);
            coup_strength = coup_mat(r,c);
            ar_mats(c,r,delay_mats(r,c,net_ind),net_ind) = coup_strength;
        end

    end
end

%Create the time series
padded_len =tr_len_samp+2*max_lag;
vect_ts = zeros(n_nodes,padded_len+max_lag,tr_num,n_nets); %add left padding to compensate for time horizon due to coupling
for net_ind = 1:n_nets
    subnet_mats = squeeze(ar_mats(:,:,:,net_ind)); %subnet parameters
    for tr_ind = 1:tr_num %DEV: vectorize if worth it
        tr_series = zeros(n_nodes,padded_len+max_lag); %preallocate vectors for trial
        tr_series(:,1:max_lag) = noise_level(:,net_ind).*randn(n_nodes,max_lag); %noise in time horizon
        for time_ind = max_lag+1:padded_len+max_lag
            %Get lagged vectors and add dummy dim for page-wise
            %multiplication
            lag_vects = tr_series(:,time_ind-1:-1:time_ind-max_lag);
            lag_vects = reshape(lag_vects,n_nodes,1,size(lag_vects,2));
            %weight lags, sum, add noise
            ar_products = pagemtimes(subnet_mats,lag_vects);
            ep = noise_level(:,net_ind).*randn(n_nodes,1);
            tr_series(:,time_ind) = sum(ar_products,3)+ep;
        end
        vect_ts(:,:,tr_ind,net_ind) = tr_series;
    end
end
%Remove padding and get final output
data = squeeze(mean(vect_ts(:,3*max_lag+1:end,:,:),4));
data = permute(data,[2,3,1]); %get times by trials by channels


