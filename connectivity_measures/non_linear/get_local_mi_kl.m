function [mi] = get_local_mi_kl(data, k, nrm)

%FUNCTION: Returns the local (time-resolved) mutual information among a
%collection of state variables via the Kozachenko-Leonenko estimator using
%the generalization of Kraskov et al. to multivariate systems. MI is
%localized in time using the trial-wise search method of Gómez-Herrero et
%al. (2015). 
% 
%If size(data,4) = 2 then this is the simple pairwise mutual
%information between two channels. Additional channels may be added to
%generalize to I([x(1),x(2),...,x(N)]) as measures of global concordance
%among >2 processes. Note that the curse of dimensionality becomes a
%concern with increasing channels, especially if the optimal delay
%dimension is high.
%
%INPUT: 
%
%   data: A dimensions by times by trials by channels array of times-series
%         where the rows contain the elements of the delay vectors obtained via
%         Takens embedding
%
%   k:    The number of nearest neighbors to include in the knn searches
%
%   nrm:  The vector norm order to define distance (Lp norm)
%
%OUTPUT:
%
%   mi: mutual information
%
%REF: 
% 
%   A. Kraskov, H. Stögbauer, and P. Grassberger (2004) Physical review E,
%   69(6):066138,
%
%   Gómez-Herrero et al (2015) "Assessing Coupling Dynamics from an
%   Ensemble of Time Series" Entropy 17 (4) 1958-1970
%
%Notes: 
%
%       
%       
%       This program can handle scalar time series. Let A be a time
%       by trials by channels array of observations. data(1,:,:,:) = A will
%       return the local MI among the channels in A. 
% 
%       inf (maximum absolute value) is the recommended norm. Theoretically
%       any norm - or any metric - should be suitable for MI estimation. In
%       practice, I've found that the Euclidean norm 
%       underesimtates the known MI of a theoretically tractable system
%       (two unit variance Gaussian processes with 0<cov(x,y)<1). 

%data parameters
n_trials = size(data,3);
n_channels = size(data,4);

%knn search output
joint_radii = get_joint_radii(data, k, nrm);
nei_count = get_neighbor_counts(data, joint_radii, nrm);

%mutual information
psi_k = psi(k);
psi_n = (n_channels-1)*psi(n_trials); %local H is indexed over trials
psi_count = psi(nei_count+1); %include center point in marginal space->N+1
psi_sum = squeeze(sum(psi_count,3));
mi = psi_k+psi_n-mean(psi_sum,2); %Kraskov et al Eq23

