function [auto_spect_mat] = get_block_matrix(s_xx,s_yy)

%FUNCTION: Generates block diagonal matrices for multivariate measures
%
%INPUT: 
%       s_xx and s_yy: spectral density arrays in channels by channels by
%       frequencies arrangement
%
%OUTPUT:
%       An array of block diagonal matrices for each frequency 
%       M  = [s_xx O; O' s_yy]
%
%A. Nakhnikian 2024

%channel number can differ but frequency binning must be identical
n_chans_x = size(s_xx,1);
n_chans_y = size(s_yy,1);
n_freqs   = size(s_xx,3);
%get block auto-spectral matrices, blkdiag doesn't work by ND arrays so
%have to make the zero blocks by hand
off_diag_block_low  = zeros(n_chans_x,n_chans_y, n_freqs);
off_diag_block_high = zeros(n_chans_y,n_chans_x, n_freqs);
auto_spect_mat      = [[s_yy;off_diag_block_low],[off_diag_block_high; s_xx]];
