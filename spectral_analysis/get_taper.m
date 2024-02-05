function [taper,vals] = get_taper(n_times,type,opt)


%FUNCTION: 
%
%   Generates the requested taper type. Advanced users can modify this code
%   to incorporate custom taper functions. 
%
%INPUT:
%
%   n_times: Taper duration
%
%   type:   'hann'
%           'hamming'
%           'blackman_harris'
%           'dpss'
%
%   name-value pairs: These are used to specify multitaper parameters
%
%            NW: Time-bandwidth parameter. Increasing this value returns
%            more tapers with energy concentration >90%. Default 5/2 (4
%            tapers)
%
%            eigenvals: Set to 1 to return the eigenvalues associated with
%            each discrete prolate spheroidal sequence. Defaut 0
%
%OUTPUT:
%
%      taper: A periodic taper of the specified type
%
%       vals: Eigenvalues associated with the multitaper functions
%
%
%Dependencies: sig proc toolbox
%
% Notes: Common choices of NW are 2, 5/2, and 4. The program returns the
% first 2*NW-1 tapers, which concentrate most of the energy in the mainlobe
% for a given NW. A well-chosen NW reduces bias by attenuating the
% influence of the Fejér kernel, and reduces variance by smoothing over
% multiple eigenspectra.
%
%A. Nakhnikian 2024

arguments
    n_times
    type {mustBeMember(type,["hann","hamming","blackman_harris","dpss"])} 
    opt.NW = 5/2
    opt.eigenvals = 0
end

%check output configuration
vals = []; 
if opt.eigenvals && ~strcmp(type,'dpss')
    error('Eigenvalues are only defined for dpss')
end

switch type
    case 'hann'
        taper = hann(n_times,'periodic');
    case 'hamming'
        taper = hamming(n_times,'periodic');
    case 'blackman_harris'
        taper = blackmanharris(n_times,'periodic');
    case 'dpss'
        NW = opt.NW;
        if opt.eigenvals
            [taper,vals] = dpss(n_times,NW,[1,2*NW-1]);
        else
            [taper,~] = dpss(n_times,NW,[1,2*NW-1]);
        end
end

