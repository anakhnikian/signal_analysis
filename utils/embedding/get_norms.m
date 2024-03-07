function [norms] = get_norms(x,nrm)
%Computes L^nrm norm for each row of x. Set nrm to 'inf' for max norm

if isinf(nrm)
    norms = max(abs(x),[],2); %maximum norm
else
    norms = sum(abs(x).^nrm,2).^(1/nrm); %L^p norm
end