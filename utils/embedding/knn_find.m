function [knn, neighbor_inds, R] = knn_find(point, data, k, nrm)

norms = get_norms(point-data,nrm);
norms(norms == 0) = []; %remove self-distance if point in data

[sort_norms, sort_inds] = sort(norms);

neighbor_inds = sort_inds(1:k);
knn = data(neighbor_inds,:);
R = sort_norms(k);
