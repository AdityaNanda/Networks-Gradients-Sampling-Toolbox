function c0 = net_corr(v,m)

% This script returns the network connectivity c0
% corresponding to time series (n x t) v and
% network affiliation vector m (n x 1)

% if v is zscored along each row,
% then c0_ij is equal to the
% average pairwise correlation between
% all nodes in the networks i and j

% mmax = max(m);  % total number of networks in m
% c0  is the correlation between the mean network timeseries 
% size(c0)= max(m) x max(m)

t= size(v,2);
 vtmp= zeros(max(m),t); % mean activity for each network
for jj = max(m):-1:1
    % mean activity of jj network
    vtmp(jj,:) = mean(v(m==jj,:),1);
end
c0= corr(vtmp');
