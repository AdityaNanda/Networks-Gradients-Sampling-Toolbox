function c = net_corr(v, m)
% NET_CORR: compute average network correlations
%
%   c = net_corr(v, m)
%
% Inputs:
%    v:     regional timeseries matrix (n x t)
%    m:     network affiliation vector (n x 1)
%
% Outputs:
%    c:     network correlations

k = max(m);
c0 = corr(v');
c = zeros(k);

% mean network correlations
for ii = k:-1:1
    for jj = ii:-1:1
        c(ii, jj) = mean(c0(m==ii, m==jj), 'all');
        c(jj, ii) = c(ii, jj);
    end
end
