function v1 = network_sampler(v0, m0, str)
% NETWORK_SAMPLER: sampling timeseries with network constraints
%
%   v1 = network_sampler(v0, m0, str)
%
% Inputs:
%    v0:    regional timeseries matrix (n x t)
%    m0:    network affiliation vector (n x 1)
%    str:   type of model: 'intra' or 'all' (default)
%               'intra': within-network correlations
%               'all': within and between-network correlations (default)
%
% Outputs:
%    v1:    sampled timeseries

% add required paths
addpath(genpath(fileparts(mfilename('fullpath'))));

[n, t] = size(v0);
k = max(m0);

% check inputs
if ~exist('str', 'var') || isempty(str)
    str = 'all';
end
if (length(m0) ~= n) || (k ~= length(unique(m0)))
    error('Incorrect module affiliation vector.');
end

%% generate synthetic mean network timeseries

vn0 = zeros(k, t);
tstd = vn0;  % standard deviation
for jj = 1:k
    ix = (m0==jj);  % all nodes for module jj
    vn0(jj,:) = mean(v0(ix,:), 1);       % mean activity
    tstd(jj,:) = std(v0(ix,:), [], 1);   % std
end
cn0 = cov(vn0');  % compute network covariances

switch str
    case 'all'
        % constrain all-to-all network correlations
        [ve, de] = eig(cn0);    % eigendcomposition
        [de, ix] = sort(diag(de), 'descend');
        ve = ve(:, ix);
        vnt = diag_cov_sampler(de,t);
        vn1 = ve * vnt;         % rotate back to original space

    case 'intra'
        % constrain network variances only
        A = ones(1,t);
        b = 0;
        for jj = k:-1:1
            sig1 = sqrt(cn0(jj, jj));
            vn1(jj,:) = nullspace_sampler(A, b, sig1);
        end
end

% sample regional activity for each network
v1 = zeros(size(v0));
for jj = 1:k
    ix = m0==jj; % all nodes for module jj
    A = ones(1,nnz(ix));
    b = 0;
    sig1 = 1;
    vtmp = nullspace_sampler(A, b, sig1, [], t);
    v1(ix,:) = (vtmp .* tstd(jj,:)) + vn1(jj,:);
end
