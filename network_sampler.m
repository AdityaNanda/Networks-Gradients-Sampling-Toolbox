function v1 = network_sampler(v0, m0, str)

%% Inputs
% v0 is an (n x t) time-series for n-regions
% m0 is an (n x 1) module affiliation vector
%   must have values between 1 and ns (ns is the network size)
% str is the type of model: 'intra' or 'all' (default)
%   'intra' preserves only within-network correlations
%   'all' (default) preserves within and between-network correlations

% add required paths 
[folder, ~,~]= fileparts(mfilename('fullpath')); 
addpath(genpath(folder));

[n, t] = size(v0);
assert(length(m0)==n);
v1 = zeros(size(v0));

% check inputs
if ~exist('str', 'var') || isempty(str)
    str = 'all';
end
if ~exist('m0', 'var') || isempty(m0) || length(m0) ~=n
    error('Incorrect dimensionality of module affiliation vector.');
end

% vn0 is synthetic mean network timeseries
vn0 = zeros(max(m0), t);
tstd = vn0;  % tstd is the std. dev

for jj = 1:max(m0)
    ix = (m0==jj);  % all nodes for module jj

    vn0(jj,:) = mean(v0(ix,:), 1);       % mean activity
    tstd(jj,:) = std(v0(ix,:), [], 1);   % std
end
cn0 = cov(vn0');  % compute network covariances 

switch str
    case 'all'   % constrain all-to-all network connectivityies
        [ve, de] = eig(cn0); % eigendcomposition
        [de, ix] = sort(diag(de), 'descend');
        ve = ve(:, ix);
        vnt = diag_cov_sampler(de,t);
        vn1 = ve*vnt; % rotate back to original space

    case 'intra'  % constrain network variances only
        A = ones(1,t);
        b = 0;
        for jj = 1:max(m0)
            sig1 = sqrt(cn0(jj,jj));  % std= sqrt(variance)
            vn1(jj,:) = nullspace_sampler(A, b, sig1); 
        end
end

% sample regional activity for each network
for jj = 1:max(m0)
    ix = m0==jj; % all nodes for module jj
    A = ones(1,nnz(ix));
    b = 0;
    sig1 = 1;
    vtmp = nullspace_sampler(A, b, sig1, [], t); % t is the number of samples
    v1(ix,:)= (vtmp.*tstd(jj,:)) + vn1(jj,:);
end

end
