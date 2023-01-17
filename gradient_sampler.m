function v1 = gradient_sampler(v0, k)
% GRADIENT_SAMPLER: sampling timeseries with gradient constraints
%
%   v1 = gradient_sampler(v0, k)
%
% Inputs:
%    v0:    regional timeseries matrix (n x t)
%    k:     number of principal components to be constrained
%
% Outputs:
%    v1:    sampled timeseries

% add required paths
addpath(genpath(fileparts(mfilename('fullpath'))));

[n, t] = size(v0);

c0 = cov(v0');          % n x n covariance matrix
[ve, de] = eig(c0);     % eigendecomposition of c0

% sort eigenvalues and eigenvectors
[de, ix] = sort(diag(de), 'descend');
ve = ve(:,ix);

% only sample eigenvalues that are large enough
ix = de > de(1) ./ 1e8;
de = de(ix);

% sample new eigenvectors ev1
pca0 = ve(:,1:k);
ev1 = pca_sampler(pca0, length(de));

% sample new timeseries with cov/corr c1
v2 = diag_cov_sampler(de, t);

v1 = ev1*v2;
