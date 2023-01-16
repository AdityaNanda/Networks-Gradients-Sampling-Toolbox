function v1 = gradient_sampler(v0, k, varargin)

%% Inputs
% v0 is a n x t time-series for n-regions
% k is the number of PCA components (gradients) to be constrained
% str is used to toggle between 'covariance' and 'product' modes
% str = 'cov' (default) computes components on the covariance matrix
% str = 'prod' uses the product matrix (v0' x v0)
% str = 'corr' uses the corr matrix (corr(v0'))

%% Intermediates
% c0 is the coavariance/corr/product matrix of v0
% pca0 is the set of constrained components
% all eigenvalues of c0 are the same as those of c1 (only eigenvectors are sampled)

% outputs
% c1 is a n x n  model covariance matrix
% v1 - n x t time-series with cov c1

[folder, ~,~]= fileparts(mfilename('fullpath'));
addpath(genpath(folder));

[n, t] = size(v0);

c0 = cov(v0'); % n x n covariance matrix
[ve, de] = eig(c0);  % eigen decomposition of c0
[de, ix] = sort(diag(de), 'descend');  % sort eigen vals
ve = ve(:,ix);

% only sample eigenvalues that are large enough
ix = de > de(1)./1e8;
de = de(ix);

%sample new eigenvectors ev1
pca0 = ve(:,1:k);
ev1 = pca_sampler(pca0, length(de));

% sample new timeseries with cov/corr c1
v2 = diag_cov_sampler(de, t);

v1= ev1*v2;

end