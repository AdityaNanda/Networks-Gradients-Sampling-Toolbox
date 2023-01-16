function v1 = diag_cov_sampler(d_cov, t)

% This function samples n time-series
% with a diagonal covairance matrix D. 
% d_cov is the diagonal variance values of D

% INPUTS
% d_cov is vector of variances to be satisfied
% t is the length of time series to be sampled

[n, ~] = size(d_cov);
assert(all(d_cov >=0)); % check eigenvalues are positive or 0

v1 = zeros(n,t);
Z1 = [];
for ii = 1:n
    % take all rows until ii-1 th step and modes
    A = [v1(1:ii-1,:); ones(1,t)];
    b = zeros(ii,1);    % define b

    sig1 = sqrt(d_cov(ii));  % std = sqrt(variance)

    [v1(ii, :), Z1] = nullspace_sampler(A, b, sig1, Z1);
end

end
