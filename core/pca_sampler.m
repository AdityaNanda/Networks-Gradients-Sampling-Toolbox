function ev1 = pca_sampler(ev0, p)
% PCA SAMPLER: sample an orthonormal eigenvector basis
%
%   ev1 = pca_sampler(ev0, p)
%
% Inputs:
%   ev0:    a matrix of orthonormal eigenvectors (n x k)
%     p:    number of eigenvectors to be sampled
%
% Outputs:
%   ev1:    orthonormal eigenvector basis (n x k)

[n,k] = size(ev0);
assert(n>k);
ev1 = zeros(n,p);           % new eigenvec basis
ev1(:, 1:k) = ev0;          % first few eigenvecs are simply ev0
z = [];                     % initialize null matrix

for ii = k+1:p              % sample eigenvectors from k+1 to p
    A = ev1(:, 1:ii-1)';    % all eigenectors from 1 to ii-1

    if isempty(z)           % if null is empty compute
        z = null(A);
    else                % update null matrix if previous iteration exists
        z = null_expander(z,q);
    end
    mz = length(z(1,:));    % dimension of q

    Nvar = normrnd(0, 1, [mz,1]);       % sample normal dist
    q = Nvar ./ (vecnorm(Nvar,2,1));    % uniformly sample q
    ev1(:,ii) = z*q;        % min-norm solution is 0
end
