function ev1= pca_sampler(ev0, p)

% This script takes as input ev0
% ev0 -> (n x k) matrix of k eigenvectors
% it returns ev1 an n x p orthonormal eigenvector basis


[n,k]= size(ev0);
assert(n>k);
ev1= zeros(n,p);   % new eigenvec basis
ev1(:, 1:k)= ev0;   % first few eigenvecs are simply ev0
z=[];  % initialize null matrix

for ii=k+1:p  % sample eigenvectors from k+1 to p

    A= ev1(:, 1:ii-1)';  % all eigenectors from 1 to ii-1

    if isempty(z)  % if null is empty compute
        z= null(A);
    else   % if null from previous iteration exists, update null matrix
        z=null_expander(z,q);
    end
    mz = length(z(1,:));                % dimension of q

    Nvar = normrnd(0,1,[mz,1]);      % sample normal dist
    q = Nvar ./ (vecnorm(Nvar,2,1));    % uniformly sample q
    ev1(:,ii)= z*q;  % min-norm solution is 0
end


