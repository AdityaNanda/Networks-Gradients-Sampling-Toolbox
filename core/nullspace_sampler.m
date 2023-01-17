function [x1, varargout] = nullspace_sampler(A, b, sig1, Z, samp)
% NULLSPACE_SAMPLER: general variant of nullspace sampling
%
%   [x1, Z] = nullspace_sampler(A, b, sig1, Z, samp)
%
% Inputs:
%    A:     constraint matrix (m x t)
%
%    b:     constraint vector (m x 1)
%
%    sig1:  constrains the std (positive) or norm (negative)
%           sig1 =  1  =>   std(x1) = 1
%           sig1 = -1  =>  norm(x1) = 1
%
%    Z:     nullspace for matrix A.  (t x [t-m])
%           if Z is empty, compute it with "null.m" (can be slow)
%
%    samp   number of samples (default is 1)
%
% Outputs:
%    x1:    sampled vector that satisfies A * x1 = b
%
%    Z1:    nullspace of an updated constraint matrix [A; x1']

t = size(A,2);
idx = find(isnan(b));   % ignore nans
b(idx) = [];
A(idx,:) = [];

% compute minimum norm solution
xmn = lsqminnorm(A,b);

% if Z is not supplied, compute Z
if ~exist('Z', 'var') || isempty(Z)
    Z=null(A);
end

if ~exist('samp', 'var') || isempty(samp)
    samp=1;
end

% compute d to satisfy std or norm constraints
mean_xn = mean(xmn);
if sig1 >= 0 % constrain std
    d = sqrt((t-1)*sig1.^2 + (t)*mean_xn^2 - norm(xmn).^2);
elseif sig1<0  % constrain the norm
    norm1 = -sig1;
    d = sqrt(norm1^2- norm(xmn)^2);
end

if ~isreal(d)
    warning('std or norm constraints not satisfied. Proceed with caution.');
    x1 = reshape(1,xmn,length(xmn));
    return;
end

mz = length(Z(1,:));                    % dimension of q
Nvar = normrnd(0, 1, [mz,samp]);        % sample normal dist
q = Nvar ./ (vecnorm(Nvar,2,1));        % uniformly sample q

x1 = reshape(xmn+d.*Z*q, [], samp);     % x1 = xmn +dZq

% update nullspace for next step
if nargout==2 && size(Z,2) > 1
    varargout{1} = null_expander(Z,q);
end
