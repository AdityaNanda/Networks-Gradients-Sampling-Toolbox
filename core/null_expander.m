function Z = null_expander(Z1, q)
% Internal: Efficient computation of nullspace matrix
% For an initial nullspace Z1 of some A1, compute
% nullspace Z for an updated A1 = [A; xm + Z1 * q];

u = [q(1)-1; q(2:end)];   % see Householder's theorem
Z = Z1(:,2:end)- 2/norm(u)^2 .* (Z1*u)*u(2:end)';
