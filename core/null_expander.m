function Z = null_expander(Z1, q)

% Given an initial nullspace Z1 for some A1 m times t matrix

% We return Z, the updated nullspace after you add 
% xm + Z1* q to the rows of A1 
% Thus, you do not have to use SVD to recompute the new nullspace 

 u = [q(1)-1;q(2:end)];   % see Householder's theorem
 Z = Z1(:,2:end)- 2/norm(u)^2 .* (Z1*u)*u(2:end)';