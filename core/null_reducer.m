function Z1 = null_reducer(Z, A, ii)
% Internal: Efficient computation of nullspace matrix
% For an initial nullspace Z of some A, compute
% nullspace Z1 for an updated A1 = A([1:ii-1 ii+1:end], :);

assert(size(Z,1) == size(A,2));
if size(A,2) - size(Z,2) < size(A,1)  % not full rank
    Z1 = Z;  % removing one row has no effect
else
    xi = A(ii,:);
    n = size(A,1);
    ix = setdiff(1:n, ii);
    xi = xi - (xi / A(ix,:)) * A(ix,:);  % regress
    Z1 = [Z, xi' ./ norm(xi)];
end
