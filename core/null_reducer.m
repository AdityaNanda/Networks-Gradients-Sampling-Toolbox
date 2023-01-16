function Z1 = null_reducer(Z, A,ii)

% this function returns the nullspace of all rows in A except ii
% Z is the nullspace of A

assert(size(Z,1)==size(A,2));

if size(A,2)-size(Z,2) < size(A,1)  % not full rank
    Z1 = Z;  % removing one row has no effect

    return;
else
    xi = A(ii,:);
    n = size(A,1);
    ix = setdiff(1:n, ii);
    xi = xi- (xi/A(ix,:))* A(ix,:);  % regress
    Z1 = [Z, xi'./norm(xi)];
end

end
