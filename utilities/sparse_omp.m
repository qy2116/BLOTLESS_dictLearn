function X_hat = sparse_omp(Y,D,k)
[~,n] = size(Y);
[~,l] = size(D);
X_hat = zeros(l,n);
for ind = 1:n
    X_hat(:,ind) = omp(D'*Y(:,ind), D'*D, k);
end
end