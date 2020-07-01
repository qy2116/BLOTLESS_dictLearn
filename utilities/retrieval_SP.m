function SP = retrieval_SP(X)
[l,n] = size(X);
SP = cell(1,n);
for ind = 1:n
    x = X(:,ind);
    x = thr(x);
    pos = find(x~= 0);
    k = length(pos);
    temp_mat = zeros(l,k);
    for ind2 = 1:k
        temp_mat(pos(ind2),ind2) = 1;
    end
    SP{ind} = temp_mat;
end
end