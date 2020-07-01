function I = isappdiag(A)
[m,n] = size(A);
for ind1 = 1:m
    for ind2 = 1:n
        A(ind1,ind2) = max(abs(A(ind1,ind2))-1e-4,0);
    end
end
I = isdiag(A);
end