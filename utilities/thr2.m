function x = thr2(A)
[m,n] = size(A);
for ind1 = 1:m
    for ind2 = 1:n
    x(ind1,ind2) = sign(A(ind1,ind2)) * max(abs(A(ind1,ind2))-0.3,0);
    end
end    
end