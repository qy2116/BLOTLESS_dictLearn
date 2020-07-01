function s = shrink(X,thr)
[m,n] = size(X);
s = zeros(m,n);
for ind1 = 1:m
    for ind2 = 1:n
    s(ind1,ind2) = sign(X(ind1,ind2)) * max(abs(X(ind1,ind2))-thr,0);
    end
end

end