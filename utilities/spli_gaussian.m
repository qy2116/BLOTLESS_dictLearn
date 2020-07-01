function M = spli_gaussian(m,n,t)
M0 = randn(m,n);
M = zeros(m,n);
for ind1 = 1:m
    for ind2 = 1:n
        if M0(ind1,ind2)>0
            M(ind1,ind2) = M0(ind1,ind2) + t;
        else
            M(ind1,ind2) = M0(ind1,ind2) - t;
        end
    end
end
end